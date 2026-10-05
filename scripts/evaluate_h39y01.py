#!/usr/bin/env python3
"""Run the single locked H39Y-01 spatial holdout and SPRT.

The script intentionally has no candidate, feature, seed, grid, emitter, or
boundary sweep. ``--plan-only`` checks the preregistered geometry and truth-count
eligibility before the candidate field is computed. The normal run scores only
H39Y-01, stops at the first SPRT boundary, and creates a submission artifact
only after H1 is accepted.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import platform
import sys
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
import scipy
import sklearn
from scipy.ndimage import binary_dilation, distance_transform_edt, label

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import emission, grid, radiometric, stack  # noqa: E402
from gems39.io39 import audit as audit_submission, sha256 as sha256_file, write as write_submission  # noqa: E402
from gems39.metric import dti, dti_binary  # noqa: E402
from gems39.sprt_select import PairwiseSPRT  # noqa: E402

DATA = ROOT / "data"
HYPOTHESIS_ID = "H39Y-01"
GRID_N = 8
COLLAR_PX = 15
CELL_EROSION_PX = 15
MIN_TRUTH_PX = 50
MIN_ELIGIBLE_CELLS = 9
BASE_SEED = 20261005
MAX_TRAIN_ROWS_PER_CLASS = 90_000
MIN_TRAIN_ROWS_PER_CLASS = 200
EMISSION_BUDGET = 24_000
MIN_DISTANCE_PX = 2.8
KNOWN_BUFFER_PX = 2
SPRT_SETTINGS = dict(p0=0.50, p1=0.70, alpha=0.05, beta=0.10)
NEGATIVE_BACKGROUND_DILATION = 8
HARD_NEGATIVE_MIN_PX = 1.0
HARD_NEGATIVE_MAX_PX = 6.0


@dataclass
class Cell:
    ordinal: int
    row: int
    col: int
    y0: int
    y1: int
    x0: int
    x1: int
    active: np.ndarray
    truth: np.ndarray
    hidden_ids: np.ndarray

    @property
    def n_active(self) -> int:
        return int(self.active.sum())

    @property
    def n_truth(self) -> int:
        return int(self.truth.sum())

    @property
    def eligible(self) -> bool:
        return self.n_truth >= MIN_TRUTH_PX

    @property
    def bounds(self) -> tuple[slice, slice]:
        return slice(self.y0, self.y1), slice(self.x0, self.x1)

    @property
    def seed(self) -> int:
        return BASE_SEED + self.ordinal

    def as_dict(self) -> dict[str, Any]:
        return dict(
            ordinal=self.ordinal,
            row=self.row,
            col=self.col,
            bbox=[self.y0, self.y1, self.x0, self.x1],
            seed=self.seed,
            active_pixels=self.n_active,
            truth_pixels=self.n_truth,
            eligible=self.eligible,
            hidden_component_count=int(self.hidden_ids.size),
        )


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot JSON-encode {type(value).__name__}")


def _read_inputs(data_dir: Path) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    sample_path = data_dir / "sample_submission.tif"
    labels_path = data_dir / "labels.tif"
    features_path = data_dir / "training_features.tif"
    ratio_path = data_dir / "external" / "geodawn_extensions_u8.tif"
    required = [sample_path, labels_path, features_path, ratio_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required local inputs: " + ", ".join(missing))

    with rasterio.open(sample_path) as sample:
        sample_array = sample.read(1)
        footprint = np.isfinite(sample_array)
        reference = dict(
            height=sample.height,
            width=sample.width,
            crs=str(sample.crs),
            transform=tuple(sample.transform)[:6],
            bounds=tuple(sample.bounds),
        )
    if not footprint.any():
        raise ValueError("Competition footprint is empty")

    with rasterio.open(labels_path) as labels_src:
        if (
            labels_src.shape != footprint.shape
            or str(labels_src.crs) != reference["crs"]
            or tuple(labels_src.transform)[:6] != reference["transform"]
        ):
            raise ValueError("Catalogue labels do not match the sample submission grid")
        raw_labels = labels_src.read(1)
    outside_labels = (raw_labels >= 1) & ~footprint
    if outside_labels.any():
        raise ValueError(f"Found {int(outside_labels.sum())} catalogue labels outside the sample footprint")
    catalogue = (raw_labels >= 1) & footprint
    return footprint, catalogue, dict(
        sample=sample_path,
        labels=labels_path,
        features=features_path,
        ratios=ratio_path,
        reference=reference,
        raw_labels_outside_footprint=int(outside_labels.sum()),
        footprint_pixels=int(footprint.sum()),
        catalogue_pixels=int(catalogue.sum()),
    )


def _make_cells(footprint: np.ndarray, catalogue: np.ndarray) -> tuple[list[Cell], dict[str, Any], np.ndarray]:
    """Build the exact registered 8×8 bounding-box grid and even/even cells."""
    y_present = np.flatnonzero(footprint.any(axis=1))
    x_present = np.flatnonzero(footprint.any(axis=0))
    if y_present.size == 0 or x_present.size == 0:
        raise ValueError("No finite footprint rows or columns")
    y_edges = np.linspace(int(y_present[0]), int(y_present[-1]) + 1, GRID_N + 1).astype(int)
    x_edges = np.linspace(int(x_present[0]), int(x_present[-1]) + 1, GRID_N + 1).astype(int)
    if np.any(np.diff(y_edges) <= 0) or np.any(np.diff(x_edges) <= 0):
        raise ValueError("Footprint bounding box is too small for the locked 8×8 grid")

    components, component_count = label(catalogue, structure=np.ones((3, 3), dtype=np.uint8))
    cells: list[Cell] = []
    ordinal = 0
    for row in range(0, GRID_N, 2):
        for col in range(0, GRID_N, 2):
            y0, y1 = int(y_edges[row]), int(y_edges[row + 1])
            x0, x1 = int(x_edges[col]), int(x_edges[col + 1])
            sub_foot = footprint[y0:y1, x0:x1]
            # A false-valued 15-pixel pad makes the rectangular cell exterior
            # part of the Euclidean inward-erosion distance calculation.
            padded = np.pad(sub_foot, CELL_EROSION_PX, mode="constant", constant_values=False)
            distance = distance_transform_edt(padded)
            inner_distance = distance[
                CELL_EROSION_PX:-CELL_EROSION_PX,
                CELL_EROSION_PX:-CELL_EROSION_PX,
            ]
            active = sub_foot & (inner_distance > CELL_EROSION_PX)
            truth = catalogue[y0:y1, x0:x1] & active

            gy0 = max(0, y0 - COLLAR_PX)
            gy1 = min(footprint.shape[0], y1 + COLLAR_PX)
            gx0 = max(0, x0 - COLLAR_PX)
            gx1 = min(footprint.shape[1], x1 + COLLAR_PX)
            hidden_ids = np.unique(components[gy0:gy1, gx0:gx1])
            hidden_ids = hidden_ids[hidden_ids > 0].astype(np.int32, copy=False)
            cells.append(
                Cell(
                    ordinal=ordinal,
                    row=row,
                    col=col,
                    y0=y0,
                    y1=y1,
                    x0=x0,
                    x1=x1,
                    active=active,
                    truth=truth,
                    hidden_ids=hidden_ids,
                )
            )
            ordinal += 1
    grid_record = dict(
        grid_n=GRID_N,
        anchor="bounding rectangle of finite sample-submission footprint",
        row_edges=y_edges.tolist(),
        col_edges=x_edges.tolist(),
        tested_cells="row and column indices both even, row-major",
        erosion=dict(method="Euclidean distance transform", pixels=CELL_EROSION_PX,
                     active_rule="footprint cells at distance > 15 px from cell exterior or invalid footprint"),
        component_connectivity=8,
        collar_px=COLLAR_PX,
        full_catalogue_component_count=int(component_count),
    )
    return cells, grid_record, components


def _fold_masks(
    cell: Cell,
    footprint: np.ndarray,
    catalogue: np.ndarray,
    components: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return train mask, training-only catalogue, and hidden component mask."""
    hidden_lookup = np.zeros(int(components.max()) + 1, dtype=bool)
    hidden_lookup[cell.hidden_ids] = True
    hidden_map = hidden_lookup[components]
    train_catalogue = catalogue & ~hidden_map

    train_mask = footprint.copy()
    gy0 = max(0, cell.y0 - COLLAR_PX)
    gy1 = min(footprint.shape[0], cell.y1 + COLLAR_PX)
    gx0 = max(0, cell.x0 - COLLAR_PX)
    gx1 = min(footprint.shape[1], cell.x1 + COLLAR_PX)
    train_mask[gy0:gy1, gx0:gx1] = False
    train_mask &= ~hidden_map
    return train_mask, train_catalogue, hidden_map


def _training_pools(
    cell: Cell,
    footprint: np.ndarray,
    catalogue: np.ndarray,
    components: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    train_mask, train_catalogue, hidden_map = _fold_masks(cell, footprint, catalogue, components)
    distance = distance_transform_edt(~train_catalogue)
    hard = (
        (distance >= HARD_NEGATIVE_MIN_PX)
        & (distance <= HARD_NEGATIVE_MAX_PX)
        & footprint
        & ~train_catalogue
    )
    background = (
        footprint
        & ~train_catalogue
        & ~binary_dilation(train_catalogue, iterations=NEGATIVE_BACKGROUND_DILATION)
    )
    negative = hard | background
    n_positive = int((train_catalogue & train_mask).sum())
    n_negative = int((negative & train_mask).sum())
    counts = dict(
        train_positive=n_positive,
        train_negative=n_negative,
        train_hard_negative=int((hard & train_mask).sum()),
        train_background_negative=int((background & train_mask).sum()),
    )
    return train_mask, train_catalogue, negative, counts


def _check_source_band_order(features_path: Path, reference: dict[str, Any]) -> list[str]:
    with rasterio.open(features_path) as src:
        if (
            (src.height, src.width) != (reference["height"], reference["width"])
            or str(src.crs) != reference["crs"]
            or tuple(src.transform)[:6] != reference["transform"]
        ):
            raise ValueError("Competition feature stack does not match sample grid")
        by_band: list[str] = []
        for band in range(1, src.count + 1):
            description = src.descriptions[band - 1] or ""
            by_band.append(description.split(" - ", 1)[0].strip().casefold())
    expected = [name.casefold() for name in stack.BAND_NAMES]
    if by_band != expected:
        raise ValueError(
            "Feature-band description order differs from the checked-in H39-A stack: "
            f"actual={by_band}; expected={expected}"
        )
    return by_band


def _prepare_stack(features_path: Path, footprint: np.ndarray) -> tuple[np.ndarray, list[str]]:
    print("[3/6] Build the checked-in H39-A feature stack once", flush=True)
    fields = stack.build_stack(features_path, footprint, quantise=True)
    columns = sorted(fields)
    if not columns:
        raise ValueError("The H39-A feature stack is empty")
    cube = np.empty((*footprint.shape, len(columns)), dtype=np.uint8)
    for j, name in enumerate(columns):
        value = fields[name]
        if value.shape != footprint.shape or value.dtype != np.uint8:
            raise ValueError(f"Unexpected baseline feature array for {name!r}: {value.shape}/{value.dtype}")
        cube[..., j] = value
    del fields
    print(f"  features={len(columns)}; cube={cube.nbytes / (1024 ** 2):.1f} MiB", flush=True)
    return cube, columns


def _fit_predict_h39a(
    cube: np.ndarray,
    cell: Cell,
    train_mask: np.ndarray,
    train_catalogue: np.ndarray,
    negative_pool: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    from sklearn.ensemble import HistGradientBoostingClassifier

    seed = cell.seed
    positive_y, positive_x = np.nonzero(train_catalogue & train_mask)
    negative_y, negative_x = np.nonzero(negative_pool & train_mask)
    m = int(min(positive_y.size, negative_y.size, MAX_TRAIN_ROWS_PER_CLASS))
    if m < MIN_TRAIN_ROWS_PER_CLASS:
        raise ValueError(
            f"Fold {cell.ordinal} has only {m} balanced rows; preregistered minimum is "
            f"{MIN_TRAIN_ROWS_PER_CLASS} per class"
        )
    rng = np.random.default_rng(seed)
    positive_sample = rng.choice(positive_y.size, size=m, replace=False)
    negative_sample = rng.choice(negative_y.size, size=m, replace=False)
    x_train = np.empty((2 * m, cube.shape[2]), dtype=np.float32)
    x_train[:m] = cube[positive_y[positive_sample], positive_x[positive_sample]]
    x_train[m:] = cube[negative_y[negative_sample], negative_x[negative_sample]]
    x_train *= np.float32(1.0 / 255.0)
    y_train = np.concatenate((np.ones(m, dtype=np.uint8), np.zeros(m, dtype=np.uint8)))

    model = HistGradientBoostingClassifier(
        max_iter=250,
        learning_rate=0.08,
        max_leaf_nodes=31,
        l2_regularization=1.0,
        early_stopping=False,
        random_state=seed,
    )
    model.fit(x_train, y_train)
    del x_train, y_train, positive_sample, negative_sample

    probability = np.zeros(cell.active.shape, dtype=np.float32)
    active_y, active_x = np.nonzero(cell.active)
    step = 50_000
    for start in range(0, active_y.size, step):
        stop = min(active_y.size, start + step)
        ys = active_y[start:stop] + cell.y0
        xs = active_x[start:stop] + cell.x0
        features = cube[ys, xs].astype(np.float32)
        features *= np.float32(1.0 / 255.0)
        probability[active_y[start:stop], active_x[start:stop]] = model.predict_proba(features)[:, 1]
        del features
    if not np.isfinite(probability).all() or np.any((probability < 0.0) | (probability > 1.0)):
        raise ValueError("H39-A comparator emitted non-finite or out-of-range probabilities")
    model_record = dict(
        seed=seed,
        rows_per_class=m,
        max_iter=250,
        learning_rate=0.08,
        max_leaf_nodes=31,
        l2_regularization=1.0,
        early_stopping=False,
        probability_min=float(probability[cell.active].min()) if cell.n_active else 0.0,
        probability_max=float(probability[cell.active].max()) if cell.n_active else 0.0,
    )
    del model
    return probability, model_record


def _score_tile(
    field: np.ndarray,
    comparator: np.ndarray,
    cell: Cell,
    visible_catalogue: np.ndarray,
    total_footprint: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    yslice, xslice = cell.bounds
    local_known = visible_catalogue[yslice, xslice]
    target_n = int(round(EMISSION_BUDGET * cell.n_active / total_footprint))
    candidate_points = emission.emit_fast(
        field[yslice, xslice],
        cell.active,
        target_n=target_n,
        min_dist=MIN_DISTANCE_PX,
        catalogue=local_known,
        cat_buffer_px=KNOWN_BUFFER_PX,
    )
    comparator_points = emission.emit_fast(
        comparator,
        cell.active,
        target_n=target_n,
        min_dist=MIN_DISTANCE_PX,
        catalogue=local_known,
        cat_buffer_px=KNOWN_BUFFER_PX,
    )
    candidate_result = dti_binary(
        candidate_points,
        cell.truth,
        valid=cell.active,
        known=local_known,
    )
    comparator_result = dti_binary(
        comparator_points,
        cell.truth,
        valid=cell.active,
        known=local_known,
    )
    if candidate_result["n_truth"] != cell.n_truth or comparator_result["n_truth"] != cell.n_truth:
        raise ValueError("DTI truth count differs from the precomputed hidden-cell truth count")
    record = dict(
        ordinal=cell.ordinal,
        row=cell.row,
        col=cell.col,
        seed=cell.seed,
        active_pixels=cell.n_active,
        truth_pixels=cell.n_truth,
        per_tile_target_pixels=target_n,
        candidate_emitted_pixels=int(candidate_points.sum()),
        comparator_emitted_pixels=int(comparator_points.sum()),
        candidate_dti=candidate_result,
        comparator_dti=comparator_result,
        candidate_strict_win=bool(candidate_result["dti"] > comparator_result["dti"]),
    )
    return record, dict(candidate=candidate_points, comparator=comparator_points)


def _pooled_dti(records: list[dict[str, Any]], key: str) -> float | None:
    if not records:
        return None
    tp = sum(float(record[key]["tp"]) for record in records)
    fp = sum(float(record[key]["fp"]) for record in records)
    fn = sum(float(record[key]["fn"]) for record in records)
    return dti(tp, fp, fn)


def _write_progress(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=_json_default, allow_nan=False) + "\n")
    tmp.replace(path)


def _emit_submission(
    candidate_field: np.ndarray,
    footprint: np.ndarray,
    catalogue: np.ndarray,
    sample_path: Path,
    out_dir: Path,
) -> dict[str, Any]:
    print("[6/6] H1 accepted: generating a fresh full-footprint candidate TIFF", flush=True)
    mask = emission.emit_fast(
        candidate_field,
        footprint,
        target_n=EMISSION_BUDGET,
        min_dist=MIN_DISTANCE_PX,
        catalogue=catalogue,
        cat_buffer_px=KNOWN_BUFFER_PX,
    )
    if int(mask.sum()) != EMISSION_BUDGET:
        raise ValueError(f"Full-footprint emitter returned {int(mask.sum())}, expected {EMISSION_BUDGET}")
    packed_digest = hashlib.sha256(np.packbits(mask, bitorder="little").tobytes()).hexdigest()
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    name = f"gemsdoe39-h39y01-radiometric-cond-{stamp}-{packed_digest[:8]}"
    out_dir.mkdir(parents=True, exist_ok=True)
    tif_path = out_dir / f"{name}-zeros.tif"
    zip_path = out_dir / f"{name}-zeros.zip"
    manifest_path = out_dir / f"{name}-manifest.json"
    if any(path.exists() for path in (tif_path, zip_path, manifest_path)):
        raise FileExistsError(f"Refusing to overwrite candidate output beginning {name!r}")

    # The all-finite, no-nodata encoding directly avoids the reported portal
    # [0,1] failure modes; predictions are still exactly binary in-footprint.
    write_submission(mask.astype(np.float32), sample_path, tif_path, footprint, outside="zeros")
    audit = audit_submission(tif_path, footprint, sample_path)
    if not audit.get("ok", False):
        tif_path.unlink(missing_ok=True)
        raise ValueError(f"Independent GeoTIFF audit failed: {audit}")
    with rasterio.open(tif_path) as src:
        written = src.read(1)
        if not np.isfinite(written).all() or np.any((written < 0) | (written > 1)):
            raise ValueError("Candidate GeoTIFF does not have all-finite [0,1] values")
        if not np.array_equal(written[footprint], mask[footprint].astype(np.float32)):
            raise ValueError("Candidate GeoTIFF in-footprint pixels changed during writing")

    # Compare exact in-footprint prediction arrays against prior artifacts local
    # to this checkout. The broader sibling-repository raster corpus was not
    # restored in this run, so this is explicitly a limited-scope audit.
    prior_paths = sorted((ROOT / "docs" / "downloads").glob("*.tif"))
    exact_duplicates = []
    compared = []
    for prior_path in prior_paths:
        if prior_path.resolve() == tif_path.resolve():
            continue
        try:
            with rasterio.open(prior_path) as prior:
                if prior.shape != footprint.shape:
                    continue
                previous = prior.read(1)[footprint]
            previous = np.nan_to_num(previous, nan=-999.0, posinf=-999.0, neginf=-999.0)
            compared.append(prior_path.name)
            if np.array_equal(previous, mask[footprint].astype(np.float32)):
                exact_duplicates.append(prior_path.name)
        except Exception:
            continue
    uniqueness = dict(
        scope="exact array comparison against same-grid GeoTIFFs already present in this repository's docs/downloads",
        compared_files=compared,
        compared_count=len(compared),
        exact_duplicates=exact_duplicates,
        unique_within_checked_local_artifacts=not exact_duplicates,
        limitation="The prior 241-file sibling-repository corpus was not restored or rerun in this validation.",
    )
    if exact_duplicates:
        tif_path.unlink(missing_ok=True)
        raise ValueError(f"Candidate duplicates a checked local artifact: {exact_duplicates}")

    with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(tif_path, arcname=tif_path.name)
    with zipfile.ZipFile(zip_path) as archive:
        bad_member = archive.testzip()
        if bad_member is not None or archive.namelist() != [tif_path.name]:
            raise ValueError(f"ZIP audit failed: bad_member={bad_member!r}, names={archive.namelist()}")

    note = (
        "GEMSDOE39 H39Y-01 radiometric-conductivity cross-gradient | preregistered holdout SPRT accept-H1; "
        "24,000 binary dots; all values in [0,1]"
    )
    if len(note) > 200:
        raise AssertionError(f"Submission note exceeds 200 characters ({len(note)})")
    audit_record = dict(audit)
    audit_record["path"] = str(tif_path.relative_to(ROOT))
    audit_record["sha256"] = sha256_file(tif_path)
    manifest = dict(
        schema_version=1,
        name=name,
        submission_name=name.upper(),
        submission_note=note,
        hypothesis_id=HYPOTHESIS_ID,
        preregistration="docs/preregistration-h39y-20261005.md plus both dated pre-score addenda",
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        method=dict(
            detector="sqrt(R*C) times absolute cosine of strongest GeoDAWN ratio gradient and cond_surf gradient",
            ratio_bands=["ThK", "UK", "UTh"],
            conductivity_band="cond_surf",
            gaussian_sigma_px=2.0,
            gradient_quantiles=[0.02, 0.995],
            budget=EMISSION_BUDGET,
            min_distance_px=MIN_DISTANCE_PX,
            visible_catalogue_emission_buffer_px=KNOWN_BUFFER_PX,
            encoding="binary float32 0/1 inside footprint, 0 outside, no nodata tag",
            source_prediction_copied=False,
        ),
        artifact=dict(
            tif=audit_record,
            zip=dict(path=str(zip_path.relative_to(ROOT)), sha256=sha256_file(zip_path),
                     bytes=zip_path.stat().st_size, one_geotiff_only=True),
            prediction_pixels=int(mask.sum()),
            packed_prediction_sha256=packed_digest,
        ),
        uniqueness=uniqueness,
        official_competition_score=None,
        note="A local holdout result is not a DrivenData leaderboard score.",
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, default=_json_default, allow_nan=False) + "\n")
    return dict(
        name=name,
        submission_name=name.upper(),
        note=note,
        tif_path=str(tif_path.relative_to(ROOT)),
        zip_path=str(zip_path.relative_to(ROOT)),
        manifest_path=str(manifest_path.relative_to(ROOT)),
        tif_sha256=sha256_file(tif_path),
        zip_sha256=sha256_file(zip_path),
        pixels=int(mask.sum()),
        audit=audit,
        uniqueness=uniqueness,
    )


def _build_report(
    status: str,
    input_info: dict[str, Any],
    grid_record: dict[str, Any],
    cells: list[Cell],
    component_count: int,
    train_counts: dict[int, dict[str, int]],
    results: list[dict[str, Any]],
    sprt: PairwiseSPRT,
    started_utc: str,
    candidate_artifact: dict[str, Any] | None,
    failure: str | None = None,
    feature_names: list[str] | None = None,
    candidate_field_sha256: str | None = None,
) -> dict[str, Any]:
    eligible_cells = [cell for cell in cells if cell.eligible]
    code_paths = {
        "evaluator": Path(__file__).resolve(),
        "radiometric": ROOT / "src" / "gems39" / "radiometric.py",
        "stack": ROOT / "src" / "gems39" / "stack.py",
        "emission": ROOT / "src" / "gems39" / "emission.py",
        "metric": ROOT / "src" / "gems39" / "metric.py",
        "sprt": ROOT / "src" / "gems39" / "sprt_select.py",
    }
    return dict(
        schema_version=1,
        candidate_id=HYPOTHESIS_ID,
        status=status,
        started_utc=started_utc,
        completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        registration_files=[
            "docs/preregistration-h39y-20261005.md",
            "docs/preregistration-addendum-h39y-20261005.md",
            "docs/preregistration-addendum2-h39y-20261005.md",
        ],
        implementation_sha256={name: _sha256(path) for name, path in code_paths.items()},
        inputs=dict(
            raster_grid=input_info["reference"],
            footprint_pixels=input_info["footprint_pixels"],
            catalogue_pixels=input_info["catalogue_pixels"],
            outside_footprint_label_pixels=input_info["raw_labels_outside_footprint"],
            sha256={key: _sha256(value) for key, value in input_info["paths"].items()},
        ),
        geometry={
            **grid_record,
            "eligible_truth_min": MIN_TRUTH_PX,
            "eligible_count": len(eligible_cells),
            "required_eligible_count": MIN_ELIGIBLE_CELLS,
            "eligible_order": [cell.ordinal for cell in eligible_cells],
            "cells": [cell.as_dict() for cell in cells],
            "full_catalogue_component_count": component_count,
        },
        comparator=dict(
            id="H39-A-model",
            feature_fields="all sorted fields returned by checked-in stack.build_stack(quantise=True)",
            feature_names=feature_names or [],
            feature_field_count=len(feature_names or []),
            model="HistGradientBoostingClassifier",
            max_iter=250,
            learning_rate=0.08,
            max_leaf_nodes=31,
            l2_regularization=1.0,
            early_stopping=False,
            max_rows_per_class=MAX_TRAIN_ROWS_PER_CLASS,
            minimum_rows_per_class=MIN_TRAIN_ROWS_PER_CLASS,
            negative_pools="training-only catalogue 1-6 px hard negatives OR outside its 8-iteration dilation",
            training_counts={str(k): v for k, v in train_counts.items()},
        ),
        candidate=dict(
            id=HYPOTHESIS_ID,
            score="sqrt(R*C)*abs(cos(theta_R-theta_C))",
            ratios=["ThK", "UK", "UTh"],
            conductivity="cond_surf",
            sigma_px=2.0,
            gradient_quantiles=[0.02, 0.995],
            computed_field_sha256=candidate_field_sha256,
        ),
        emitter=dict(
            type="checked-in best-first Poisson disk",
            budget_full_footprint=EMISSION_BUDGET,
            min_distance_px=MIN_DISTANCE_PX,
            visible_catalogue_buffer_px=KNOWN_BUFFER_PX,
            per_tile_target="round(24,000 * active tile pixels / total footprint pixels)",
        ),
        scorer=dict(
            implementation="src/gems39/metric.py:dti_binary",
            known_fault_scoring_mask="exact visible catalogue pixels only; no buffer, per DrivenData staff clarification",
            radius_px=3.0,
            alpha_fp=0.2,
            beta_fn=0.8,
            strict_win="candidate cell DTI > comparator cell DTI; ties are losses",
        ),
        sprt=sprt.as_dict(),
        evaluated_tiles=results,
        pooled_dti=dict(
            candidate=_pooled_dti(results, "candidate_dti"),
            comparator=_pooled_dti(results, "comparator_dti"),
        ),
        promotion_gate=dict(
            weekly_slot_recommendation=(status == "accept_H1" and candidate_artifact is not None),
            candidate_artifact_created=candidate_artifact is not None,
            decision="recommend only after registered SPRT accept_H1; otherwise keep incumbent and no slot",
            artifact=candidate_artifact,
        ),
        limitation=(
            "The public USGS/INGENIOUS catalogue is a development proxy, not the private new-fault labels. "
            "Spatially separated tiles do not prove independent Bernoulli outcomes; the nominal SPRT error "
            "guarantees are conditional on the preregistered independence/conditional-supermartingale model. "
            "No DrivenData leaderboard score is inferred."
        ),
        failure=failure,
        software=dict(
            python=platform.python_version(),
            numpy=np.__version__,
            scipy=scipy.__version__,
            scikit_learn=sklearn.__version__,
            rasterio=rasterio.__version__,
        ),
    )


def _write_report(out_path: Path, report: dict[str, Any]) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, default=_json_default, allow_nan=False) + "\n")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "docs" / "downloads")
    parser.add_argument("--plan-only", action="store_true", help="Check cell geometry and truth counts only")
    parser.add_argument("--preflight-only", action="store_true", help="Check training pools and data grids; never compute candidate scores")
    args = parser.parse_args()

    started = dt.datetime.now(dt.timezone.utc).isoformat()
    print("[1/6] Read and cross-check the registered grids and catalogue", flush=True)
    footprint, catalogue, input_info = _read_inputs(args.data_dir)
    input_info["paths"] = {
        "sample_submission.tif": input_info["sample"],
        "labels.tif": input_info["labels"],
        "training_features.tif": input_info["features"],
        "geodawn_extensions_u8.tif": input_info["ratios"],
    }
    print(
        f"  footprint={input_info['footprint_pixels']:,}; catalogue={input_info['catalogue_pixels']:,}; "
        f"grid={footprint.shape}; CRS={input_info['reference']['crs']}",
        flush=True,
    )

    print("[2/6] Construct the frozen 8×8 cells; score nothing unless eligibility passes", flush=True)
    cells, grid_record, components = _make_cells(footprint, catalogue)
    eligible_cells = [cell for cell in cells if cell.eligible]
    print("  truth counts by tested-cell ordinal:", [cell.n_truth for cell in cells], flush=True)
    print(f"  eligible={len(eligible_cells)}/{len(cells)}; required={MIN_ELIGIBLE_CELLS}", flush=True)
    if args.plan_only:
        print(json.dumps(dict(grid=grid_record, cells=[c.as_dict() for c in cells]), indent=2), flush=True)
        return 0

    train_counts: dict[int, dict[str, int]] = {}
    if len(eligible_cells) >= MIN_ELIGIBLE_CELLS:
        print("  checking every eligible fold's training-only label pools before candidate computation", flush=True)
        pool_failure = None
        for cell in eligible_cells:
            _, _, _, counts = _training_pools(cell, footprint, catalogue, components)
            train_counts[cell.ordinal] = counts
            print(
                f"  cell {cell.ordinal:02d} ({cell.row},{cell.col}): "
                f"truth={cell.n_truth}; training pos={counts['train_positive']:,}; "
                f"neg={counts['train_negative']:,}",
                flush=True,
            )
            if min(counts["train_positive"], counts["train_negative"]) < MIN_TRAIN_ROWS_PER_CLASS:
                pool_failure = (
                    f"cell ordinal {cell.ordinal} has fewer than {MIN_TRAIN_ROWS_PER_CLASS} "
                    "training examples in one class"
                )
        if pool_failure:
            status = "data_failure"
            sprt = PairwiseSPRT(**SPRT_SETTINGS)
            report = _build_report(
                status, input_info, grid_record, cells, int(components.max()), train_counts,
                [], sprt, started, None, failure=pool_failure,
            )
            stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            report_path = ROOT / "docs" / "reports" / f"h39y01-validation-{stamp}.json"
            _write_report(report_path, report)
            print(f"DATA FAILURE: {pool_failure}; report={report_path.relative_to(ROOT)}", flush=True)
            return 2

    if len(eligible_cells) < MIN_ELIGIBLE_CELLS:
        sprt = PairwiseSPRT(**SPRT_SETTINGS)
        report = _build_report(
            "not_testable", input_info, grid_record, cells, int(components.max()),
            train_counts, [], sprt, started, None,
            failure=f"Only {len(eligible_cells)} cells meet the {MIN_TRUTH_PX}-truth-pixel rule",
        )
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        report_path = ROOT / "docs" / "reports" / f"h39y01-validation-{stamp}.json"
        _write_report(report_path, report)
        print(f"NOT TESTABLE: fewer than {MIN_ELIGIBLE_CELLS} eligible cells; no candidate was scored", flush=True)
        return 0

    source_band_order = _check_source_band_order(input_info["features"], input_info["reference"])
    if args.preflight_only:
        print("[PREFLIGHT] Read and validate aligned detector rasters without computing a candidate score", flush=True)
        ratio_fields = radiometric.read_ratio_bands(
            input_info["ratios"], footprint, reference_path=input_info["sample"]
        )
        conductivity = radiometric.read_conductivity(
            input_info["features"], footprint, reference_path=input_info["sample"]
        )
        ratio_ranges = {
            name: [float(layer[footprint].min()), float(layer[footprint].max())]
            for name, layer in ratio_fields.items()
        }
        conductivity_range = [
            float(conductivity[footprint].min()),
            float(conductivity[footprint].max()),
        ]
        print(
            f"  ratios={sorted(ratio_fields)}; conductivity finite={bool(np.isfinite(conductivity[footprint]).all())}; "
            f"eligible cells={len(eligible_cells)}; training-pool checks passed",
            flush=True,
        )
        del ratio_fields, conductivity
        source_hashes = {key: _sha256(path) for key, path in input_info["paths"].items()}
        code_paths = {
            "evaluator": Path(__file__).resolve(),
            "radiometric": ROOT / "src" / "gems39" / "radiometric.py",
            "stack": ROOT / "src" / "gems39" / "stack.py",
            "emission": ROOT / "src" / "gems39" / "emission.py",
            "metric": ROOT / "src" / "gems39" / "metric.py",
            "sprt": ROOT / "src" / "gems39" / "sprt_select.py",
        }
        preflight = dict(
            candidate_id=HYPOTHESIS_ID,
            checked_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            candidate_score_computed=False,
            holdout_dti_computed=False,
            candidate_tif_generated=False,
            registration_files=[
                "docs/preregistration-h39y-20261005.md",
                "docs/preregistration-addendum-h39y-20261005.md",
                "docs/preregistration-addendum2-h39y-20261005.md",
            ],
            implementation_sha256={name: _sha256(path) for name, path in code_paths.items()},
            input_sha256=source_hashes,
            grid=input_info["reference"],
            footprint_pixels=input_info["footprint_pixels"],
            catalogue_pixels=input_info["catalogue_pixels"],
            row_edges=grid_record["row_edges"],
            col_edges=grid_record["col_edges"],
            cells=[cell.as_dict() for cell in cells],
            eligible_cell_ordinals=[cell.ordinal for cell in eligible_cells],
            eligible_count=len(eligible_cells),
            required_eligible_count=MIN_ELIGIBLE_CELLS,
            full_catalogue_component_count=int(components.max()),
            training_pools={str(k): value for k, value in train_counts.items()},
            checked_feature_band_order=source_band_order,
            checked_ratio_ranges_0_1=ratio_ranges,
            checked_raw_conductivity_range=conductivity_range,
            status="preflight_passed_no_candidate_scores",
        )
        preflight_path = ROOT / "docs" / "reports" / "h39y01-preflight-20261005.json"
        _write_report(preflight_path, preflight)
        print(f"  pre-score plan recorded at {preflight_path.relative_to(ROOT)}", flush=True)
        return 0

    print("[3/6] Read the frozen layers and compute H39Y-01 once", flush=True)
    ratio_fields = radiometric.read_ratio_bands(
        input_info["ratios"], footprint, reference_path=input_info["sample"]
    )
    conductivity = radiometric.read_conductivity(
        input_info["features"], footprint, reference_path=input_info["sample"]
    )
    candidate_field = radiometric.radiometric_conductivity_score(
        ratio_fields, conductivity, footprint
    )
    del ratio_fields, conductivity
    print(
        f"  score range in footprint={float(candidate_field[footprint].min()):.6g}.."
        f"{float(candidate_field[footprint].max()):.6g}; positive={int((candidate_field[footprint] > 0).sum()):,}",
        flush=True,
    )

    cube, feature_names = _prepare_stack(input_info["features"], footprint)
    candidate_field_sha256 = hashlib.sha256(
        np.ascontiguousarray(candidate_field, dtype=np.float32).tobytes()
    ).hexdigest()
    sprt = PairwiseSPRT(**SPRT_SETTINGS)
    results: list[dict[str, Any]] = []
    progress_path = args.data_dir / "h39y01-progress.json"
    for cell in eligible_cells:
        print(
            f"[4/6] Fit comparator for cell {cell.ordinal:02d} ({cell.row},{cell.col}), "
            f"seed={cell.seed}, held-out truth={cell.n_truth}",
            flush=True,
        )
        fold_start = time.monotonic()
        train_mask, train_catalogue, negative_pool, current_counts = _training_pools(
            cell, footprint, catalogue, components
        )
        if current_counts != train_counts[cell.ordinal]:
            raise ValueError(f"Training pool counts changed for cell {cell.ordinal}")
        comparator_field, fit_record = _fit_predict_h39a(
            cube, cell, train_mask, train_catalogue, negative_pool
        )
        del train_mask, negative_pool
        tile_result, masks = _score_tile(
            candidate_field, comparator_field, cell, train_catalogue, input_info["footprint_pixels"]
        )
        sprt_state = sprt.update(tile_result["candidate_strict_win"])
        tile_result["comparator_fit"] = fit_record
        tile_result["sprt_after_tile"] = sprt_state
        tile_result["elapsed_seconds"] = time.monotonic() - fold_start
        results.append(tile_result)
        print(
            f"  candidate DTI={tile_result['candidate_dti']['dti']:.6f}; "
            f"H39-A DTI={tile_result['comparator_dti']['dti']:.6f}; "
            f"win={tile_result['candidate_strict_win']}; LLR={sprt_state['llr']:+.6f}; "
            f"decision={sprt_state['decision']}",
            flush=True,
        )
        progress = dict(
            candidate_id=HYPOTHESIS_ID,
            status=sprt_state["decision"],
            scored_tiles=results,
            sprt=sprt_state,
            note="Progress only; final report is written after the predeclared stop.",
        )
        _write_progress(progress_path, progress)
        del comparator_field, masks, train_catalogue
        if sprt.decision != "continue":
            print("  SPRT boundary crossed; no later eligible tile will be fit or scored.", flush=True)
            break

    if sprt.decision == "accept_H1":
        status = "accept_H1"
    elif sprt.decision == "accept_H0":
        status = "accept_H0"
    else:
        status = "inconclusive"
    artifact = None
    if status == "accept_H1":
        artifact = _emit_submission(
            candidate_field,
            footprint,
            catalogue,
            input_info["sample"],
            args.out_dir,
        )
    else:
        print(f"[6/6] {status}: retaining the incumbent; no candidate TIFF or weekly slot", flush=True)

    report = _build_report(
        status,
        input_info,
        grid_record,
        cells,
        int(components.max()),
        train_counts,
        results,
        sprt,
        started,
        artifact,
        feature_names=feature_names,
        candidate_field_sha256=candidate_field_sha256,
    )
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_path = ROOT / "docs" / "reports" / f"h39y01-validation-{stamp}.json"
    _write_report(report_path, report)
    print(f"REPORT: {report_path.relative_to(ROOT)}", flush=True)
    print(
        f"FINAL: {status}; wins={sprt.wins}/{sprt.n}; LLR={sprt.llr:+.6f}; "
        f"candidate pooled DTI={report['pooled_dti']['candidate']}; "
        f"comparator pooled DTI={report['pooled_dti']['comparator']}; "
        f"weekly_slot_recommendation={report['promotion_gate']['weekly_slot_recommendation']}",
        flush=True,
    )
    if artifact:
        print(f"TIFF: {artifact['tif_path']} SHA256={artifact['tif_sha256']}", flush=True)
        print(f"ZIP:  {artifact['zip_path']} SHA256={artifact['zip_sha256']}", flush=True)
        print(f"NAME: {artifact['submission_name']}", flush=True)
        print(f"NOTE: {artifact['note']} ({len(artifact['note'])}/200 chars)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
