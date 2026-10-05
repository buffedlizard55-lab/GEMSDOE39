"""Deterministic spatially blocked proxy holdout for the GEMS label raster.

This is a development proxy: public catalogue faults are not the competition's
private newly identified faults. Entire 8-connected components are hidden from
the catalogue mask, folds are geographically separated by an inward collar,
and no held-out component is used to generate that fold's candidate predictions.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_erosion, label

from .metric import dti_binary

FOLD_NAMES = ("NW", "NE", "SW", "SE")
COLLAR_PX = 15
DEFAULT_HIDE_FRACTION = 0.20
DEFAULT_SEED = 20261005


@dataclass
class Fold:
    name: str
    bbox: tuple[slice, slice]
    domain: np.ndarray
    truth: np.ndarray
    known: np.ndarray
    train_catalogue: np.ndarray
    hidden_component_ids: np.ndarray
    n_truth: int


@dataclass
class Holdout:
    footprint: np.ndarray
    labels: np.ndarray
    folds: list[Fold]
    component_count: int
    collar_px: int
    hide_fraction: float
    seed: int


def _quadrant_ids(footprint: np.ndarray) -> tuple[np.ndarray, int, int]:
    rows, cols = np.nonzero(footprint)
    if rows.size == 0:
        raise ValueError("Competition footprint is empty")
    y_mid, x_mid = int(np.median(rows)), int(np.median(cols))
    H, W = footprint.shape
    y, x = np.ogrid[:H, :W]
    quadrants = np.full(footprint.shape, -1, dtype=np.int8)
    quadrants[(y < y_mid) & (x < x_mid) & footprint] = 0
    quadrants[(y < y_mid) & (x >= x_mid) & footprint] = 1
    quadrants[(y >= y_mid) & (x < x_mid) & footprint] = 2
    quadrants[(y >= y_mid) & (x >= x_mid) & footprint] = 3
    return quadrants, y_mid, x_mid


def load_holdout(
    data_dir: str | Path,
    hide_fraction: float = DEFAULT_HIDE_FRACTION,
    seed: int = DEFAULT_SEED,
    collar_px: int = COLLAR_PX,
) -> Holdout:
    """Create four fixed component-hiding folds before scoring any candidate."""
    data_dir = Path(data_dir)
    if not 0.0 < hide_fraction < 1.0:
        raise ValueError("hide_fraction must be between 0 and 1")
    if collar_px < 1:
        raise ValueError("collar_px must be at least 1 for spatial blocking")
    with rasterio.open(data_dir / "sample_submission.tif") as sample, rasterio.open(
        data_dir / "labels.tif"
    ) as label_ds:
        if (
            sample.width != label_ds.width
            or sample.height != label_ds.height
            or sample.crs != label_ds.crs
            or tuple(sample.transform)[:6] != tuple(label_ds.transform)[:6]
        ):
            raise ValueError("Sample and labels do not share the same raster grid")
        footprint = np.isfinite(sample.read(1))
        labels = (label_ds.read(1) >= 1) & footprint

    structure = np.ones((3, 3), dtype=np.uint8)
    components, n_components = label(labels, structure=structure)
    quadrants, y_mid, x_mid = _quadrant_ids(footprint)

    # Assign each complete component to exactly one fold by its centroid. This
    # avoids reusing one connected trace as two nominally separate outcomes.
    y_label, x_label = np.nonzero(labels)
    ids_at_label = components[y_label, x_label]
    counts = np.bincount(ids_at_label, minlength=n_components + 1)
    sum_y = np.bincount(ids_at_label, weights=y_label, minlength=n_components + 1)
    sum_x = np.bincount(ids_at_label, weights=x_label, minlength=n_components + 1)
    ids = np.arange(1, n_components + 1, dtype=np.int32)
    centroid_y = np.zeros(n_components + 1, dtype=np.float64)
    centroid_x = np.zeros(n_components + 1, dtype=np.float64)
    nonzero = counts > 0
    centroid_y[nonzero] = sum_y[nonzero] / counts[nonzero]
    centroid_x[nonzero] = sum_x[nonzero] / counts[nonzero]
    component_fold = np.full(n_components + 1, -1, dtype=np.int8)
    component_fold[(centroid_y < y_mid) & (centroid_x < x_mid)] = 0
    component_fold[(centroid_y < y_mid) & (centroid_x >= x_mid)] = 1
    component_fold[(centroid_y >= y_mid) & (centroid_x < x_mid)] = 2
    component_fold[(centroid_y >= y_mid) & (centroid_x >= x_mid)] = 3

    folds: list[Fold] = []
    for fold_id, fold_name in enumerate(FOLD_NAMES):
        quadrant = quadrants == fold_id
        domain_full = binary_erosion(quadrant, iterations=collar_px, border_value=0) & footprint
        rr, cc = np.nonzero(domain_full)
        if rr.size == 0:
            raise ValueError(f"Spatial fold {fold_name} is empty after the collar")
        bbox = (slice(int(rr.min()), int(rr.max()) + 1), slice(int(cc.min()), int(cc.max()) + 1))
        domain = domain_full[bbox]
        label_domain = domain & labels[bbox]
        component_sizes = np.bincount(components[bbox][label_domain], minlength=n_components + 1)
        eligible = ids[(component_fold[ids] == fold_id) & (component_sizes[ids] > 0)]
        if eligible.size == 0:
            raise ValueError(f"No label components available in holdout fold {fold_name}")
        target_pixels = max(1, int(round(hide_fraction * int(label_domain.sum()))))
        rng = np.random.default_rng(seed + fold_id)
        shuffled = eligible[rng.permutation(eligible.size)]
        cumulative = np.cumsum(component_sizes[shuffled])
        n_take = min(shuffled.size, int(np.searchsorted(cumulative, target_pixels, side="left")) + 1)
        hidden_ids = shuffled[:n_take]
        selected = np.zeros(n_components + 1, dtype=bool)
        selected[hidden_ids] = True
        hidden_full = selected[components]
        train_catalogue = labels & ~hidden_full
        hidden_sub = hidden_full[bbox]
        known_sub = labels[bbox] & ~hidden_sub
        truth = hidden_sub & domain
        folds.append(
            Fold(
                name=fold_name,
                bbox=bbox,
                domain=domain,
                truth=truth,
                known=known_sub,
                train_catalogue=train_catalogue,
                hidden_component_ids=hidden_ids.copy(),
                n_truth=int(truth.sum()),
            )
        )

    return Holdout(
        footprint=footprint,
        labels=labels,
        folds=folds,
        component_count=int(n_components),
        collar_px=int(collar_px),
        hide_fraction=float(hide_fraction),
        seed=int(seed),
    )


def evaluate_fold(prediction: np.ndarray, fold: Fold) -> dict:
    """Score predictions only inside one held-out spatial domain."""
    prediction = np.asarray(prediction, bool)
    pred_sub = prediction[fold.bbox]
    return dti_binary(pred_sub, fold.truth, valid=fold.domain, known=fold.known)
