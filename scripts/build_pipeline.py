#!/usr/bin/env python3
"""Reproduce and audit the exploratory H39X-01 candidate.

The exact detector formula was recorded in the hypothesis register only after
its 2026-10-05 holdout result. This script can reproduce that exposed analysis,
but it can NEVER approve H39X-01 for a weekly slot, even if a future rerun
changes the SPRT state. Fixed test settings do not retroactively preregister a
candidate. This script neither uploads nor consumes a submission slot.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import emission, features, grid, holdout  # noqa: E402
from gems39.metric import dti  # noqa: E402
from gems39.sprt_select import PairwiseSPRT  # noqa: E402

FEATURE_BANDS = [
    "geod_2ndinv",
    "geod_shearrate",
    "geod_dilaterate",
    "det_elev",
    "tmi",
    "det_elev_slope",
]
HYPOTHESIS_ID = "H39X-01-strain-dilatation-corridor"
DEFAULT_BUDGET = 44_090
MIN_DISTANCE_PX = 2.7
HIDE_FRACTION = 0.20
HOLDOUT_SEED = 20261005
COLLAR_PX = 15
SPRT_P0 = 0.50
SPRT_P1 = 0.70
SPRT_ALPHA = 0.05
SPRT_BETA = 0.10


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _manifest_hashes() -> dict[str, str]:
    path = ROOT / "registry" / "data_manifest.json"
    if not path.exists():
        return {}
    manifest = json.loads(path.read_text())
    return {item["id"]: item["sha256"] for item in manifest.get("files", [])}


def _check_no_duplicate_pixels(candidate: Path, footprint: np.ndarray) -> dict:
    """Check local and restored historical GeoTIFFs for identical in-footprint values."""
    search_roots = [ROOT / "docs" / "downloads", ROOT / "artifacts", ROOT / "data" / "scored"]
    compared = []
    duplicate_of = []
    with rasterio.open(candidate) as new:
        candidate_values = new.read(1)[footprint]
    candidate_values = np.nan_to_num(candidate_values, nan=-999.0, posinf=-999.0, neginf=-999.0)
    for folder in search_roots:
        if not folder.exists():
            continue
        for path in sorted(folder.glob("*.tif")):
            if path.resolve() == candidate.resolve():
                continue
            try:
                with rasterio.open(path) as prior:
                    if prior.shape != footprint.shape:
                        continue
                    prior_values = prior.read(1)[footprint]
                prior_values = np.nan_to_num(prior_values, nan=-999.0, posinf=-999.0, neginf=-999.0)
                compared.append(str(path.relative_to(ROOT)))
                if np.array_equal(candidate_values, prior_values):
                    duplicate_of.append(str(path.relative_to(ROOT)))
            except Exception:
                # A corrupt/unreadable prior artifact is not treated as proof of uniqueness.
                continue
    return {
        "compared_artifacts": compared,
        "compared_count": len(compared),
        "duplicate_of": duplicate_of,
        "unique_vs_checked_local_artifacts": not duplicate_of,
        "scope_note": "This local check does not inspect every file in all sibling GEMSDOE repositories.",
    }


def _pooled_dti(results: list[dict]) -> float:
    tp = sum(float(r["tp"]) for r in results)
    fp = sum(float(r["fp"]) for r in results)
    fn = sum(float(r["fn"]) for r in results)
    return dti(tp, fp, fn)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=str(ROOT / "data"))
    parser.add_argument("--out-dir", default=str(ROOT / "docs" / "downloads"))
    parser.add_argument("--tag", default=None, help="Optional unique filename suffix.")
    parser.add_argument(
        "--replay-exposed-holdout",
        action="store_true",
        help="Recompute the exact locked exploratory run for audit only; never slot-eligible.",
    )
    args = parser.parse_args()
    if not args.replay_exposed_holdout:
        parser.error(
            "The H39X-01 holdout is exhausted. Refusing to re-score by default; "
            "the manifest already records the result. Use --replay-exposed-holdout "
            "only for exact deterministic reproduction, never for variant selection."
        )

    print("NOTE: this is a locked replay of an exposed exploratory run; it cannot approve a submission.")

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    required = ["training_features.tif", "labels.tif", "sample_submission.tif"]
    missing = [name for name in required if not (data_dir / name).is_file()]
    if missing:
        raise SystemExit(
            "Missing competition inputs in " + str(data_dir) + ": " + ", ".join(missing)
            + ". Restore the SHA-256-pinned owner mirror with scripts/restore_data.py --group core, "
            + "or place the authenticated DrivenData downloads there."
        )

    print("[1/6] Read sample footprint and explicitly named feature bands")
    footprint = grid.read_footprint_from_sample(data_dir / "sample_submission.tif")
    print(f"  footprint={int(footprint.sum()):,}/{footprint.size:,} cells")
    bands = grid.read_bands(data_dir / "training_features.tif", footprint, names=FEATURE_BANDS)
    candidate_field = features.build_h39x01_strain(bands, footprint)
    incumbent_field = features.build_incumbent_ridge(bands, footprint)
    del bands

    print("[2/6] Reconstruct the exposed exploratory spatial holdout (not a confirmation set)")
    context = holdout.load_holdout(
        data_dir,
        hide_fraction=HIDE_FRACTION,
        seed=HOLDOUT_SEED,
        collar_px=COLLAR_PX,
    )
    print(
        f"  components={context.component_count:,}; folds={len(context.folds)}; "
        f"collar={context.collar_px} px ({context.collar_px * 100:,} m); "
        f"hide_fraction={context.hide_fraction:.2f}; seed={context.seed}"
    )

    print("[3/6] Evaluate candidate and fixed incumbent fold-by-fold; update SPRT immediately")
    sprt = PairwiseSPRT(
        p0=SPRT_P0,
        p1=SPRT_P1,
        alpha=SPRT_ALPHA,
        beta=SPRT_BETA,
    )
    fold_records = []
    candidate_fold_scores = []
    incumbent_fold_scores = []
    candidate_components = []
    incumbent_components = []
    for fold in context.folds:
        candidate_mask = emission.emit_fast(
            candidate_field,
            context.footprint,
            target_n=DEFAULT_BUDGET,
            min_dist=MIN_DISTANCE_PX,
            catalogue=fold.train_catalogue,
            cat_buffer_px=0,
        )
        incumbent_mask = emission.emit_fast(
            incumbent_field,
            context.footprint,
            target_n=DEFAULT_BUDGET,
            min_dist=MIN_DISTANCE_PX,
            catalogue=fold.train_catalogue,
            cat_buffer_px=0,
        )
        candidate_result = holdout.evaluate_fold(candidate_mask, fold)
        incumbent_result = holdout.evaluate_fold(incumbent_mask, fold)
        win = float(candidate_result["dti"]) > float(incumbent_result["dti"]) + 1e-12
        sprt_state = sprt.update(win)
        candidate_fold_scores.append(candidate_result["dti"])
        incumbent_fold_scores.append(incumbent_result["dti"])
        candidate_components.append(candidate_result)
        incumbent_components.append(incumbent_result)
        record = {
            "fold": fold.name,
            "hidden_components": int(fold.hidden_component_ids.size),
            "hidden_truth_pixels": fold.n_truth,
            "candidate": candidate_result,
            "incumbent": incumbent_result,
            "candidate_wins_fold": bool(win),
            "sprt_after_fold": sprt_state,
            "candidate_emitted_pixels": int(candidate_mask.sum()),
            "incumbent_emitted_pixels": int(incumbent_mask.sum()),
        }
        fold_records.append(record)
        print(
            f"  {fold.name}: candidate={candidate_result['dti']:.6f}, "
            f"incumbent={incumbent_result['dti']:.6f}, win={win}, "
            f"LLR={sprt_state['llr']:.4f}/{sprt_state['upper']:.4f}, "
            f"decision={sprt_state['decision']}"
        )
        del candidate_mask, incumbent_mask
        if sprt_state["decision"] != "continue":
            print("  SPRT boundary crossed; no additional folds were evaluated.")
            break

    candidate_pooled = _pooled_dti(candidate_components)
    incumbent_pooled = _pooled_dti(incumbent_components)
    candidate_mean = float(np.mean(candidate_fold_scores)) if candidate_fold_scores else 0.0
    incumbent_mean = float(np.mean(incumbent_fold_scores)) if incumbent_fold_scores else 0.0
    beats_pooled = candidate_pooled > incumbent_pooled + 1e-12
    nominal_screen_pass = sprt.decision == "accept_H1" and beats_pooled
    # The exact detector formula was documented post-hoc. Never convert a
    # rerun on the same exposed folds into a submission recommendation.
    eligible = False
    print("[4/6] Decision")
    print(f"  pooled DTI proxy: candidate={candidate_pooled:.6f}, incumbent={incumbent_pooled:.6f}")
    print(f"  unweighted block mean: candidate={candidate_mean:.6f}, incumbent={incumbent_mean:.6f}")
    print(f"  SPRT={sprt.decision}; nominal_screen_pass={nominal_screen_pass}; "
          f"recipe_registration=posthoc; eligible_for_submission_slot={eligible}")

    print("[5/6] Emit a fresh full-footprint candidate (no prior prediction pixels are reused)")
    labels = grid.read_labels(data_dir / "labels.tif") & footprint
    final_mask = emission.emit_fast(
        candidate_field,
        footprint,
        target_n=DEFAULT_BUDGET,
        min_dist=MIN_DISTANCE_PX,
        catalogue=labels,
        cat_buffer_px=0,
    )
    pixel_digest = _sha256_bytes(np.packbits(final_mask, bitorder="little").tobytes())
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = re.sub(r"[^A-Za-z0-9-]+", "-", args.tag).strip("-") if args.tag else timestamp
    if not suffix:
        raise SystemExit("--tag must contain at least one alphanumeric character")
    # Hash suffix avoids collisions if a rerun occurs in the same UTC second.
    name = f"gemsdoe39-h39x01-strain-corridor-{suffix}-{pixel_digest[:8]}"
    output_path = out_dir / f"{name}-nan.tif"
    if output_path.exists():
        raise SystemExit(f"Refusing to overwrite existing candidate: {output_path}")
    prediction = final_mask.astype(np.float32)
    grid.write_submission(
        prediction,
        data_dir / "sample_submission.tif",
        output_path,
        footprint,
        outside="nan",
    )
    audit = grid.audit_submission(output_path, footprint, data_dir / "sample_submission.tif")
    audit["path"] = _display_path(output_path)
    if not audit["ok"]:
        raise SystemExit(f"Generated GeoTIFF failed its independent audit: {audit}")
    duplicate_audit = _check_no_duplicate_pixels(output_path, footprint)
    if duplicate_audit["duplicate_of"]:
        output_path.unlink(missing_ok=True)
        raise SystemExit(f"Candidate prediction pixels duplicate a checked prior artifact: {duplicate_audit}")
    print(f"  path={_display_path(output_path)}")
    print(f"  pixels={audit['positive_px']:,}; sha256={audit['sha256']}; bytes={audit['bytes']:,}")
    print(
        f"  duplicate_check={duplicate_audit['compared_count']} local files; "
        f"unique_vs_checked_local_artifacts={duplicate_audit['unique_vs_checked_local_artifacts']}"
    )

    note = (
        "GEMSDOE39 H39X-01 geodetic strain/dilatation corridor | "
        f"budget {DEFAULT_BUDGET}; exact known-fault mask; exploratory post-hoc candidate; "
        "not eligible for a weekly submission slot."
    )
    manifest = {
        "schema_version": 2,
        "name": name,
        "note_for_drivendata": note,
        "hypothesis_id": HYPOTHESIS_ID,
        "hypothesis_registration": {
            "status": "post_hoc_exploratory",
            "exact_detector_recipe_recorded_after_holdout_scoring": True,
            "confirmatory": False,
            "weekly_slot_eligible": False,
            "note": "Do not describe fixed SPRT settings or the current register as retroactive preregistration.",
        },
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "method": {
            "detector": "robust-rank geodetic second-invariant and shear amplitude; multi-scale strain-invariant line response; paired local dilatation or dilatation-gradient support",
            "feature_bands": FEATURE_BANDS,
            "candidate_feature_bands": ["geod_2ndinv", "geod_shearrate", "geod_dilaterate"],
            "incumbent_feature_bands": ["det_elev", "tmi", "det_elev_slope"],
            "band_preprocessing": {
                "names_from": "GeoTIFF band_name tag or band description; no numeric-order inference",
                "valid_cells": "finite and greater than the float32 nodata sentinel threshold",
                "invalid_in_footprint": "nearest-valid fill",
                "clip_quantiles": [0.001, 0.999],
                "outside_footprint_value": 0.0,
            },
            "feature_transform_parameters": {
                "robust_quantiles": [0.02, 0.995],
                "dilatation_gaussian_sigma_px": 1.5,
                "dilatation_pair_window_px": 7,
                "dilatation_gradient_sigma_px": 3.0,
                "strain_invariant_line_scales_px": [2.0, 4.0],
                "linear_weights": {"strain_amplitude": 0.50, "line_response": 0.30, "dilatation_support": 0.20},
            },
            "emitter": {"type": "best-first Poisson-disk", "target_pixels": DEFAULT_BUDGET, "min_distance_px": MIN_DISTANCE_PX, "known_fault_buffer_px": 0},
            "prediction_encoding": "binary float32 0/1 inside footprint; NaN outside",
            "source_prediction_copied": False,
        },
        "inputs": {
            "provenance": "SHA-256-pinned owner mirror; hashes verify mirror integrity, not organizer authentication",
            "sha256": _manifest_hashes(),
        },
        "holdout": {
            "proxy_only": True,
            "description": (
                "Exposed exploratory proxy holdout: four geographic quadrants; hide whole "
                f"8-connected catalogue components totalling approximately {HIDE_FRACTION:.3f} "
                f"of fold label pixels; {COLLAR_PX * 100:g} m inward collar; seed {HOLDOUT_SEED}; "
                "not the private expert-label test set. Exact H39X-01 formula was recorded after "
                "scoring; do not reuse for confirmation."
            ),
            "fold_order": [record["fold"] for record in fold_records],
            "seed": HOLDOUT_SEED,
            "hide_fraction": HIDE_FRACTION,
            "collar_px": COLLAR_PX,
            "candidate_pooled_dti": candidate_pooled,
            "incumbent_pooled_dti": incumbent_pooled,
            "candidate_unweighted_fold_mean": candidate_mean,
            "incumbent_unweighted_fold_mean": incumbent_mean,
            "candidate_wins": sprt.wins,
            "candidate_losses": sprt.losses,
            "exploratory_only": True,
            "confirmatory": False,
            "fold_results": fold_records,
            "sprt": sprt.as_dict(),
            "statistical_assumption": "SPRT alpha/beta control is conditional on independent Bernoulli block outcomes or a valid conditional supermartingale; spatial collar does not prove this assumption.",
        },
        "multiplicity": {
            "current_pipeline_candidate_count": 1,
            "development_variants_fully_registered": False,
            "current_result_exploratory": True,
            "future_policy": "Use one candidate selected on development data with a fresh holdout, or preallocate family-wise alpha across every candidate on a fresh holdout; never reuse these exposed folds.",
        },
        "submission_gate": {
            "eligible_for_weekly_slot": bool(eligible),
            "reason": (
                "Not eligible: the exact H39X-01 formula was recorded after holdout scoring; the exposed exploratory folds cannot approve this candidate."
                + (" The SPRT was also not accepted." if sprt.decision != "accept_H1" else "")
                + (" Pooled proxy DTI did not beat the comparator." if not beats_pooled else "")
            ),
            "leaderboard_score": None,
        },
        "artifact": audit,
        "uniqueness_audit": duplicate_audit,
        "official_sources": {
            "problem_description": "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/",
            "organizer_known_fault_masking": "https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2",
            "usgs_geodawn": "https://doi.org/10.5066/P93LGLVQ",
        },
    }
    manifest_path = out_dir / f"{name}-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print("[6/6] Complete")
    print(f"  audit_manifest={_display_path(manifest_path)}")
    print(f"  submission_note={note}")
    print(f"  submission_slot_status={'APPROVED' if eligible else 'NOT APPROVED / candidate-only'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
