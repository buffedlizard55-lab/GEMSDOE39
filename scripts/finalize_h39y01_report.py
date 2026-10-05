#!/usr/bin/env python3
"""Finalize the H39Y-01 report from saved sequential outcomes, without rescoring.

The locked evaluator stopped correctly at its SPRT boundary but its final
report assembly hit a duplicate-key TypeError after the boundary. This recovery
script verifies the saved per-tile outcomes and SPRT arithmetic, reconstructs
only the already-registered geometry/metadata, and writes the report. It never
fits a model, computes candidate DTI, or generates a TIFF.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from gems39.sprt_select import PairwiseSPRT  # noqa: E402
from scripts.evaluate_h39y01 import _build_report, _make_cells, _read_inputs, _sha256, _write_report  # noqa: E402

PRECHECK_PATH = ROOT / "docs" / "reports" / "h39y01-preflight-20261005.json"
PROGRESS_PATH = ROOT / "data" / "h39y01-progress.json"
REPORT_PATH = ROOT / "docs" / "reports" / "h39y01-validation-20261005.json"
MARKDOWN_PATH = ROOT / "docs" / "reports" / "h39y01-validation-20261005.md"


def _iso_from_timestamp(path: Path) -> str:
    return dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc).isoformat()


def _markdown(report: dict[str, Any]) -> str:
    rows = [
        "# H39Y-01 spatial holdout result",
        "",
        f"**Decision: {report['status']} — no weekly submission slot recommended.**",
        "",
        "H39Y-01 was the only candidate scored. The frozen sequential test crossed its lower boundary after five strict losses. The incumbent H39-A-model is retained. No H39Y TIFF or ZIP was generated.",
        "",
        "## Pre-score eligibility and stop",
        "",
        f"- Eligible cells: **{report['geometry']['eligible_count']}** of 16 (threshold: {report['geometry']['eligible_truth_min']} hidden-truth pixels; minimum required: {report['geometry']['required_eligible_count']}).",
        f"- Evaluated: **{len(report['evaluated_tiles'])}** cells, in the registered row-major even/even sequence; the test stopped at the first boundary. No later cell was fit or scored.",
        f"- SPRT: `p0={report['sprt']['p0']:.2f}`, `p1={report['sprt']['p1']:.2f}`, `alpha={report['sprt']['alpha']:.2f}`, `beta={report['sprt']['beta']:.2f}`; boundaries `{report['sprt']['lower']:.6f}` / `{report['sprt']['upper']:.6f}`.",
        f"- Final LLR: **{report['sprt']['llr']:.6f}**, below the lower boundary; **{report['sprt']['wins']} wins, {report['sprt']['losses']} losses**.",
        "",
        "## Paired tile results",
        "",
        "| Ordinal | Grid cell (row, col) | Hidden truth px | Active px | Budget | H39Y-01 DTI | H39-A DTI | Strict win? |",
        "|---:|:---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for tile in report["evaluated_tiles"]:
        rows.append(
            f"| {tile['ordinal']} | ({tile['row']}, {tile['col']}) | {tile['truth_pixels']:,} | {tile['active_pixels']:,} | {tile['per_tile_target_pixels']:,} | {tile['candidate_dti']['dti']:.6f} | {tile['comparator_dti']['dti']:.6f} | No |"
        )
    pooled = report["pooled_dti"]
    rows.extend([
        "",
        f"Pooled proxy DTI over the **five evaluated cells only**: H39Y-01 **{pooled['candidate']:.6f}**; H39-A-model **{pooled['comparator']:.6f}**.",
        "",
        "## Interpretation and caveats",
        "",
        "- This is a local, component-hidden public-catalogue proxy, not the private newly identified faults and not an organizer leaderboard measurement.",
        "- The known-fault scoring mask was pixel-exact, matching the DrivenData staff clarification; the 2-pixel known-fault buffer was applied to emitted predictions only.",
        "- The SPRT error guarantees are conditional on independent Bernoulli cell outcomes (or the registered conditional-supermartingale model). A 15-pixel collar reduces spatial leakage but cannot prove geological independence.",
        "- H39Y-01 is not promoted. Do not tune or re-score it on these exposed cells; any future candidate needs a fresh holdout or a preregistered family-wise alpha allocation.",
        "",
        "## Reproducibility and reporting irregularity",
        "",
        f"- Pre-score data/grid/training-pool lock: [`h39y01-preflight-20261005.json`](h39y01-preflight-20261005.json).",
        f"- Machine-readable paired outcomes: [`h39y01-validation-20261005.json`](h39y01-validation-20261005.json).",
        f"- The evaluator reached `accept_H0` and stopped as registered. Its post-stop report assembly then raised a duplicate-key `TypeError`; this report was reconstructed from the saved progress outcomes and pre-score lock. No fold was refit, no DTI was recomputed, and no candidate artifact was generated during recovery.",
        f"- The original evaluator did not persist its start timestamp or continuous candidate-field digest in the progress file; those fields are reported as unavailable rather than reconstructed from outcome data.",
        "",
        "## Source for scoring-mask geometry",
        "",
        "[DrivenData staff clarification, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4): known-fault mask is pixel-exact; nearby predictions are still scored normally; new-fault truth may lie within 300 m of known traces.",
        "",
    ])
    return "\n".join(rows)


def main() -> int:
    if not PRECHECK_PATH.is_file() or not PROGRESS_PATH.is_file():
        raise SystemExit("Required preflight or saved progress file is missing; refusing to infer outcomes")
    preflight = json.loads(PRECHECK_PATH.read_text())
    progress = json.loads(PROGRESS_PATH.read_text())
    if preflight.get("candidate_score_computed") or preflight.get("holdout_dti_computed"):
        raise SystemExit("Preflight record is not a pre-score lock")
    if progress.get("candidate_id") != "H39Y-01" or progress.get("status") != "accept_H0":
        raise SystemExit("Unexpected or incomplete H39Y-01 progress record")

    sprt = PairwiseSPRT(
        p0=preflight.get("sprt", {}).get("p0", 0.50),
        p1=preflight.get("sprt", {}).get("p1", 0.70),
        alpha=preflight.get("sprt", {}).get("alpha", 0.05),
        beta=preflight.get("sprt", {}).get("beta", 0.10),
    )
    outcomes = progress["scored_tiles"]
    for index, tile in enumerate(outcomes):
        state = sprt.update(tile["candidate_strict_win"])
        saved = tile["sprt_after_tile"]
        if state["n"] != saved["n"] or abs(state["llr"] - saved["llr"]) > 1e-12:
            raise SystemExit(f"SPRT replay check failed at saved tile {index}")
        if state["decision"] != saved["decision"]:
            raise SystemExit(f"SPRT decision mismatch at saved tile {index}")
        if state["decision"] != "continue" and index != len(outcomes) - 1:
            raise SystemExit("Progress contains outcomes after a registered stopping boundary")
    if sprt.as_dict() != progress["sprt"]:
        raise SystemExit("Final saved SPRT state does not match replay of saved outcomes")

    footprint, catalogue, input_info = _read_inputs(ROOT / "data")
    input_info["paths"] = {
        "sample_submission.tif": input_info["sample"],
        "labels.tif": input_info["labels"],
        "training_features.tif": input_info["features"],
        "geodawn_extensions_u8.tif": input_info["ratios"],
    }
    current_input_hashes = {
        key: _sha256(path) for key, path in input_info["paths"].items()
    }
    if current_input_hashes != preflight["input_sha256"]:
        raise SystemExit("Input bytes differ from the committed pre-score lock")
    cells, grid_record, components = _make_cells(footprint, catalogue)
    if [cell.as_dict() for cell in cells] != preflight["cells"]:
        raise SystemExit("Reconstructed geometry differs from the committed pre-score cell plan")

    train_counts = {int(key): value for key, value in preflight["training_pools"].items()}
    report = _build_report(
        status="accept_H0",
        input_info=input_info,
        grid_record=grid_record,
        cells=cells,
        component_count=int(components.max()),
        train_counts=train_counts,
        results=outcomes,
        sprt=sprt,
        started_utc="unavailable: original progress record did not persist start time",
        candidate_artifact=None,
        feature_names=[],
        candidate_field_sha256=None,
    )
    # The checked-in evaluator changed only its post-stop report-assembly code
    # after this run. Preserve the exact source hashes committed before scoring.
    report["implementation_sha256"] = preflight["implementation_sha256"]
    report["started_utc"] = None
    report["completed_utc"] = _iso_from_timestamp(PROGRESS_PATH)
    report["comparator"]["feature_names"] = []
    report["comparator"]["feature_field_count"] = 40
    report["comparator"]["feature_names_note"] = "The evaluator did not persist the sorted feature-name list; the exact stack implementation is SHA-256 pinned."
    report["report_finalization"] = dict(
        original_report_assembly_error="TypeError: dict() got multiple values for keyword argument 'full_catalogue_component_count'",
        recovery="Reconstructed solely from the committed pre-score plan and data/h39y01-progress.json; no fit, candidate scoring, DTI computation, or TIFF emission was rerun.",
        outcome_replay_verified=True,
        start_time_persisted=False,
        candidate_field_digest_persisted=False,
        finalizer_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    if REPORT_PATH.exists() or MARKDOWN_PATH.exists():
        raise SystemExit("Final report already exists; refusing to overwrite it")
    _write_report(REPORT_PATH, report)
    MARKDOWN_PATH.write_text(_markdown(report))
    print(f"Replayed {len(outcomes)} saved outcomes; SPRT={sprt.decision}, LLR={sprt.llr:.9f}")
    print(f"No scoring or TIFF generation performed. Report: {REPORT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
