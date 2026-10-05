#!/usr/bin/env python3
"""Generate a simple, explicitly non-trained baseline GeoTIFF.

This is a convenience baseline only; it is not the exploratory H39X-01
candidate and must not be used for a leaderboard slot without a fresh, pre-specified holdout review.
It averages robustly normalized selected feature bands on the competition
footprint and writes NaN outside the sample bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39 import grid  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="sample_submission.tif")
    parser.add_argument("--features", required=True, help="training_features.tif")
    parser.add_argument("--output", required=True, help="output GeoTIFF path")
    parser.add_argument("--bands", default="", help="comma-separated 1-based bands; default all")
    args = parser.parse_args()
    reference = Path(args.reference)
    features = Path(args.features)
    output = Path(args.output)
    if not reference.is_file() or not features.is_file():
        raise SystemExit("Both --reference and --features must exist")

    with rasterio.open(reference) as ref:
        footprint = np.isfinite(ref.read(1))
        reference_signature = (ref.crs, ref.width, ref.height, tuple(ref.transform)[:6])
    with rasterio.open(features) as src:
        source_signature = (src.crs, src.width, src.height, tuple(src.transform)[:6])
        if source_signature != reference_signature:
            raise SystemExit("Reference/features grid, CRS, or transform mismatch")
        bands = [int(item) for item in args.bands.split(",") if item.strip()]
        bands = bands or list(range(1, src.count + 1))
        if any(item < 1 or item > src.count for item in bands):
            raise SystemExit(f"Band index must be in 1..{src.count}")
        total = np.zeros(footprint.shape, dtype=np.float32)
        count = np.zeros(footprint.shape, dtype=np.uint8)
        for band in bands:
            values = src.read(band, masked=True).astype(np.float32).filled(np.nan)
            valid = footprint & np.isfinite(values)
            if int(valid.sum()) < 2:
                continue
            lo, hi = np.percentile(values[valid], [2.0, 98.0])
            scale = float(hi - lo)
            if not np.isfinite(scale) or scale <= 0:
                continue
            normalized = np.clip((values[valid] - lo) / scale, 0.0, 1.0)
            total[valid] += normalized.astype(np.float32)
            count[valid] += 1
        prediction = np.zeros(footprint.shape, dtype=np.float32)
        valid_count = count > 0
        prediction[valid_count] = total[valid_count] / count[valid_count]
        prediction[~footprint] = 0.0  # writer replaces outside pixels with NaN

    grid.write_submission(prediction, reference, output, footprint, outside="nan")
    audit = grid.audit_submission(output, footprint, reference)
    if not audit["ok"]:
        raise SystemExit(f"Generated baseline failed format audit: {audit}")
    sidecar = {
        "schema_version": 1,
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "method": "mean of per-band 2nd/98th percentile normalized values",
        "band_indices_1_based": bands,
        "format_audit": audit,
        "warning": "Untrained baseline only; no holdout or leaderboard improvement is claimed.",
    }
    Path(str(output) + ".json").write_text(json.dumps(sidecar, indent=2, allow_nan=False) + "\n")
    print(f"output={output} sha256={sidecar['sha256']} positive={audit['positive_px']:,} format_pass={audit['ok']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
