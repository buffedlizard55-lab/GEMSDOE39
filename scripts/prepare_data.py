#!/usr/bin/env python3
"""Inspect competition rasters and report geometry, CRS, finite value ranges.

Run after placing official files in data/ (or after scripts/restore_data.py).
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39 import grid  # noqa: E402


def main():
    ddir = ROOT / "data"
    required = ["training_features.tif", "labels.tif", "sample_submission.tif"]
    for n in required:
        p = ddir / n
        if not p.exists():
            raise SystemExit(
                f"MISSING data/{n}. Run `python scripts/restore_data.py` first "
                f"(uses gh-authenticated GitHub mirrors with sha256 integrity pins), "
                f"or place the official DrivenData download in data/."
            )
        print(f"FOUND data/{n}: {p.stat().st_size:,} bytes")
    # Feature inspection
    with rasterio.open(ddir / "training_features.tif") as s:
        print(f"training_features: shape={s.shape}, crs={s.crs}, count={s.count}, "
              f"dtype={s.dtypes[0]}, nodata={s.nodata}")
        b1 = s.read(1)
        foot = np.isfinite(b1) & (b1 > grid.NODATA_F32 * np.float32(0.5))
        print(f"  footprint (finite band1): {int(foot.sum()):,} / {b1.size:,} pixels")
    with rasterio.open(ddir / "labels.tif") as s:
        a = s.read(1)
        vals, cnts = np.unique(a, return_counts=True)
        print(f"labels: unique values/counts = {list(zip(vals.tolist(), cnts.tolist()))}")
    with rasterio.open(ddir / "sample_submission.tif") as s:
        a = s.read(1)
        print(f"sample_submission: shape={s.shape}, crs={s.crs}, nodata={s.nodata}")
        fin = np.isfinite(a)
        print(f"  finite: {int(fin.sum()):,}, nan-outside: {int(np.isnan(a[~fin]).sum()) if (~fin).any() else 0}, "
              f"min/max of finite: {float(np.nanmin(a[fin])):.4f} / {float(np.nanmax(a[fin])):.4f}")
    # Write a small JSON data receipt
    receipt = {
        "training_features": str(ddir / "training_features.tif"),
        "labels": str(ddir / "labels.tif"),
        "sample_submission": str(ddir / "sample_submission.tif"),
        "footprint_pixels": int(foot.sum()),
    }
    (ddir / "prepared_manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("OK: data/ prepared. Run `python scripts/build_pipeline.py` next.")


if __name__ == "__main__":
    main()
