#!/usr/bin/env python3
"""Inspect restored competition rasters and write a reproducible data receipt."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    data_dir = ROOT / "data"
    required = ["training_features.tif", "labels.tif", "sample_submission.tif"]
    for name in required:
        path = data_dir / name
        if not path.is_file():
            raise SystemExit(
                f"MISSING data/{name}. Run `python scripts/restore_data.py --group core` "
                "or place the authenticated DrivenData download in data/."
            )
        print(f"FOUND data/{name}: {path.stat().st_size:,} bytes; sha256={sha256_file(path)}")

    sample_path = data_dir / "sample_submission.tif"
    with rasterio.open(sample_path) as sample:
        footprint = np.isfinite(sample.read(1))
        signature = (sample.crs, sample.width, sample.height, tuple(sample.transform)[:6])
        print(
            f"sample_submission: shape={(sample.height, sample.width)}, crs={sample.crs}, "
            f"finite_footprint={int(footprint.sum()):,}/{footprint.size:,}, nodata={sample.nodata}"
        )

    feature_path = data_dir / "training_features.tif"
    with rasterio.open(feature_path) as features:
        feature_signature = (features.crs, features.width, features.height, tuple(features.transform)[:6])
        if feature_signature != signature:
            raise SystemExit("training_features.tif does not match sample_submission.tif grid")
        print(
            f"training_features: shape={(features.height, features.width)}, crs={features.crs}, "
            f"bands={features.count}, dtype={features.dtypes[0]}, nodata={features.nodata}"
        )
        band_info = []
        for index in range(1, features.count + 1):
            tags = features.tags(index)
            values = features.read(index, out_shape=(min(512, features.height), min(512, features.width)))
            valid_sample = np.isfinite(values) & (values > np.float32(-3.4028234663852886e38) * 0.5)
            band_info.append(
                {
                    "index_1_based": index,
                    "name": tags.get("band_name") or features.descriptions[index - 1],
                    "category": tags.get("data_category"),
                    "valid_sampled_values": int(valid_sample.sum()),
                    "description": tags.get("description") or features.descriptions[index - 1],
                }
            )
        print("feature band names:", [item["name"] for item in band_info])

    labels_path = data_dir / "labels.tif"
    with rasterio.open(labels_path) as labels_ds:
        labels_signature = (labels_ds.crs, labels_ds.width, labels_ds.height, tuple(labels_ds.transform)[:6])
        if labels_signature != signature:
            raise SystemExit("labels.tif does not match sample_submission.tif grid")
        label_array = labels_ds.read(1)
        values, counts = np.unique(label_array, return_counts=True)
        known_positive = int(((label_array >= 1) & footprint).sum())
        print(f"labels: unique values/counts={list(zip(values.tolist(), counts.tolist()))}")
        print(f"  positive labelled footprint cells={known_positive:,}")

    receipt = {
        "schema_version": 2,
        "observed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "provenance": "Source files were read from local data/. When restored by restore_data.py, SHA-256 pins prove owner-mirror integrity, not organizer authentication.",
        "grid": {
            "crs": str(signature[0]),
            "width": signature[1],
            "height": signature[2],
            "transform": list(signature[3]),
            "footprint_pixels": int(footprint.sum()),
        },
        "files": {
            name: {
                "bytes": (data_dir / name).stat().st_size,
                "sha256": sha256_file(data_dir / name),
            }
            for name in required
        },
        "label_values": {str(value): int(count) for value, count in zip(values, counts)},
        "positive_label_pixels_inside_footprint": known_positive,
        "feature_bands": band_info,
    }
    receipt_path = data_dir / "prepared_manifest.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(f"OK: data receipt written to {receipt_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
