#!/usr/bin/env python3
"""Validate a GEMS submission against the sample GeoTIFF and official format.

The competition problem page requires one float32 layer on the matching
EPSG:32611 / 100 m grid, finite probabilities in [0,1] in the footprint, and
null/NaN outside bounds. This validator uses the supplied sample raster as the
source of truth for the grid and footprint; it will not silently validate a
file against its own mask if the reference raster is missing.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", help="GeoTIFF submission to check")
    parser.add_argument(
        "--reference",
        default=str(ROOT / "data" / "sample_submission.tif"),
        help="Competition sample GeoTIFF (default: data/sample_submission.tif)",
    )
    parser.add_argument(
        "--outside",
        choices=("auto", "zeros", "nan"),
        default="auto",
        help=(
            "Outside-footprint encoding to validate. 'auto' detects it from the file. "
            "'zeros' is the all-finite encoding that avoids both documented causes of the "
            "portal's 'Predicted values must be in range [0, 1]' rejection; 'nan' is the "
            "encoding the problem page literally describes."
        ),
    )
    args = parser.parse_args()
    path = Path(args.file)
    ref_path = Path(args.reference)
    checks: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))
        print(("PASS" if ok else "FAIL"), name, detail)

    if not path.is_file():
        check("file_exists", False, str(path))
        return 1
    if not ref_path.is_file():
        check("reference_exists", False, f"{ref_path} (restore core data or pass --reference)")
        return 1

    try:
        with rasterio.open(ref_path) as ref:
            ref_meta = {
                "crs": ref.crs,
                "width": ref.width,
                "height": ref.height,
                "transform": tuple(ref.transform)[:6],
                "bounds": tuple(ref.bounds),
            }
            sample = ref.read(1)
            footprint = np.isfinite(sample)
        with rasterio.open(path) as src:
            array = src.read(1)
            check("readable_geotiff", src.driver == "GTiff", f"driver={src.driver}")
            check("single_band", src.count == 1, f"count={src.count}")
            check("float32", src.dtypes == ("float32",), f"dtypes={src.dtypes}")
            check("epsg_32611", src.crs is not None and src.crs.to_epsg() == 32611, f"crs={src.crs}")
            check("crs_matches_sample", src.crs == ref_meta["crs"], f"crs={src.crs}")
            check(
                "shape_matches_sample",
                src.width == ref_meta["width"] and src.height == ref_meta["height"],
                f"shape={(src.height, src.width)} expected={(ref_meta['height'], ref_meta['width'])}",
            )
            check(
                "transform_matches_sample",
                tuple(src.transform)[:6] == ref_meta["transform"],
                f"transform={tuple(src.transform)[:6]}",
            )
            check("bounds_match_sample", tuple(src.bounds) == ref_meta["bounds"], f"bounds={tuple(src.bounds)}")
            check(
                "resolution_100m",
                abs(src.transform.a) == 100.0 and abs(src.transform.e) == 100.0,
                f"pixel_size=({src.transform.a},{src.transform.e})",
            )
            nodata = src.nodata

        # Validate all prediction cells. NaN is allowed only outside the sample footprint.
        inside = array[footprint]
        finite_inside = np.isfinite(inside)
        in_range = finite_inside & (inside >= 0.0) & (inside <= 1.0)
        outside = array[~footprint]
        check("finite_inside_footprint", bool(finite_inside.all()), f"{int(finite_inside.sum()):,}/{int(footprint.sum()):,}")
        check(
            "range_0_to_1_inside_footprint",
            bool(in_range.all()),
            f"min={float(np.nanmin(inside)):.8g} max={float(np.nanmax(inside)):.8g}",
        )
        # Which outside-footprint encoding is this file using?
        out_nan = int(np.isnan(outside).sum())
        out_zero = int((outside == 0).sum())
        mode = args.outside
        if mode == "auto":
            mode = "nan" if out_nan > 0 else "zeros"
        print(f"  outside_footprint_encoding: {mode}")

        if mode == "nan":
            check(
                "nodata_nan",
                nodata is not None and bool(np.isnan(nodata)),
                f"nodata={nodata}",
            )
            check(
                "nan_only_outside_footprint",
                out_nan == outside.size,
                f"outside_nan={out_nan:,}/{outside.size:,}",
            )
        else:
            # All-finite encoding: the strongest possible guard against the portal's
            # range rejection. Both documented causes are excluded by construction:
            # no -3.4e38 sentinel can be written through and no NaN nodata tag exists.
            sentinel_ok = (
                nodata is None
                or (np.isfinite(nodata) and abs(float(nodata)) <= 1.0)
            )
            check("no_nodata_sentinel", bool(sentinel_ok), f"nodata={nodata}")
            check(
                "zeros_only_outside_footprint",
                out_zero == outside.size,
                f"outside_zero={out_zero:,}/{outside.size:,}",
            )
            check(
                "all_finite_everywhere",
                bool(np.isfinite(array).all()),
                f"non_finite={int((~np.isfinite(array)).sum()):,}/{array.size:,}",
            )
        check("no_infinity_anywhere", not bool(np.isinf(array).any()), f"infinities={int(np.isinf(array).sum())}")
        positive = int(np.count_nonzero((array > 0.0) & footprint))
        print(f"  positive_pixels_inside: {positive:,}")
        print(f"  sha256: {sha256_file(path)}")
        print(f"  bytes: {path.stat().st_size:,}")
    except Exception as exc:
        check("read_and_audit", False, f"{type(exc).__name__}: {exc}")
        return 1

    passed = sum(ok for _, ok, _ in checks)
    failed = len(checks) - passed
    print(f"Result: {'PASS' if failed == 0 else 'FAIL'} ({passed}/{len(checks)} checks passed)")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
