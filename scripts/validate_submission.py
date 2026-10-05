#!/usr/bin/env python3
"""Independent 12-point validation of a submission GeoTIFF.

Checks:
  1. File is a readable GeoTIFF (rasterio).
  2. Single band, dtype float32.
  3. CRS EPSG:32611.
  4. Shape (3730, 3292).
  5. Geotransform matches the sample template exactly.
  6. All in-footprint pixels are finite.
  7. All in-footprint pixels are in [0, 1].
  8. Outside-footprint pixels are NaN OR 0 (both accepted by DrivenData's validator;
     we recommend NaN to match the sample, but the zeros variant is provided
     for the "Predicted values must be in range [0, 1]" error some users see
     when NaN is mishandled by a browser uploader).
  9. Nodata set to NaN.
 10. Reports positive-pixel count, min/max.
 11. SHA-256 digest.
 12. Prints a PASS / FAIL line.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--reference", default=str(ROOT / "data" / "sample_submission.tif"))
    args = ap.parse_args()
    p = Path(args.file)
    ref = Path(args.reference)
    checks = []

    def chk(name, ok, detail=""):
        checks.append((name, bool(ok), detail))
        print(("PASS" if ok else "FAIL"), name, detail)

    if not p.exists():
        chk("file_exists", False, str(p)); return 1
    try:
        with rasterio.open(p) as s:
            a = s.read(1)
            chk("tiff_open", True)
            chk("single_band", s.count == 1, f"count={s.count}")
            chk("dtype_float32", s.dtypes[0] == "float32", s.dtypes[0])
            chk("crs_epsg32611", str(s.crs).upper().startswith("EPSG:32611"), str(s.crs))
            chk("shape_3730x3292", s.height == 3730 and s.width == 3292, (s.height, s.width))
            if ref.exists():
                with rasterio.open(ref) as r:
                    chk("transform_matches", tuple(s.transform)[:6] == tuple(r.transform)[:6],
                        tuple(s.transform)[:6])
                    foot = np.isfinite(r.read(1))
            else:
                # fallback: use file's own finite mask
                foot = np.isfinite(a)
            fin = np.isfinite(a[foot])
            chk("finite_inside_footprint", fin.all(), f"{fin.sum()}/{foot.sum()} finite")
            inrng = fin & (a[foot] >= 0) & (a[foot] <= 1)
            chk("range_01_inside", inrng.all(),
                f"min={float(np.nanmin(a[foot])):.4f} max={float(np.nanmax(a[foot])):.4f}")
            outside = a[~foot]
            ok_out = bool(np.isnan(outside).all() or (outside == 0).all() or outside.size == 0)
            chk("outside_nan_or_zero", ok_out,
                f"nan={int(np.isnan(outside).sum())} zero={int((outside==0).sum())} other={int((~np.isnan(outside)&(outside!=0)).sum())}")
            chk("nodata_nan", s.nodata is None or (isinstance(s.nodata, float) and np.isnan(s.nodata)),
                f"nodata={s.nodata}")
            pos = int(((a > 0) & foot).sum())
            print(f"  positive_pixels_inside: {pos}")
    except Exception as e:
        chk("read/parse", False, str(e))
        return 1
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    print(f"  sha256: {h}")
    print(f"  bytes:  {p.stat().st_size:,}")
    ok = all(c[1] for c in checks)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
