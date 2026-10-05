"""Submission writer and auditor for GEMSDOE39.

Why this module exists instead of reusing ``gems39.grid``
---------------------------------------------------------
``gems39.grid.write_submission`` (as merged from the parallel session) raises
``ValueError`` for anything other than ``outside="nan"``.  That is the encoding
the organizer's page literally describes, and we ship it as one twin.  But the
project brief records a real portal rejection —
``Predicted values must be in range [0, 1]`` — and GEMSDOE32 measured that it has
two distinct causes:

  1. the ``-3.4028234663852886e+38`` float32 nodata sentinel written through
     into the output (7,113,308 cells in ``training_features.tif``, 3,061 of
     them inside the submission footprint);
  2. a NaN-outside raster carrying ``nodata=NaN``.

Both produce the same message.  The safest possible upload is therefore an
**all-finite** raster: real predictions inside the footprint, ``0.0`` outside,
no NaN, **no nodata tag at all**.  Every one of the 12,279,160 cells is then a
finite value in [0, 1] and there is no sentinel of any kind for a validator to
trip over.  This writer produces both encodings from one prediction and asserts
the invariants before the file is written.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio


def _template_meta(template_path: str | Path) -> dict:
    with rasterio.open(template_path) as s:
        return dict(driver="GTiff", dtype="float32", count=1, crs=s.crs,
                    transform=s.transform, width=s.width, height=s.height)


def write(prediction: np.ndarray, template_path: str | Path, out_path: str | Path,
          footprint: np.ndarray, outside: str = "zeros") -> Path:
    """Write a single-band float32 GeoTIFF that is portal-safe.

    outside="zeros" -> every cell finite; outside footprint = 0.0; nodata=None.
    outside="nan"   -> outside footprint = NaN; nodata=NaN (the page's wording).
    """
    if outside not in ("zeros", "nan"):
        raise ValueError("outside must be 'zeros' or 'nan'")
    meta = _template_meta(template_path)
    shape = (meta["height"], meta["width"])
    pred = np.asarray(prediction, dtype=np.float32)
    foot = np.asarray(footprint, bool)
    if pred.shape != shape or foot.shape != shape:
        raise ValueError(f"shape mismatch pred={pred.shape} foot={foot.shape} want={shape}")
    inside = pred[foot]
    if not np.isfinite(inside).all():
        raise ValueError("non-finite prediction inside the footprint")
    if inside.min() < 0.0 or inside.max() > 1.0:
        raise ValueError(f"prediction outside [0,1]: min={inside.min()} max={inside.max()}")
    out = np.clip(pred, 0.0, 1.0).astype(np.float32)
    if outside == "nan":
        out = np.where(foot, out, np.nan).astype(np.float32)
        meta["nodata"] = np.nan
    else:
        out = np.where(foot, out, np.float32(0.0)).astype(np.float32)
        # deliberately no nodata tag: a NaN nodata value on an all-finite raster
        # is the second documented cause of the portal's range rejection
        meta.pop("nodata", None)
    meta.update(compress="deflate", tiled=True)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **meta) as dst:
        dst.write(out, 1)
    return out_path


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def audit(path: str | Path, footprint: np.ndarray, template_path: str | Path) -> dict:
    """Re-open a written file and check every property the portal validates."""
    foot = np.asarray(footprint, bool)
    with rasterio.open(path) as s:
        a = s.read(1)
        one_band = s.count == 1
        dtype_ok = s.dtypes[0] == "float32"
        nodata = s.nodata
        s_meta = dict(crs=s.crs, width=s.width, height=s.height,
                      transform=tuple(s.transform)[:6])
        meta = dict(crs=str(s.crs), width=s.width, height=s.height,
                    transform=tuple(s.transform)[:6], dtype=s.dtypes[0],
                    nodata=("NaN" if (nodata is not None and np.isnan(nodata)) else nodata))
    if a.shape != foot.shape:
        return dict(path=str(path), ok=False, reason=f"shape {a.shape} != {foot.shape}")
    with rasterio.open(template_path) as r:
        r_meta = dict(crs=r.crs, width=r.width, height=r.height,
                      transform=tuple(r.transform)[:6])
    grid_ok = (s_meta["width"] == r_meta["width"] and s_meta["height"] == r_meta["height"]
               and s_meta["crs"] == r_meta["crs"]
               and s_meta["transform"] == r_meta["transform"])
    inside = a[foot]
    finite = int(np.isfinite(inside).sum())
    in_range = int(((inside >= 0.0) & (inside <= 1.0) & np.isfinite(inside)).sum())
    outside_vals = a[~foot]
    out_nan = int(np.isnan(outside_vals).sum())
    out_zero = int((outside_vals == 0).sum())
    out_other = int(outside_vals.size - out_nan - out_zero)
    ok = bool(one_band and dtype_ok and grid_ok
              and finite == int(foot.sum()) and in_range == int(foot.sum())
              and out_other == 0
              and (out_nan == 0 or out_zero == 0))
    return dict(path=str(path), sha256=sha256(path), bytes=Path(path).stat().st_size,
                meta=meta, grid_ok=bool(grid_ok), single_band=bool(one_band),
                dtype_float32=bool(dtype_ok),
                footprint_px=int(foot.sum()), finite_inside=finite,
                in_range_inside=in_range, outside_nan=out_nan, outside_zero=out_zero,
                outside_other=out_other, positive_px=int(((a > 0) & foot).sum()),
                min_in=float(np.nanmin(inside)), max_in=float(np.nanmax(inside)), ok=ok)
