"""Competition grid geometry, footprint, IO constants.

Verified from sample_submission.tif (sha256 pinned via registry/data_manifest.json).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio

NODATA_F32 = np.float32(-3.4028234663852886e38)


def open_template(path: str | Path) -> dict:
    with rasterio.open(path) as s:
        return dict(
            driver="GTiff", dtype="float32", count=1,
            crs=s.crs, transform=s.transform,
            width=s.width, height=s.height, nodata=np.nan,
        )


def read_footprint_from_sample(sample_path: str | Path) -> np.ndarray:
    """Boolean domain from the official sample submission's finite pixels (authoritative)."""
    with rasterio.open(sample_path) as s:
        a = s.read(1)
    return np.isfinite(a)


def read_footprint(features_path: str | Path, band: int = 1) -> np.ndarray:
    """Boolean domain from finite pixels of feature band 1."""
    with rasterio.open(features_path) as s:
        a = s.read(band)
    return np.isfinite(a) & (a > NODATA_F32 * np.float32(0.5))


def read_labels(path: str | Path) -> np.ndarray:
    with rasterio.open(path) as s:
        a = s.read(1)
    return a >= 1


def read_all_bands(features_path: str | Path, foot: np.ndarray) -> dict[str, np.ndarray]:
    """Return clean float32 bands. Nodata sentinel (-3.4e38) and non-finite pixels
    are replaced by nearest-valid fill, then zeroed outside the footprint."""
    bands = {}
    with rasterio.open(features_path) as s:
        raws = [s.read(i).astype(np.float32) for i in range(1, s.count + 1)]
    names_numeric = ["mag_anom", "rtp", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
                     "b11", "det_elev", "iso_grav_anom", "tmi", "depth_to_base_surf",
                     "b16", "cond_surf", "tc", "det_elev_slope"]
    for i, name in enumerate(names_numeric, start=1):
        a = raws[i - 1]
        valid = foot & np.isfinite(a) & (a > NODATA_F32 * np.float32(0.999))
        # nearest-fill the invalid in-footprint pixels
        if not valid.all():
            from scipy.ndimage import distance_transform_edt
            idx = distance_transform_edt(~valid, return_distances=False, return_indices=True)
            filled = a.copy()
            filled[~valid] = a[tuple(ax[~valid] for ax in idx)]
        else:
            filled = a
        filled[~foot] = np.float32(0.0)
        # Clip to p0.1-p99.9 per band inside footprint to suppress outliers
        v = filled[foot]
        lo, hi = np.percentile(v, [0.1, 99.9])
        filled = np.clip(filled, lo, hi)
        bands[name] = filled.astype(np.float32)
    return bands


def write_submission(mask: np.ndarray, template_path: str | Path, out_path: str | Path,
                     footprint: np.ndarray, outside: str = "nan") -> Path:
    """Write a single-band float32 GeoTIFF strictly in [0,1] inside footprint."""
    meta = open_template(template_path)
    out = np.asarray(mask, dtype=np.float32)
    assert out.shape == (meta["height"], meta["width"]), (out.shape, meta)
    vals = out[footprint]
    assert np.isfinite(vals).all(), "non-finite inside footprint"
    assert (vals >= -1e-6).all() and (vals <= 1.0 + 1e-6).all(), "out of [0,1] inside footprint"
    out = np.clip(out, 0.0, 1.0)
    out = np.where(footprint, out, np.nan if outside == "nan" else np.float32(0.0))
    meta["nodata"] = np.nan
    meta["compress"] = "deflate"
    meta["tiled"] = True
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


def audit_submission(path: str | Path, footprint: np.ndarray) -> dict:
    with rasterio.open(path) as s:
        a = s.read(1)
        meta = dict(crs=str(s.crs), width=s.width, height=s.height,
                    transform=tuple(s.transform)[:6], dtype=s.dtypes[0], nodata=s.nodata)
    fin = np.isfinite(a[footprint])
    rng = fin & (a[footprint] >= 0) & (a[footprint] <= 1)
    outside = a[~footprint]
    return dict(
        path=str(path), sha256=sha256(path), bytes=Path(path).stat().st_size, meta=meta,
        footprint_px=int(footprint.sum()),
        finite_inside=int(fin.sum()), in_range_inside=int(rng.sum()),
        outside_nan=int(np.isnan(outside).sum()), outside_zero=int((outside == 0).sum()),
        positive_px=int(((a > 0) & footprint).sum()),
        min_in=float(np.nanmin(a[footprint])), max_in=float(np.nanmax(a[footprint])),
        ok=bool(rng.sum() == fin.sum() == int(footprint.sum())),
    )
