"""Pre-registered radiometric–conductivity cross-gradient for H39Y-01.

The score is a hypothesis-only detector. It uses GeoDAWN Th/K, U/K, U/Th
ratio rasters and the competition's surface-conductivity band. It does not read
labels, mapped-fault distance fields, or leaderboard outcomes.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage as ndi

from .grid import NODATA_F32, _band_names

RATIO_BANDS = ("thk", "uk", "uth")
CONDUCTIVITY_BAND = "cond_surf"
SIGMA_PX = 2.0


def _robust_unit(values: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    values = np.asarray(values, np.float32)
    valid = np.asarray(footprint, bool) & np.isfinite(values)
    if not valid.any():
        raise ValueError("Cannot scale a field with no finite footprint values")
    lo, hi = np.quantile(values[valid].astype(np.float64), [0.02, 0.995])
    out = np.zeros(values.shape, np.float32)
    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        out[valid] = np.clip((values[valid] - lo) / (hi - lo), 0.0, 1.0)
    return out


def _check_same_grid(src, *, shape, crs, transform, path: str | Path) -> None:
    if (
        (src.height, src.width) != tuple(shape)
        or src.crs != crs
        or tuple(src.transform)[:6] != tuple(transform)[:6]
    ):
        raise ValueError(f"Raster grid mismatch for {path}")


def read_ratio_bands(
    path: str | Path,
    footprint: np.ndarray,
    *,
    reference_path: str | Path,
) -> dict[str, np.ndarray]:
    """Read the three ratio layers by raster metadata, never guessed band order."""
    footprint = np.asarray(footprint, bool)
    with rasterio.open(reference_path) as ref:
        shape, crs, transform = (ref.height, ref.width), ref.crs, ref.transform
    with rasterio.open(path) as src:
        _check_same_grid(src, shape=shape, crs=crs, transform=transform, path=path)
        index = _band_names(src)
        missing = [name for name in RATIO_BANDS if name not in index]
        if missing:
            raise ValueError(f"Missing required ratio bands {missing}; found {sorted(index)}")
        if footprint.shape != (src.height, src.width):
            raise ValueError("Footprint and ratio raster shape mismatch")
        wrong_dtype = [name for name in RATIO_BANDS if src.dtypes[index[name] - 1] != "uint8"]
        if wrong_dtype:
            raise ValueError(f"Expected uint8 ratio bands, found non-uint8: {wrong_dtype}")
        out = {}
        for name in RATIO_BANDS:
            band = src.read(index[name]).astype(np.float32) / np.float32(255.0)
            if not np.isfinite(band[footprint]).all():
                raise ValueError(f"Non-finite values inside footprint in {name}")
            band[~footprint] = 0.0
            out[name] = band
    return out


def read_conductivity(
    path: str | Path,
    footprint: np.ndarray,
    *,
    reference_path: str | Path,
) -> np.ndarray:
    """Read the raw ``cond_surf`` band by metadata, without quantile clipping.

    Invalid cells inside the competition footprint are nearest-filled from
    valid in-footprint measurements; cells outside the footprint are set to
    zero. This keeps the preregistered input independent of the baseline
    feature-stack clipping transform.
    """
    footprint = np.asarray(footprint, bool)
    with rasterio.open(reference_path) as ref:
        shape, crs, transform = (ref.height, ref.width), ref.crs, ref.transform
    with rasterio.open(path) as src:
        _check_same_grid(src, shape=shape, crs=crs, transform=transform, path=path)
        index = _band_names(src)
        if CONDUCTIVITY_BAND not in index:
            raise ValueError(f"Missing {CONDUCTIVITY_BAND!r} band; found {sorted(index)}")
        if footprint.shape != (src.height, src.width):
            raise ValueError("Footprint and conductivity raster shape mismatch")
        raw = src.read(index[CONDUCTIVITY_BAND]).astype(np.float32, copy=False)
        valid = footprint & np.isfinite(raw) & (raw > NODATA_F32 * np.float32(0.999))
        if src.nodata is not None:
            if np.isnan(src.nodata):
                valid &= ~np.isnan(raw)
            else:
                valid &= raw != np.float32(src.nodata)
        if not valid.any():
            raise ValueError("No valid conductivity values inside the footprint")
        out = np.array(raw, dtype=np.float32, copy=True)
        missing_inside = footprint & ~valid
        if missing_inside.any():
            nearest = ndi.distance_transform_edt(
                ~valid, return_distances=False, return_indices=True
            )
            out[missing_inside] = raw[tuple(axis[missing_inside] for axis in nearest)]
        out[~footprint] = 0.0
        if not np.isfinite(out[footprint]).all():
            raise ValueError("Non-finite conductivity remains inside the footprint")
        return out


def radiometric_conductivity_score(
    ratios: dict[str, np.ndarray],
    conductivity: np.ndarray,
    footprint: np.ndarray,
) -> np.ndarray:
    """Compute the fixed H39Y-01 co-gradient score in [0,1].

    For each ratio and conductivity, gradients are taken after a 2-pixel
    Gaussian derivative. Ratio-gradient magnitudes are robustly scaled by the
    footprint's fixed 2nd/99.5th percentiles; the strongest ratio gradient is
    paired with the conductivity gradient. The score is sqrt(R*C) times the
    absolute gradient-direction cosine. Outside-footprint values are zero.
    """
    footprint = np.asarray(footprint, bool)
    conductivity = np.asarray(conductivity, np.float32)
    shape = footprint.shape
    if conductivity.shape != shape:
        raise ValueError("Conductivity field and footprint shape mismatch")
    if set(ratios) != set(RATIO_BANDS):
        raise ValueError(f"Expected exactly the ratio bands {RATIO_BANDS}; got {tuple(ratios)}")
    if not footprint.any():
        raise ValueError("Competition footprint is empty")
    if not np.isfinite(conductivity[footprint]).all():
        raise ValueError("Conductivity contains non-finite footprint values")

    ratio_strength = np.zeros(shape, np.float32)
    gx = np.zeros(shape, np.float32)
    gy = np.zeros(shape, np.float32)
    for name in RATIO_BANDS:
        field = np.asarray(ratios[name], np.float32)
        if field.shape != shape or not np.isfinite(field[footprint]).all():
            raise ValueError(f"Invalid ratio field {name!r}")
        ratio_x = ndi.gaussian_filter(field, sigma=SIGMA_PX, order=(0, 1), mode="nearest")
        ratio_y = ndi.gaussian_filter(field, sigma=SIGMA_PX, order=(1, 0), mode="nearest")
        strength = _robust_unit(np.hypot(ratio_x, ratio_y), footprint)
        winner = strength > ratio_strength
        np.copyto(ratio_strength, strength, where=winner)
        np.copyto(gx, ratio_x, where=winner)
        np.copyto(gy, ratio_y, where=winner)
        del field, ratio_x, ratio_y, strength, winner

    cgx = ndi.gaussian_filter(conductivity, sigma=SIGMA_PX, order=(0, 1), mode="nearest")
    cgy = ndi.gaussian_filter(conductivity, sigma=SIGMA_PX, order=(1, 0), mode="nearest")
    conductivity_strength = _robust_unit(np.hypot(cgx, cgy), footprint)
    cosine = np.abs(gx * cgx + gy * cgy) / (
        np.hypot(gx, gy) * np.hypot(cgx, cgy) + np.float32(1e-12)
    )
    np.clip(cosine, 0.0, 1.0, out=cosine)
    score = np.sqrt(ratio_strength * conductivity_strength) * cosine
    score = np.where(footprint, score, 0.0).astype(np.float32, copy=False)
    if not np.isfinite(score).all() or np.any((score < 0.0) | (score > 1.0)):
        raise ValueError("H39Y-01 score is not finite in [0, 1]")
    return score
