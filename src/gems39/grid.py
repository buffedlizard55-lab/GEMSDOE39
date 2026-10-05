"""Competition-grid geometry, feature-band IO, and submission GeoTIFF handling.

The submission grid is read from the provided sample raster at runtime. Band
names are taken from each source GeoTIFF band's ``band_name`` tag or description;
we deliberately do not infer feature meaning from an undocumented numeric order.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt

NODATA_F32 = np.float32(-3.4028234663852886e38)


def open_template(path: str | Path) -> dict:
    with rasterio.open(path) as src:
        return dict(
            driver="GTiff",
            dtype="float32",
            count=1,
            crs=src.crs,
            transform=src.transform,
            width=src.width,
            height=src.height,
            nodata=np.nan,
        )


def read_footprint_from_sample(sample_path: str | Path) -> np.ndarray:
    """Return the finite-data domain from the competition's sample raster."""
    with rasterio.open(sample_path) as src:
        sample = src.read(1)
    return np.isfinite(sample)


def read_footprint(features_path: str | Path, band: int = 1) -> np.ndarray:
    """Return the valid-data domain from one feature band (fallback only)."""
    with rasterio.open(features_path) as src:
        a = src.read(band)
    return np.isfinite(a) & (a > NODATA_F32 * np.float32(0.5))


def read_labels(path: str | Path) -> np.ndarray:
    with rasterio.open(path) as src:
        return src.read(1) >= 1


def _band_names(src) -> dict[str, int]:
    """Map unambiguous names from raster metadata to 1-based band numbers."""
    names: dict[str, int] = {}
    for band in range(1, src.count + 1):
        tagged = src.tags(band).get("band_name", "").strip()
        described = (src.descriptions[band - 1] or "").strip()
        described_name = described.split(" - ", 1)[0].strip()
        if tagged and described_name and tagged.casefold() != described_name.casefold():
            raise ValueError(
                f"Band {band} name disagreement: tag={tagged!r}, description={described!r}"
            )
        name = tagged or described_name
        if not name:
            raise ValueError(
                f"Band {band} has no band_name tag or usable description; refusing to guess"
            )
        key = name.casefold()
        if key in names:
            raise ValueError(f"Duplicate feature-band name in raster: {name!r}")
        names[key] = band
    return names


def read_bands(
    features_path: str | Path,
    foot: np.ndarray,
    names: list[str] | tuple[str, ...] | None = None,
) -> dict[str, np.ndarray]:
    """Read and clean selected bands, using explicit GeoTIFF metadata names.

    Only requested bands are held in memory. Invalid in-footprint values are
    nearest-filled from valid cells, outside-footprint cells are set to zero,
    and extreme valid values are clipped to the 0.1st/99.9th percentiles.
    """
    foot = np.asarray(foot, bool)
    wanted = None if names is None else [str(n).casefold() for n in names]
    out: dict[str, np.ndarray] = {}
    with rasterio.open(features_path) as src:
        if (src.height, src.width) != foot.shape:
            raise ValueError(
                f"Feature grid {(src.height, src.width)} does not match footprint {foot.shape}"
            )
        index = _band_names(src)
        if wanted is None:
            selected = list(index.items())
        else:
            missing = [name for name in wanted if name not in index]
            if missing:
                raise ValueError(
                    f"Requested feature bands absent from {features_path}: {missing}; "
                    f"available={sorted(index)}"
                )
            selected = [(name, index[name]) for name in wanted]
        for key, band in selected:
            raw = src.read(band).astype(np.float32, copy=False)
            valid = foot & np.isfinite(raw) & (raw > NODATA_F32 * np.float32(0.999))
            if not valid.any():
                raise ValueError(f"Feature band {key!r} has no valid cells in the footprint")
            filled = np.array(raw, dtype=np.float32, copy=True)
            if not valid.all():
                nearest = distance_transform_edt(
                    ~valid, return_distances=False, return_indices=True
                )
                filled[~valid] = raw[tuple(axis[~valid] for axis in nearest)]
            values = filled[foot]
            lo, hi = np.percentile(values, [0.1, 99.9])
            if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                np.clip(filled, lo, hi, out=filled)
            filled[~foot] = np.float32(0.0)
            out[key] = filled
    return out


def read_all_bands(features_path: str | Path, foot: np.ndarray) -> dict[str, np.ndarray]:
    """Read all available features by metadata name (large; prefer ``read_bands``)."""
    return read_bands(features_path, foot, names=None)


def write_submission(
    prediction: np.ndarray,
    template_path: str | Path,
    out_path: str | Path,
    footprint: np.ndarray,
    outside: str = "nan",
) -> Path:
    """Write a single-band float32 GeoTIFF on the template grid.

    The official problem description requires null/NaN outside the provided
    bounds. Zero-filled outside variants are intentionally rejected here.
    """
    if outside != "nan":
        raise ValueError("Official submission output requires NaN outside the footprint")
    footprint = np.asarray(footprint, bool)
    meta = open_template(template_path)
    pred = np.asarray(prediction, dtype=np.float32)
    expected_shape = (meta["height"], meta["width"])
    if pred.shape != expected_shape or footprint.shape != expected_shape:
        raise ValueError(f"Expected shape {expected_shape}, got pred={pred.shape}, foot={footprint.shape}")
    inside = pred[footprint]
    if not np.isfinite(inside).all():
        raise ValueError("Non-finite predictions inside the competition footprint")
    if np.any((inside < 0.0) | (inside > 1.0)):
        raise ValueError("Predictions inside the competition footprint must be in [0, 1]")
    output = np.full(expected_shape, np.nan, dtype=np.float32)
    output[footprint] = inside
    meta.update(compress="deflate", tiled=True, nodata=np.nan)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **meta) as dst:
        dst.write(output, 1)
    return path


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_submission(path: str | Path, footprint: np.ndarray, template_path: str | Path) -> dict:
    """Re-open a written file and audit format, grid, nodata, and value range."""
    with rasterio.open(path) as src, rasterio.open(template_path) as ref:
        a = src.read(1)
        grid_ok = (
            src.count == 1
            and src.dtypes == ("float32",)
            and src.crs == ref.crs
            and src.width == ref.width
            and src.height == ref.height
            and tuple(src.transform)[:6] == tuple(ref.transform)[:6]
        )
        nodata_ok = src.nodata is not None and np.isnan(src.nodata)
        outside = a[~footprint]
        outside_ok = bool(np.isnan(outside).all())
        inside = a[footprint]
        finite = np.isfinite(inside)
        in_range = finite & (inside >= 0.0) & (inside <= 1.0)
        nodata_value = src.nodata
        if nodata_value is not None and np.isnan(nodata_value):
            nodata_value = "NaN"
        metadata = dict(
            crs=str(src.crs), width=src.width, height=src.height,
            transform=tuple(src.transform)[:6], dtype=src.dtypes[0], nodata=nodata_value,
        )
    return dict(
        path=str(path), sha256=sha256(path), bytes=Path(path).stat().st_size,
        meta=metadata, grid_ok=bool(grid_ok), nodata_nan=bool(nodata_ok),
        footprint_px=int(footprint.sum()), finite_inside=int(finite.sum()),
        in_range_inside=int(in_range.sum()), outside_nan=int(np.isnan(outside).sum()),
        positive_px=int(((a > 0) & footprint).sum()),
        min_in=float(np.nanmin(inside)), max_in=float(np.nanmax(inside)),
        ok=bool(grid_ok and nodata_ok and outside_ok and finite.all() and in_range.all()),
    )
