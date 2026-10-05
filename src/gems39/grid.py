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


EXPECTED_BAND_COUNT = 19

# Authoritative band order, read from the GeoTIFF's own GDAL band metadata
# (``band_name`` / ``description`` tags) on 2026-10-05 from the SHA-256 pinned
# mirror (registry/data_manifest.json id=training_features,
# sha256 4371c82e3b83...).  ``read_band_names`` re-derives it at run time and
# ``read_all_bands`` asserts the file agrees; this tuple exists only so a
# mismatch is loud instead of silent.
VERIFIED_BAND_ORDER = (
    "mag_anom",              # 1
    "rtp",                   # 2
    "tmi_hg",                # 3
    "geod_2ndinv",           # 4
    "iso_grav_anom_slope",   # 5
    "tc",                    # 6
    "geod_shearrate",        # 7
    "geod_dilaterate",       # 8
    "tmi_vg",                # 9
    "deq_n100a15",           # 10
    "iso_grav_anom_vg",      # 11
    "det_elev",              # 12
    "iso_grav_anom",         # 13
    "tmi",                   # 14
    "depth_to_base_surf",    # 15
    "ieq_n100a15",           # 16
    "cond_surf",             # 17
    "iso_grav_anom_hg",      # 18
    "det_elev_slope",        # 19
)


def read_band_names(features_path: str | Path) -> list[str]:
    """Band names straight from the raster's own metadata -- never guessed.

    Preference order: GDAL ``band_name`` metadata item, then the TIFF
    ``ImageDescription``.  Both are present in the competition file.
    """
    with rasterio.open(features_path) as s:
        out = []
        for i in range(1, s.count + 1):
            tags = s.tags(i)
            name = tags.get("band_name")
            if not name:
                desc = s.descriptions[i - 1] or ""
                name = desc.split(" - ")[0].strip()
            if not name:
                raise ValueError(f"band {i} of {features_path} carries no name; refusing to guess")
            out.append(name)
        return out


def read_all_bands(features_path: str | Path, foot: np.ndarray) -> dict[str, np.ndarray]:
    """Return clean float32 bands keyed by their authoritative metadata names.

    Nodata sentinel (-3.4e38) and non-finite pixels are replaced by
    nearest-valid fill, then zeroed outside the footprint.

    DEFECT FIXED 2026-10-05: an earlier revision of this function hardcoded a
    placeholder name list in which band 6 was unnamed and ``tc`` was placed at
    index 18.  The raster's own tags put ``tc`` at band 6 and
    ``iso_grav_anom_hg`` at band 18, so every detector that asked for ``tc``
    silently received the isostatic-gravity horizontal gradient instead.  Names
    are now read from the file and asserted against ``VERIFIED_BAND_ORDER``.
    """
    bands = {}
    names_numeric = read_band_names(features_path)
    if len(names_numeric) != EXPECTED_BAND_COUNT:
        raise ValueError(f"expected {EXPECTED_BAND_COUNT} bands, file has {len(names_numeric)}")
    if tuple(names_numeric) != VERIFIED_BAND_ORDER:
        raise ValueError(
            "band order differs from the verified order; refusing to build features on a "
            f"guessed mapping.\n  file: {names_numeric}\n  expected: {list(VERIFIED_BAND_ORDER)}")
    from scipy.ndimage import distance_transform_edt
    # Streamed one band at a time: the previous revision materialised all 19 raw
    # bands (936 MB) *and* the 19 cleaned copies at once, which together with the
    # external layers exceeded the 3.9 GB sandbox and was OOM-killed (exit 137).
    with rasterio.open(features_path) as src:
        for i, name in enumerate(names_numeric, start=1):
            a = src.read(i).astype(np.float32)
            valid = foot & np.isfinite(a) & (a > NODATA_F32 * np.float32(0.999))
            if not valid.all():
                idx = distance_transform_edt(~valid, return_distances=False,
                                             return_indices=True)
                filled = a
                filled[~valid] = a[tuple(ax[~valid] for ax in idx)]
                del idx
            else:
                filled = a
            del a
            filled[~foot] = np.float32(0.0)
            v = filled[foot]
            lo, hi = np.percentile(v, [0.1, 99.9])
            np.clip(filled, lo, hi, out=filled)
            bands[name] = filled
            del v
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
    """Independent re-read of a written artifact against the organizer's format rule.

    Format rule, quoted from
    https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/:
    "your submission should have the same bounds as the training data, and data
    outside the bounds is null or nan" and "a single layer of float32 values
    between 0 and 1".

    Every check is recomputed from the file on disk -- nothing is carried over
    from the array that was written -- so a silent writer bug cannot pass.
    """
    with rasterio.open(path) as s:
        a = s.read(1)
        meta = dict(crs=str(s.crs), width=s.width, height=s.height,
                    transform=tuple(s.transform)[:6], dtype=s.dtypes[0],
                    nodata=s.nodata, count=s.count)
    fin = np.isfinite(a[footprint])
    rng = fin & (a[footprint] >= 0) & (a[footprint] <= 1)
    outside = a[~footprint]
    nfp = int(footprint.sum())
    nod = meta["nodata"]
    # every value is a plain Python bool: np.bool_ is not JSON serialisable and
    # the manifest must round-trip through json without a lossy default= hook.
    checks = {
        "single_band": bool(meta["count"] == 1),
        "dtype_float32": bool(meta["dtype"] == "float32"),
        "crs_epsg_32611": bool(meta["crs"].upper().startswith("EPSG:32611")),
        "dimensions_3730x3292": bool((meta["height"], meta["width"]) == (3730, 3292)),
        "in_footprint_all_finite": bool(int(fin.sum()) == nfp),
        "in_footprint_range_0_1": bool(int(rng.sum()) == nfp),
        "in_footprint_zero_sentinel": bool(not (a[footprint] < -1e30).any()),
        "outside_footprint_null_or_zero": bool(
            np.isnan(outside).all() or (outside == 0).all() or outside.size == 0),
        "nodata_is_nan": bool(nod is None or (isinstance(nod, float) and np.isnan(nod))),
        "values_are_binary_0_or_1": bool(
            np.isin(np.unique(a[np.isfinite(a)]), [0.0, 1.0]).all()),
        "sha256_is_64_hex": bool(len(sha256(path)) == 64),
    }
    return dict(
        path=str(path), sha256=sha256(path), bytes=Path(path).stat().st_size, meta=meta,
        footprint_px=nfp,
        finite_inside=int(fin.sum()), in_range_inside=int(rng.sum()),
        outside_nan=int(np.isnan(outside).sum()), outside_zero=int((outside == 0).sum()),
        positive_px=int(((a > 0) & footprint).sum()),
        min_in=float(np.nanmin(a[footprint])), max_in=float(np.nanmax(a[footprint])),
        checks=checks, ok=bool(all(checks.values())),
    )
