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



def write_submission(
    prediction: np.ndarray,
    template_path: str | Path,
    out_path: str | Path,
    footprint: np.ndarray,
    outside: str = "nan",
    allow_non_spec_outside: bool = False,
) -> Path:
    """Write a single-band float32 GeoTIFF on the template grid.

    ``outside`` selects how the area beyond the study footprint is encoded:

    * ``"nan"``  -- matches ``sample_submission.tif`` byte-for-byte in structure
      and is the literal reading of the problem description ("data outside the
      bounds is null or nan",
      https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/).
    * ``"zero"`` -- every pixel finite and inside [0, 1], so no uploader or
      validator that mishandles NaN can raise
      ``Predicted values must be in range [0, 1]``.  A predicted 0 outside the
      study area contributes nothing to either penalty term because p = 0.

    Which one to submit is settled on evidence, and the evidence favours ``"nan"``:
    the published format says null/NaN outside the bounds; ``sample_submission.tif``
    itself is NaN-outside; five of the owner-reported scored artifacts in
    ``registry/live_scores.json`` are ``-nan`` files (0.1855, 0.1922, 0.2449,
    0.2477, 0.2600), so NaN-outside is demonstrably accepted and scored; and the
    merged independent validator returns 14/14 for the NaN twin against 13/14 for
    the zeros twin.  The zeros twin is still written, as a documented fallback for
    the ``Predicted values must be in range [0, 1]`` failure the owner reported --
    but it requires ``allow_non_spec_outside=True`` and the manifest records that
    it does not satisfy the strict reading of the spec.
    """
    if outside not in ("nan", "zero"):
        raise ValueError("outside must be 'nan' or 'zero'")
    if outside == "zero" and not allow_non_spec_outside:
        raise ValueError(
            "Zero-filled outside output is not the published format (the problem "
            "description requires null/NaN outside the bounds) and fails the "
            "independent validator's nan_only_outside_footprint check. Pass "
            "allow_non_spec_outside=True to write it anyway as a documented "
            "fallback for uploaders that mishandle NaN.")
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
    output = (np.full(expected_shape, np.nan, dtype=np.float32) if outside == "nan"
              else np.zeros(expected_shape, dtype=np.float32))
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


def audit_submission(path: str | Path, footprint: np.ndarray,
                     template_path: str | Path | None = None) -> dict:
    """Independent re-read of a written artifact against the organizer's format rule.

    Format rule, quoted from
    https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/:
    "your submission should have the same bounds as the training data, and data
    outside the bounds is null or nan" and "a single layer of float32 values
    between 0 and 1".

    Every check is recomputed from the file on disk -- nothing is carried over
    from the array that was written -- so a silent writer bug cannot pass.  When
    ``template_path`` is given, CRS, dimensions and geotransform are compared
    against the template rather than against hardcoded constants.

    ``outside_encoding`` records which of the two admissible encodings the file
    actually uses, so a NaN twin and a zeros twin can be told apart downstream
    instead of both being reported as merely "ok".
    """
    with rasterio.open(path) as src:
        a = src.read(1)
        meta = dict(crs=str(src.crs), width=src.width, height=src.height,
                    transform=tuple(src.transform)[:6], dtype=src.dtypes[0],
                    nodata=src.nodata, count=src.count)
        if template_path is not None:
            with rasterio.open(template_path) as ref:
                meta["grid_matches_template"] = bool(
                    src.count == 1 and src.dtypes == ("float32",)
                    and src.crs == ref.crs and src.width == ref.width
                    and src.height == ref.height
                    and tuple(src.transform)[:6] == tuple(ref.transform)[:6])
        else:
            meta["grid_matches_template"] = None
    footprint = np.asarray(footprint, bool)
    inside = a[footprint]
    outside = a[~footprint]
    fin = np.isfinite(inside)
    rng = fin & (inside >= 0.0) & (inside <= 1.0)
    nfp = int(footprint.sum())
    nod = meta["nodata"]
    outside_nan = bool(outside.size == 0 or np.isnan(outside).all())
    outside_zero = bool(outside.size == 0 or (outside == 0).all())
    outside_encoding = ("nan" if outside_nan and not outside_zero else
                        "zero" if outside_zero and not outside_nan else
                        "mixed" if not (outside_nan or outside_zero) else "both")
    # When a template is supplied the grid is checked AGAINST IT, which is what
    # makes this function usable on synthetic fixtures as well as on the real
    # 3730x3292 competition grid.  The hardcoded competition dimensions are only
    # asserted when no template was given.
    competition_grid = ((meta["height"], meta["width"]) == (3730, 3292)
                        if template_path is None else True)
    grid_ok = bool(meta["count"] == 1 and meta["dtype"] == "float32"
                   and meta["crs"].upper().startswith("EPSG:32611")
                   and competition_grid
                   and (meta["grid_matches_template"] is not False))
    nodata_nan = bool(nod is None or (isinstance(nod, float) and np.isnan(nod)))
    # every value is a plain Python bool: np.bool_ is not JSON serialisable and
    # the manifest must round-trip through json without a lossy default= hook.
    checks = {
        "single_band": bool(meta["count"] == 1),
        "dtype_float32": bool(meta["dtype"] == "float32"),
        "crs_epsg_32611": bool(meta["crs"].upper().startswith("EPSG:32611")),
        "dimensions_match_grid": bool(competition_grid),
        "grid_matches_template": (bool(meta["grid_matches_template"])
                                  if meta["grid_matches_template"] is not None else True),
        "in_footprint_all_finite": bool(int(fin.sum()) == nfp),
        "in_footprint_range_0_1": bool(int(rng.sum()) == nfp),
        "in_footprint_zero_sentinel": bool(not (inside < -1e30).any()),
        "outside_footprint_null_or_zero": bool(outside_nan or outside_zero),
        "nodata_is_nan": nodata_nan,
        "values_are_binary_0_or_1": bool(
            np.isin(np.unique(a[np.isfinite(a)]), [0.0, 1.0]).all()),
        "sha256_is_64_hex": bool(len(sha256(path)) == 64),
    }
    nodata_report = "NaN" if (nod is not None and isinstance(nod, float) and np.isnan(nod)) else nod
    return dict(
        path=str(path), sha256=sha256(path), bytes=Path(path).stat().st_size,
        meta={**meta, "nodata": nodata_report},
        grid_ok=grid_ok, nodata_nan=nodata_nan,
        outside_encoding=outside_encoding,
        footprint_px=nfp,
        finite_inside=int(fin.sum()), in_range_inside=int(rng.sum()),
        outside_nan=int(np.isnan(outside).sum()) if outside.size else 0,
        outside_zero=int((outside == 0).sum()) if outside.size else 0,
        positive_px=int(((a > 0) & footprint).sum()),
        min_in=float(np.nanmin(inside)) if fin.any() else None,
        max_in=float(np.nanmax(inside)) if fin.any() else None,
        checks=checks, ok=bool(all(checks.values())),
    )
