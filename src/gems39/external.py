"""Free, official external layers for the H40 play-fairway family.

Every loader in this module returns arrays on the competition grid
(3730 rows x 3292 cols, EPSG:32611, 100 m, origin 243350 E / 4508550 N) and
records how each source was georeferenced, so no coordinate system is guessed.

Sources (all free / public; manual-review links):

* Siler, E.B. (2022) *Shapefile for slip tendency and dilation tendency
  calculated for Quaternary faults in the Great Basin*. USGS ScienceBase data
  release. https://doi.org/10.5066/P9YL58W6
  (item page https://www.sciencebase.gov/catalog/item/6296974dd34ec53d276bb33d)
* DeAngelo, J.V. et al. (2022) *Heat flow maps and supporting data for the
  Great Basin, USA*. USGS ScienceBase data release. https://doi.org/10.5066/P9BZPVUC
  (item page https://www.sciencebase.gov/catalog/item/6297d2fad34ec53d276c5b28)
* Ayling, B. et al. (2022) INGENIOUS Great Basin regional dataset compilation,
  DOE Geothermal Data Repository submission 1391. https://doi.org/10.15121/1881483
  (https://gdr.openei.org/submissions/1391) - wells, springs, geothermometers,
  Quaternary fault traces with slip rate and recency.
* USGS State Geologic Map Compilation (SGMC) faults, 100 m raster
  (https://mrdata.usgs.gov/geology/state/), derived by the owner's CI.

PROVENANCE CAVEAT (stated, not hidden): the byte payloads used here were
restored from the owner's integrity-pinned GitHub mirrors because this sandbox
cannot reach sciencebase.gov / gdr.openei.org / mrdata.usgs.gov (verified
2026-10-05: HTTP 000 for all four hosts).  The mirrors are SHA-256 pinned in
``registry/data_manifest.json``; they are *not* organizer-authenticated.  The
source CRS strings and attribute inventories quoted below were read from the
mirror's own schema records, which reproduce the release's ``.prj`` verbatim.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt, gaussian_filter

# Competition grid geometry, verified against data/sample_submission.tif
# (sha256 2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc).
ORIGIN_E = 243350.0
ORIGIN_N = 4508550.0
PIX_M = 100.0
WIDTH = 3292
HEIGHT = 3730

# Source CRS of the USGS Great Basin heat-flow release, copied verbatim from the
# release's USGS_gbHeatFlowWells.prj as recorded in the mirror schema
# (docs/data/sb_heat_flow_in_footprint.json -> schema.layers).  Albers Conic
# Equal Area, NAD83, centre 37.5/-117, standard parallels 29.5/45.5, no false
# easting/northing.
USAEAC_83_117 = (
    'PROJCS["USAEAC_83_117",GEOGCS["NAD83",'
    'DATUM["North_American_Datum_1983",'
    'SPHEROID["GRS 1980",6378137,298.257222101]],'
    'PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]],'
    'PROJECTION["Albers_Conic_Equal_Area"],'
    'PARAMETER["latitude_of_center",37.5],'
    'PARAMETER["longitude_of_center",-117],'
    'PARAMETER["standard_parallel_1",29.5],'
    'PARAMETER["standard_parallel_2",45.5],'
    'PARAMETER["false_easting",0],'
    'PARAMETER["false_northing",0],'
    'UNIT["metre",1]]'
)


def utm_to_rc(easting, northing):
    """EPSG:32611 metres -> (row, col) on the competition grid."""
    col = (np.asarray(easting, float) - ORIGIN_E) / PIX_M
    row = (ORIGIN_N - np.asarray(northing, float)) / PIX_M
    return row, col


def lonlat_to_rc(lon, lat, transformer=None):
    """WGS84/NAD83 lon-lat -> (row, col).  Uses pyproj (installed dependency)."""
    if transformer is None:
        transformer = _wgs_to_utm()
    e, n = transformer.transform(np.asarray(lon, float), np.asarray(lat, float))
    return utm_to_rc(e, n)


_wgs_cache = {}


def _wgs_to_utm():
    if "t" not in _wgs_cache:
        from pyproj import Transformer
        _wgs_cache["t"] = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
    return _wgs_cache["t"]


def _eac_to_utm():
    if "e" not in _wgs_cache:
        from pyproj import Transformer
        _wgs_cache["e"] = Transformer.from_crs(USAEAC_83_117, "EPSG:32611", always_xy=True)
    return _wgs_cache["e"]


def _scatter(rows, cols, values, shape, weight=None):
    """Accumulate weighted point samples into grid cells (mean per cell)."""
    out = np.zeros(shape, np.float64)
    cnt = np.zeros(shape, np.float64)
    r = np.asarray(rows).round().astype(np.int64)
    c = np.asarray(cols).round().astype(np.int64)
    v = np.asarray(values, np.float64)
    w = np.ones_like(v) if weight is None else np.asarray(weight, np.float64)
    ok = (r >= 0) & (r < shape[0]) & (c >= 0) & (c < shape[1]) & np.isfinite(v) & np.isfinite(w)
    r, c, v, w = r[ok], c[ok], v[ok], w[ok]
    flat = r * shape[1] + c
    order = np.argsort(flat)
    flat, v, w = flat[order], v[order], w[order]
    np.add.at(out.reshape(-1), flat, v * w)
    np.add.at(cnt.reshape(-1), flat, w)
    m = cnt > 0
    out[m] /= cnt[m]
    return out, cnt


def _idw(grid_vals, grid_cnt, shape, sigma_px=6.0, fill_sigma_px=25.0):
    """Smooth sparse per-cell samples into a continuous field.

    Gaussian-weighted local mean at two scales; the wide scale only fills the
    gaps the narrow scale cannot reach.  This is a smoothing interpolator, not a
    kriging estimate -- no variogram model is claimed.
    """
    num = gaussian_filter(grid_vals * grid_cnt, sigma_px, mode="nearest")
    den = gaussian_filter(grid_cnt, sigma_px, mode="nearest")
    fine = np.where(den > 1e-9, num / np.maximum(den, 1e-9), np.nan)
    num2 = gaussian_filter(np.nan_to_num(grid_vals) * grid_cnt, fill_sigma_px, mode="nearest")
    den2 = gaussian_filter(grid_cnt, fill_sigma_px, mode="nearest")
    coarse = np.where(den2 > 1e-9, num2 / np.maximum(den2, 1e-9), np.nan)
    out = np.where(np.isfinite(fine), fine, coarse)
    if not np.isfinite(out).all():
        g = float(np.nanmean(out)) if np.isfinite(out).any() else 0.0
        out = np.where(np.isfinite(out), out, g)
    return out.astype(np.float32)


@dataclass
class ExternalSet:
    """Container + provenance log for every external layer actually loaded."""
    heat_flow: np.ndarray | None = None
    heat_flow_n_wells: int = 0
    slip_tendency: np.ndarray | None = None
    dilation_tendency: np.ndarray | None = None
    shmax_az: np.ndarray | None = None
    stress_n_segments: int = 0
    deep_temperature: np.ndarray | None = None
    spring_density: np.ndarray | None = None
    hot_spring: np.ndarray | None = None
    wells_n: int = 0
    slip_rate: np.ndarray | None = None
    young_fault: np.ndarray | None = None
    qfaults_n: int = 0
    sgmc_off: np.ndarray | None = None
    sgmc_n_px: int = 0
    log: list = field(default_factory=list)

    def as_dict(self):
        return {"heat_flow_n_wells": self.heat_flow_n_wells,
                "stress_n_segments": self.stress_n_segments,
                "wells_n": self.wells_n, "qfaults_n": self.qfaults_n,
                "sgmc_off_px": self.sgmc_n_px, "log": self.log}


# ---------------------------------------------------------------- heat flow
def load_heat_flow(ddir: Path, shape, log: list, qc_ok=("A", "B", "C"),
                   hi_clip=400.0, lo_clip=20.0):
    """DeAngelo et al. 2022 (DOI 10.5066/P9BZPVUC) measured heat flow, mW/m^2.

    The release carries one row per well twice (inputs layer
    ``USGS_gbHeatFlowWells.shp`` and outputs layer
    ``USGS_gbHeatFlowWells_wEstimates.shp``); rows are de-duplicated on
    ``unique_id``.  ``qc_code`` is the release's own quality flag; ``X`` and
    ``D`` are excluded by default.  Values outside [lo_clip, hi_clip] mW/m^2
    are winsorised because the release contains non-physical outliers
    (observed max 10 234 mW/m^2 against a Great Basin regional median of
    ~135 mW/m^2).
    """
    p = ddir / "external2" / "sb_heat_flow_in_footprint.csv"
    if not p.exists():
        p = ddir / "external" / "sb_heat_flow_in_footprint.csv"
    if not p.exists():
        log.append("heat_flow: source CSV absent -> layer skipped")
        return None, 0
    seen = {}
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("kind") != "Point":
                continue
            uid = row.get("unique_id") or (row.get("x_eac117"), row.get("y_eac117"))
            if uid in seen:
                continue
            try:
                x, y = float(row["x_eac117"]), float(row["y_eac117"])
                hf = float(row["hf_meas"])
            except (TypeError, ValueError, KeyError):
                continue
            if row.get("qc_code") not in qc_ok:
                continue
            seen[uid] = (x, y, min(max(hf, lo_clip), hi_clip))
    if not seen:
        log.append("heat_flow: 0 usable QC'd wells -> layer skipped")
        return None, 0
    xs = np.array([v[0] for v in seen.values()])
    ys = np.array([v[1] for v in seen.values()])
    hs = np.array([v[2] for v in seen.values()])
    e, n = _eac_to_utm().transform(xs, ys)
    r, c = utm_to_rc(e, n)
    gv, gc = _scatter(r, c, hs, shape)
    fld = _idw(gv, gc, shape, sigma_px=8.0, fill_sigma_px=40.0)
    log.append(f"heat_flow: {len(seen)} QC{'+'.join(qc_ok)} wells, winsorised to "
               f"[{lo_clip},{hi_clip}] mW/m^2, USAEAC_83_117 -> EPSG:32611, IDW sigma 8/40 px")
    return fld, len(seen)


# ------------------------------------------------------- slip / dilation tendency
def load_stress(ddir: Path, shape, log: list, seg_min_len_m=200.0):
    """Siler (2022) slip/dilation tendency + the release's own stress model.

    Returns (TS_norm field, TD field, SHmax azimuth field, n_segments).  Each is
    the length-weighted mean over fault segments whose midpoint falls in the
    cell, smoothed.  ``SHmaxAz`` is the release's maximum-horizontal-stress
    azimuth in degrees; it is what makes the H40-A favourability kernel
    *location specific* rather than a single regional constant.
    """
    p = ddir / "external2" / "sb_slip_tendency_in_footprint.csv"
    if not p.exists():
        p = ddir / "external" / "sb_slip_tendency_in_footprint.csv"
    if not p.exists():
        log.append("stress: slip-tendency CSV absent -> layer skipped")
        return None, None, None, 0
    lon, lat, ts, td, az, ln = [], [], [], [], [], []
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                x0, y0 = float(row["X_START"]), float(row["Y_START"])
                x1, y1 = float(row["X_END"]), float(row["Y_END"])
                t_ = float(row["TS_norm"]); d_ = float(row["TD"])
                a_ = float(row["SHmaxAz"]); L = float(row.get("length_m") or 0.0)
            except (TypeError, ValueError, KeyError):
                continue
            if L < seg_min_len_m:
                continue
            lon.append(0.5 * (x0 + x1)); lat.append(0.5 * (y0 + y1))
            ts.append(t_); td.append(d_); az.append(a_); ln.append(L)
    if not lon:
        log.append("stress: 0 usable segments -> layer skipped")
        return None, None, None, 0
    r, c = lonlat_to_rc(np.array(lon), np.array(lat))
    w = np.array(ln)
    ts_f = _idw(*_scatter(r, c, np.array(ts), shape, w), shape, sigma_px=4.0, fill_sigma_px=20.0)
    td_f = _idw(*_scatter(r, c, np.array(td), shape, w), shape, sigma_px=4.0, fill_sigma_px=20.0)
    # azimuth must be averaged on the circle, not on the raw numbers
    a = np.deg2rad(np.array(az))
    azs = _idw(*_scatter(r, c, np.sin(2 * a), shape, w), shape, sigma_px=4.0, fill_sigma_px=25.0)
    azc = _idw(*_scatter(r, c, np.cos(2 * a), shape, w), shape, sigma_px=4.0, fill_sigma_px=25.0)
    az_f = (np.rad2deg(np.arctan2(azs, azc)) / 2.0) % 180.0
    log.append(f"stress: {len(lon)} segments >= {seg_min_len_m:.0f} m from "
               "GB_Quaternary_faults_TD_TS.shp (DOI 10.5066/P9YL58W6); "
               "TS_norm/TD length-weighted; SHmaxAz averaged on the circle (doubling trick)")
    return ts_f.astype(np.float32), td_f.astype(np.float32), az_f.astype(np.float32), len(lon)


def strike_favourability(strike_deg, ts_norm, shmax_az, n_bins=36):
    """Empirical P(high slip tendency | relative strike) from the Siler release.

    For each fault segment we know its strike, the local SHmax azimuth and its
    normalised slip tendency.  We bin the *relative* strike
    ``dtheta = strike - (shmax_az + 90)`` (mod 180, folded to [0, 90]) and take
    the mean TS_norm in each bin.  A normal-fault plane that is parallel to
    sigma_Hmax is favourably oriented for slip (Morris, Ferrizzoli & Zoback
    1996, Geology 24:1107; Barton, Zoback & Moos 1995, Geology 23:913), so the
    curve peaks near dtheta = 0 and falls off towards 90 degrees.  The returned
    kernel is *empirical* -- it is whatever the release says, not an assumed
    cosine.

    Returns (bin_centres_deg, mean_ts_norm, n_per_bin).
    """
    d = ((np.asarray(strike_deg, float) - (np.asarray(shmax_az, float) + 90.0)) % 180.0)
    d = np.minimum(d, 180.0 - d)          # fold to [0, 90]
    ts = np.asarray(ts_norm, float)
    ok = np.isfinite(d) & np.isfinite(ts)
    d, ts = d[ok], ts[ok]
    edges = np.linspace(0.0, 90.0, n_bins + 1)
    idx = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, n_bins - 1)
    s = np.bincount(idx, weights=ts, minlength=n_bins)
    n = np.bincount(idx, minlength=n_bins)
    mean = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    glob = float(np.nanmean(mean)) if np.isfinite(mean).any() else 0.5
    mean = np.where(np.isfinite(mean), mean, glob)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, mean.astype(np.float64), n


def favourability_field(rel_strike_deg, centres, curve):
    """Look up the empirical favourability curve at each pixel's relative strike."""
    rel = np.asarray(rel_strike_deg, float)
    rel = np.minimum(rel % 180.0, 180.0 - (rel % 180.0))
    step = centres[1] - centres[0]
    idx = np.clip(((rel - centres[0]) / step + 0.5).astype(np.int64), 0, len(curve) - 1)
    return curve[idx].astype(np.float32)


# ------------------------------------------------------------- INGENIOUS wells
def load_wellspring(ddir: Path, shape, log: list):
    """GDR 1391 wells / springs (DOI 10.15121/1881483).

    The mirror already carries ``row``/``col`` in competition-grid coordinates.
    That mapping is *verified* here rather than trusted: the same rows also
    carry ``utm_x``/``utm_y``, so we recompute row/col from the UTM values with
    the grid transform and assert agreement.  A disagreement raises instead of
    silently shifting the layer.

    Returns (deep_temperature field, hot-spring/well point field, n_rows).
    Deep temperature is the maximum of the three reported geothermometers
    (quartz, chalcedony, Na-K-Ca) per location -- a reservoir-temperature
    estimate, i.e. the heat half of a play fairway at depth rather than at the
    surface.
    """
    p = ddir / "external" / "gdr_wellspring_in_footprint.csv"
    if not p.exists():
        log.append("wellspring: CSV absent -> layer skipped")
        return None, None, 0
    rows, cols, deep, hot = [], [], [], []
    mism = 0
    checked = 0
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                r0 = float(row["row"]); c0 = float(row["col"])
            except (TypeError, ValueError, KeyError):
                continue
            try:
                ux, uy = float(row["utm_x"]), float(row["utm_y"])
                rr, cc = utm_to_rc(ux, uy)
                checked += 1
                if abs(rr - r0) > 1.5 or abs(cc - c0) > 1.5:
                    mism += 1
            except (TypeError, ValueError):
                pass
            g = []
            for k in ("geothermquartz_c", "geothermchalc_c", "geothermcat_c"):
                try:
                    v = float(row[k])
                    if 20.0 <= v <= 400.0:
                        g.append(v)
                except (TypeError, ValueError, KeyError):
                    pass
            tcls = (row.get("thermalclass") or "").strip().lower()
            try:
                t = float(row["temp_c"])
            except (TypeError, ValueError, KeyError):
                t = float("nan")
            rows.append(r0); cols.append(c0)
            deep.append(max(g) if g else float("nan"))
            hot.append(1.0 if (tcls.startswith("hot") or (np.isfinite(t) and t >= 65.0)) else 0.0)
    if mism > max(5, 0.01 * checked):
        raise ValueError(f"wellspring row/col does not match utm_x/utm_y for {mism}/{checked} rows; "
                         "refusing to place the layer on a guessed transform")
    r = np.array(rows); c = np.array(cols)
    dv = np.array(deep, float); hv = np.array(hot, float)
    dgrid, dcnt = _scatter(r, c, np.nan_to_num(dv, nan=0.0), shape, np.isfinite(dv).astype(float))
    deep_f = _idw(dgrid, dcnt, shape, sigma_px=8.0, fill_sigma_px=45.0)
    hgrid, hcnt = _scatter(r, c, hv, shape)
    hot_f = _idw(hgrid, hcnt, shape, sigma_px=6.0, fill_sigma_px=30.0)
    log.append(f"wellspring: {len(rows)} rows ({int(np.isfinite(dv).sum())} with a geothermometer, "
               f"{int(hv.sum())} hot); row/col cross-checked against utm_x/utm_y on {checked} rows, "
               f"{mism} mismatches; IDW sigma 8/45 px (deep T) and 6/30 px (hot)")
    return deep_f.astype(np.float32), hot_f.astype(np.float32), len(rows)


# ------------------------------------------------------- Qfaults slip rate/recency
def load_qfaults(ddir: Path, shape, log: list):
    """GDR 1391 Quaternary fault traces: slip rate and recency.

    Faulds & Hinz (2015) report that characterised Great Basin geothermal
    systems sit preferentially on young, actively slipping structures
    (https://www.osti.gov/servlets/purl/1724082).  We build two fields: a
    length-weighted slip-rate field and a 'young fault' field (recency
    < 130 ka, i.e. the release's youngest two classes).
    """
    p = ddir / "external" / "gdr_qfaults_traces.csv"
    if not p.exists():
        log.append("qfaults: CSV absent -> layer skipped")
        return None, None, 0
    r, c, sr, yng, ln = [], [], [], [], []
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                r0 = float(row["centroid_row"]); c0 = float(row["centroid_col"])
            except (TypeError, ValueError, KeyError):
                continue
            if str(row.get("centroid_in_footprint", "1")) not in ("1", "1.0", "True", "true"):
                continue
            try:
                s = float(row["slip_rate"])
            except (TypeError, ValueError):
                s = float("nan")
            rec = (row.get("recency") or "").replace(",", "").replace(" ", "")
            y = 1.0 if rec in ("<15000", "<130000", "<13000", "<150") else 0.0
            try:
                L = float(row.get("clipped_length_m") or row.get("full_length_m") or 0.0)
            except (TypeError, ValueError):
                L = 0.0
            r.append(r0); c.append(c0); sr.append(s); yng.append(y); ln.append(max(L, 1.0))
    if not r:
        log.append("qfaults: 0 in-footprint traces -> layer skipped")
        return None, None, 0
    w = np.array(ln)
    srf = _idw(*_scatter(np.array(r), np.array(c), np.nan_to_num(np.array(sr)), shape,
                            w * np.isfinite(np.array(sr)).astype(float)),
                  shape, sigma_px=8.0, fill_sigma_px=40.0)
    ygf = _idw(*_scatter(np.array(r), np.array(c), np.array(yng), shape, w),
                  shape, sigma_px=8.0, fill_sigma_px=40.0)
    log.append(f"qfaults: {len(r)} in-footprint trace centroids; slip rate length-weighted; "
               f"young = recency < 130 ka ({int(np.array(yng).sum())} traces)")
    return srf.astype(np.float32), ygf.astype(np.float32), len(r)


# ------------------------------------------------------------------ SGMC off-cat
def sgmc_off_catalogue(ddir: Path, catalogue, foot, min_dist_px=3.0):
    """USGS SGMC fault pixels that are NOT in the scored catalogue.

    DrivenData staff confirmed that pixels of the supplied USGS/INGENIOUS
    catalogue are masked *pixel-exactly* from scoring
    (https://community.drivendata.org/t/11516), so the scored population is
    'fault pixels absent from the catalogue'.  SGMC faults > min_dist_px from
    the catalogue are an independent, official realisation of exactly that
    population, which is why they are used both as a holdout instrument and --
    new in H40 -- as the positive class of a supervised propensity model.
    """
    p = ddir / "external" / "derived_sgmc_faults_100m_u8.tif"
    if not p.exists():
        return None, 0
    with rasterio.open(p) as s:
        sg = (s.read(1) > 0) & foot
    d = distance_transform_edt(~catalogue)
    off = sg & (d > min_dist_px)
    return off, int(off.sum())


def load_all(ddir: Path, shape, catalogue=None, foot=None):
    log: list = []
    ext = ExternalSet(log=log)
    ext.heat_flow, ext.heat_flow_n_wells = load_heat_flow(ddir, shape, log)
    ext.slip_tendency, ext.dilation_tendency, ext.shmax_az, ext.stress_n_segments = \
        load_stress(ddir, shape, log)
    ext.deep_temperature, ext.hot_spring, ext.wells_n = load_wellspring(ddir, shape, log)
    ext.slip_rate, ext.young_fault, ext.qfaults_n = load_qfaults(ddir, shape, log)
    if catalogue is not None and foot is not None:
        ext.sgmc_off, ext.sgmc_n_px = sgmc_off_catalogue(ddir, catalogue, foot)
        log.append(f"sgmc_off_catalogue: {ext.sgmc_n_px} px more than 300 m from the scored catalogue")
    return ext
