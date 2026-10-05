"""H40 detector family -- play-fairway permeability targeting.

Design principle
----------------
The scored population is *fault pixels that are absent from the supplied
USGS/INGENIOUS catalogue*: DrivenData staff confirmed that catalogue pixels are
masked pixel-exactly from evaluation
(https://community.drivendata.org/t/11516) and that a "new fault" is any fault
pixel not already captured by USGS/INGENIOUS, including newly mapped geometry of
an existing system.  A detector for that population has to answer a different
question from "where is a lineament?" -- it has to answer "where is an
*unmapped, permeable, hot* structure?".

Each H40 channel therefore multiplies two things that prior rounds in this
project's family kept separate:

  1. STRUCTURE  -- a multi-scale, multi-physics lineament detector with an
     explicit local strike (so orientation can be used, not just amplitude).
  2. PHYSICAL FAVOURABILITY -- an independently measured reason for that
     structure to be a geothermal conduit rather than any old fracture:
     critical stress orientation (Siler 2022), heat flow (DeAngelo 2022),
     deep reservoir temperature from fluid geothermometers (GDR 1391),
     Quaternary slip rate / recency (GDR 1391), or a supervised model of
     "fault pixels that are missing from the catalogue" (USGS SGMC).

Literature anchors (manual review):
  Barton, Zoback & Moos (1995) Geology 23:913  - critically stressed faults are
      hydraulically conductive.  https://doi.org/10.1130/0091-7613(1995)023<0913:UFFASA>2.3.CO;2
  Morris, Ferrizzoli & Zoback (1996) Geology 24:1107 - slip/dilation tendency as
      a fault-permeability discriminator.  https://doi.org/10.1130/0091-7613(1996)024<1107:FPTSOT>2.3.CO;2
  Faulds & Hinz (2015) - favourable tectonic/structural settings of Great Basin
      geothermal systems.  https://www.osti.gov/servlets/purl/1724082
  Siler (2022) slip/dilation tendency release.  https://doi.org/10.5066/P9YL58W6
  DeAngelo et al. (2022) Great Basin heat flow.  https://doi.org/10.5066/P9BZPVUC
  Ayling et al. (2022) INGENIOUS / GDR 1391.    https://doi.org/10.15121/1881483
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import distance_transform_edt

from . import external as ex
from .holdout40 import STRUCT8  # noqa: F401  (shared 8-connectivity constant)
from .features import (_fill_nearest, _robust_unit, _gauss, _grad, _hessian_eigs,
                       _structure_tensor, _line_response)

R_PX = 3.0


def norm(a, foot, q_hi=0.99999, eps=1e-12):
    """Monotone, TIE-FREE normalisation for a ranking field.

    DEFECT FIXED 2026-10-05: every H40 channel ended in
    ``features._robust_unit(score, foot)``, which clips at the 99.5th percentile.
    On a 5.1 Mpx active domain that pins roughly 25 500 pixels to exactly 1.0 --
    about two thirds of a 37 654-dot budget.  The emitter then sorted a field
    whose top two thirds were numerically identical, so dot placement was
    decided by raster order instead of by geology.  Instrument hit rates of
    0.9-4.2% against the anchor's 7.65% were the symptom.

    ``norm`` divides by a very high quantile and does NOT clip, so the ordering
    is preserved everywhere and exact ties only occur where the underlying
    product is genuinely equal.  Values above 1.0 are expected and harmless:
    the emitter and the fusion both need order, not a [0,1] range.  The written
    GeoTIFF is binary {0,1} so nothing out-of-range can reach the validator.
    """
    a = np.asarray(a, np.float32)
    v = a[foot]
    if v.size == 0:
        return np.zeros(a.shape, np.float32)
    hi = float(np.quantile(v, q_hi))
    if not np.isfinite(hi) or hi <= 0:
        hi = float(np.abs(v).max()) or 1.0
    return (a / (hi + eps)).astype(np.float32)


# --------------------------------------------------------------- geometry helpers
def image_to_geographic_strike(gx, gy):
    """Convert an image-space gradient vector to a geographic *strike* axis.

    Competition grid axes: column index increases east, row index increases
    SOUTH (north-up rasters have transform[4] = -100).  So for a gradient vector
    (gx, gy) = (d/dEast, d/dSouth) the geographic azimuth of the gradient,
    measured clockwise from north, is atan2(gx, -gy).  A lineament is
    perpendicular to its gradient, so strike = azimuth + 90 deg.  Both are axes
    (mod 180), which is all a stress-orientation comparison needs.
    """
    az = np.degrees(np.arctan2(gx, -gy))
    return (az + 90.0) % 180.0


def structure_strike(fields, foot, sigmas=(1.5, 3.0)):
    """Fused multi-physics structure tensor -> (coherence, strength, strike_deg).

    The tensor is accumulated edge-energy weighted, so a strong single-physics
    edge sets the orientation where the other fields are ambiguous, and
    agreement between fields raises coherence.
    """
    Jxx = np.zeros(foot.shape, np.float64)
    Jyy = np.zeros(foot.shape, np.float64)
    Jxy = np.zeros(foot.shape, np.float64)
    strength = np.zeros(foot.shape, np.float32)
    for sig in sigmas:
        acc = np.zeros(foot.shape, np.float64)
        for f in fields:
            gx, gy = _grad(f, sig)
            m = np.hypot(gx, gy)
            Jxx += gx * gx
            Jyy += gy * gy
            Jxy += gx * gy
            acc += m
        strength += _robust_unit(acc / max(len(fields), 1), foot)
    strength = strength / len(sigmas)
    tr = Jxx + Jyy + 1e-12
    coh = np.sqrt(np.maximum((Jxx - Jyy) ** 2 / 4.0 + Jxy * Jxy, 0.0)) / tr
    # dominant gradient direction in image space
    theta = 0.5 * np.arctan2(2.0 * Jxy, Jxx - Jyy + 1e-12)
    gx = np.cos(theta) * np.sqrt(tr)
    gy = np.sin(theta) * np.sqrt(tr)
    strike = image_to_geographic_strike(gx, gy)
    return coh.astype(np.float32), strength, strike.astype(np.float32)


# ------------------------------------------------------- multi-physics edge stack
def edge_stack(bands, foot):
    """Normalised edge-magnitude stack over every physically distinct layer.

    Uses all 19 competition bands (the previous revision of this repository left
    11 of them unnamed and unused, and asked for ``tc`` at the wrong index -- see
    the defect note in ``grid.read_all_bands``).
    """
    groups = {
        "topo": ["det_elev", "det_elev_slope"],
        "mag": ["mag_anom", "rtp", "tmi", "tmi_hg", "tmi_vg", "tc"],
        "grav": ["iso_grav_anom", "iso_grav_anom_slope", "iso_grav_anom_hg", "iso_grav_anom_vg"],
        "subsurf": ["depth_to_base_surf", "cond_surf"],
        "geod": ["geod_2ndinv", "geod_shearrate", "geod_dilaterate"],
        "seis": ["deq_n100a15", "ieq_n100a15"],
    }
    out = {}
    for gname, names in groups.items():
        acc = np.zeros(foot.shape, np.float32)
        used = 0
        for nm in names:
            if nm not in bands:
                continue
            f = band(bands, nm, foot, cache)
            for sig in (1.0, 2.5):
                gx, gy = _grad(f, sig)
                acc += _robust_unit(np.hypot(gx, gy), foot)
                used += 1
        out[gname] = (acc / max(used, 1)).astype(np.float32)
    return out


STRUCTURE_BANDS = ("det_elev", "tmi", "iso_grav_anom", "cond_surf", "tc",
                   "geod_2ndinv", "depth_to_base_surf")


def band(bands, nm, foot, cache=None):
    """Band accessor that never duplicates an already-clean array.

    ``grid.read_all_bands`` nearest-fills and clips every band, so
    ``features._fill_nearest`` would only allocate a second 49 MB copy.  On a
    3.9 GB sandbox with 19 bands that copy is the difference between fitting and
    being OOM-killed.  The cache is shared with ``structure_field`` and
    ``FeatureBuilder`` so each band exists once.
    """
    if cache is None:
        cache = {}
    if nm not in cache:
        a = np.asarray(bands[nm], np.float32)
        if not bool(np.isfinite(a[foot]).all()):
            a = _fill_nearest(a, foot)
        cache[nm] = a
    return cache[nm]


def structure_field(bands, foot, ext=None, sigmas=(1.5, 3.0), cache=None):
    """Return (coherence, strength, strike_deg) for the fused physics stack.

    Memory note: this runs on a 12.3 Mpx grid in a 3.9 GB sandbox, so the
    per-band nearest-filled arrays are cached across the six H40 channels
    instead of being rebuilt by each one.
    """
    if cache is None:
        cache = {}
    fields = [band(bands, nm, foot, cache) for nm in STRUCTURE_BANDS if nm in bands]
    return structure_strike(fields, foot, sigmas=sigmas)


# ------------------------------------------------------------ H40-A: stress gate
def build_h40a(bands, foot, ext, log, cache=None):
    """Critically-stressed unmapped lineament.

    Signature targeted: orientation, not amplitude.  A lineament whose strike is
    close to the local maximum-horizontal-stress azimuth is favourably oriented
    for slip in the present stress field and is therefore far more likely to be
    hydraulically conductive (Barton et al. 1995; Morris et al. 1996).  The
    favourability curve is *measured* from the 37 811 in-footprint segments of
    the USGS Siler (2022) release -- it is not an assumed cosine -- and it is
    indexed by the release's own local SHmax azimuth, so it is location
    specific.

    Why it finds faults missing from the catalogue: the QFFD/INGENIOUS catalogue
    records faults with geomorphic expression.  Orientation favourability is a
    property of the stress field and is available everywhere, including under
    alluvium where no scarp survives.
    """
    if ext.slip_tendency is None or ext.shmax_az is None:
        log.append("H40-A: stress release unavailable -> zeros")
        return np.zeros(foot.shape, np.float32)
    # empirical curve, recomputed from the release's own segment attributes
    import csv as _csv
    p = None
    for cand in (ext_dir_candidates()):
        if cand.exists():
            p = cand
            break
    strike_l, ts_l, az_l = [], [], []
    if p is not None:
        with p.open(newline="") as fh:
            for row in _csv.DictReader(fh):
                try:
                    strike_l.append(float(row["Strike"]))
                    ts_l.append(float(row["TS_norm"]))
                    az_l.append(float(row["SHmaxAz"]))
                except (TypeError, ValueError, KeyError):
                    continue
    if len(strike_l) < 500:
        log.append("H40-A: too few stress segments -> zeros")
        return np.zeros(foot.shape, np.float32)
    centres, curve, nbins = ex.strike_favourability(
        np.array(strike_l), np.array(ts_l), np.array(az_l), n_bins=18)
    coh, strength, strike = structure_field(bands, foot, ext, cache=cache)
    dtheta = strike - (ext.shmax_az + 90.0)
    fav = ex.favourability_field(dtheta, centres, curve)
    favn = (fav - curve.min()) / max(float(curve.max() - curve.min()), 1e-9)
    score = strength * (0.25 + 0.75 * coh) * (0.30 + 0.70 * favn.astype(np.float32))
    log.append(f"H40-A: favourability curve from {len(strike_l)} segments; "
               f"TS_norm range {curve.min():.3f}-{curve.max():.3f} "
               f"(ratio {curve.max()/max(curve.min(),1e-9):.2f}x); peak at "
               f"dtheta={centres[int(np.argmax(curve))]:.1f} deg")
    return norm(score, foot)


def ext_dir_candidates():
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    return [root / "data" / "external2" / "sb_slip_tendency_in_footprint.csv",
            root / "data" / "external" / "sb_slip_tendency_in_footprint.csv"]


# ------------------------------------------------------------ H40-B: heat flow
def build_h40b(bands, foot, ext, log, cache=None):
    """Heat-flow anomaly x unmapped structure (thermal half of a play fairway).

    Signature targeted: the residual of the USGS Great Basin heat-flow field
    (DeAngelo et al. 2022, DOI 10.5066/P9BZPVUC) above its own regional trend,
    intersected with a lineament detector.  Only QC classes A/B/C of the release
    are used; class G is excluded because its distribution is a different
    population (median 218 mW/m^2 against 85 for A/B/C, maximum 10 234 mW/m^2,
    physically impossible for conductive crustal heat flow) -- verified from the
    release's own attributes, not assumed.
    """
    if ext.heat_flow is None:
        log.append("H40-B: heat-flow release unavailable -> zeros")
        return np.zeros(foot.shape, np.float32)
    hf = ext.heat_flow.astype(np.float32)
    trend = gaussian_local_mean(hf, foot, 60.0)          # ~6 km regional trend
    resid = hf - trend
    coh, strength, _ = structure_field(bands, foot, ext, cache=cache)
    score = strength * (0.25 + 0.75 * coh) * (0.35 + 0.65 * _robust_unit(resid, foot))
    log.append(f"H40-B: heat flow from {ext.heat_flow_n_wells} QC A/B/C wells; "
               f"residual vs 60 px (6 km) local mean; p99 resid="
               f"{float(np.quantile(resid[foot],0.99)):.1f} mW/m^2")
    return norm(score, foot)


def gaussian_local_mean(a, foot, sigma_px):
    v = np.where(foot, a, np.nan)
    v = np.where(np.isfinite(v), v, np.nanmean(v) if np.isfinite(v).any() else 0.0)
    return ndi.gaussian_filter(v, sigma_px, mode="nearest").astype(np.float32)


# --------------------------------------------------- H40-C: deep reservoir temp
def build_h40c(bands, foot, ext, log, cache=None):
    """Deep reservoir temperature anomaly x unmapped structure.

    Signature targeted: the maximum of the three fluid geothermometers reported
    in the INGENIOUS compilation (quartz, chalcedony, Na-K-Ca) -- an estimate of
    temperature *at reservoir depth*, not at the surface.  A hot deep reservoir
    requires a deep permeable pathway, which is a fault.  Surface heat flow
    misses this because conduction smears it; geothermometry measures the fluid
    that actually travelled the conduit.
    """
    if ext.deep_temperature is None:
        log.append("H40-C: GDR wellspring layer unavailable -> zeros")
        return np.zeros(foot.shape, np.float32)
    dt = ext.deep_temperature.astype(np.float32)
    trend = gaussian_local_mean(dt, foot, 80.0)
    resid = dt - trend
    coh, strength, _ = structure_field(bands, foot, ext, cache=cache)
    hot = ext.hot_spring if ext.hot_spring is not None else np.zeros(foot.shape, np.float32)
    score = strength * (0.25 + 0.75 * coh) * (0.35 + 0.65 * _robust_unit(resid, foot)) \
        * (0.75 + 0.50 * _robust_unit(hot, foot, lo_q=0.5, hi_q=0.999))
    log.append(f"H40-C: deep T from {ext.wells_n} GDR rows; residual vs 80 px local mean; "
               f"median deep T {float(np.median(dt[foot])):.1f} C")
    return norm(score, foot)


# ------------------------------------------- H40-D: strain-rate kink + seismic lineament
def build_h40d(bands, foot, ext, log, cache=None):
    """Geodetic strain-rate *kink* + seismicity lineament (five unused bands).

    Signature targeted: (a) the normal curvature of the second invariant of the
    strain-rate tensor -- an actively deforming zone shows as a linear *kink* in
    a smooth field, which curvature detects and amplitude thresholding does not;
    (b) elongation of the earthquake-density field (ieq_n100a15) along a
    consistent azimuth, i.e. a seismic lineament rather than a cluster;
    (c) proximity weighting by distance-to-earthquake (deq_n100a15).

    Why it finds faults missing from the catalogue: the USGS QFFD requires
    geomorphic evidence.  Strain-rate and microseismicity are measured from
    GPS/InSAR and earthquakes, so they outline faults that deform or creep
    without leaving a scarp -- precisely the blind/alluvium-covered population.
    """
    if "geod_2ndinv" not in bands or "ieq_n100a15" not in bands:
        log.append("H40-D: geodetic/seismic bands missing -> zeros")
        return np.zeros(foot.shape, np.float32)
    inv = band(bands, "geod_2ndinv", foot, cache)
    she = band(bands, "geod_shearrate", foot, cache)
    dil = band(bands, "geod_dilaterate", foot, cache)
    ieq = band(bands, "ieq_n100a15", foot, cache)
    deq = band(bands, "deq_n100a15", foot, cache)

    kink = np.zeros(foot.shape, np.float32)
    for sig in (2.0, 4.0, 8.0):
        dom, tan = _hessian_eigs(_gauss(inv, sig), sig)
        # a linear kink: large |normal curvature|, small |tangential curvature|
        k = np.abs(dom) / (np.abs(dom) + np.abs(tan) + 1e-12) * _robust_unit(np.abs(dom), foot)
        kink += _robust_unit(k, foot)
    kink /= 3.0

    # seismic lineament: coherence of the earthquake-density gradient
    coh_s, _th = _structure_tensor(np.log1p(np.maximum(ieq, 0.0)), 3.0)
    seis_line = _robust_unit(coh_s * _robust_unit(ieq, foot), foot)

    # extensional regime bonus: dilatation is extension-positive in this release
    dil_u = _robust_unit(dil, foot)
    she_u = _robust_unit(she, foot)
    prox = _robust_unit(-deq, foot)          # small distance-to-earthquake is good

    score = (0.45 * kink + 0.35 * seis_line + 0.20 * dil_u) * (0.5 + 0.5 * prox) \
        * (0.7 + 0.3 * she_u)
    log.append("H40-D: bands geod_2ndinv/geod_shearrate/geod_dilaterate/"
               "ieq_n100a15/deq_n100a15 (all unused by the previous revision)")
    return norm(score, foot)


# ------------------------------------------------- H40-E: corrected tilt angle
def build_h40e(bands, foot, ext, log, cache=None):
    """Tilt-angle zero-crossing contact locator (band 6 ``tc``).

    Signature targeted: the tilt angle arctan(VG/HG) crosses zero directly above
    a vertical contact, so its *zero crossing* is a sub-pixel contact locator,
    whereas the analytic-signal magnitude peaks off-contact.  Corroborated with
    the release's own tmi_vg / tmi_hg pair.

    DEFECT CONTEXT: the previous revision of this repository asked for band
    ``tc`` but the hardcoded name list put it at index 18, which the raster's
    own metadata shows is ``iso_grav_anom_hg``.  Every tilt-based detector built
    before this session therefore ran on the isostatic-gravity horizontal
    gradient.  This channel reads ``tc`` from the file's metadata (verified:
    band 6, description "Tilt angle or total curvature - magnetic field
    derivative for edge detection").
    """
    if "tc" not in bands:
        log.append("H40-E: tc band missing -> zeros")
        return np.zeros(foot.shape, np.float32)
    tc = band(bands, "tc", foot, cache)
    out = np.zeros(foot.shape, np.float32)
    for sig in (1.0, 2.0, 3.5):
        gx, gy = _grad(tc, sig)
        mag = np.hypot(gx, gy)
        scale = float(np.quantile(np.abs(tc[foot]), 0.90)) + 1e-9
        # zero-crossing locator: strong gradient AND small value
        zc = _robust_unit(mag, foot) * np.exp(-0.5 * (tc / scale) ** 2).astype(np.float32)
        out += zc
    out /= 3.0
    if "tmi_vg" in bands and "tmi_hg" in bands:
        vg = band(bands, "tmi_vg", foot, cache)
        hg = band(bands, "tmi_hg", foot, cache)
        ratio = np.abs(vg) / (np.abs(hg) + 1e-9)
        # tilt magnitude is monotone in vg/hg; a contact has a large ratio swing
        swing = _robust_unit(np.abs(_gauss(np.log1p(ratio), 1.5)
                                  - _gauss(np.log1p(ratio), 6.0)), foot)
        out = out * (0.6 + 0.4 * swing)
    log.append("H40-E: tc read from band metadata (band 6); zero-crossing locator "
               "with tmi_vg/tmi_hg ratio swing corroboration")
    return norm(out, foot)


# ------------------------------------- H40-F: supervised off-catalogue propensity
FEATURE_BANDS = ("mag_anom", "rtp", "tmi_hg", "geod_2ndinv", "iso_grav_anom_slope",
                 "tc", "geod_shearrate", "geod_dilaterate", "tmi_vg", "deq_n100a15",
                 "iso_grav_anom_vg", "det_elev", "iso_grav_anom", "tmi",
                 "depth_to_base_surf", "ieq_n100a15", "cond_surf", "iso_grav_anom_hg",
                 "det_elev_slope")


# One scale only.  A second smoothed scale (sigma=3 px) costs a full-grid
# Gaussian per band per prediction chunk -- 285 of them per fold, ~20 minutes of
# the run -- and the external layers (heat flow, slip tendency, deep temperature,
# slip rate) are already smoothed at 4-45 px, so the multi-scale context is
# present without recomputing it 4 times.  Dropping it also cuts the feature
# cube from 45 columns to 26.
FEATURE_SIGMAS = (0.0,)
EXTERNAL_FEATURES = ("heat_flow", "slip_tendency", "dilation_tendency",
                     "deep_temperature", "hot_spring", "slip_rate", "young_fault")


class FeatureBuilder:
    """Lazy, chunked feature extractor.

    The full cube would be 19 bands x 2 scales + 7 external layers = 45 arrays
    of 12.3 Mpx (2.2 GB) which does not fit the 3.9 GB sandbox alongside the
    rasters.  Features are therefore materialised one band at a time and only
    the requested rows/points are kept.
    """

    def __init__(self, bands, foot, ext, cache=None):
        self.bands, self.foot, self.ext = bands, foot, ext
        self.cache = cache if cache is not None else {}
        self.band_names = [nm for nm in FEATURE_BANDS if nm in bands]
        self.ext_names = [nm for nm in EXTERNAL_FEATURES if getattr(ext, nm, None) is not None]
        self.names = [f"{nm}_s{s}" for nm in self.band_names for s in FEATURE_SIGMAS] \
            + list(self.ext_names)

    def n_features(self):
        return len(self.names)

    def _base(self, nm):
        return band(self.bands, nm, self.foot, self.cache)

    def at_points(self, rows, cols):
        out = np.empty((len(self.names), rows.size), np.float32)
        i = 0
        for nm in self.band_names:
            f = self._base(nm)
            out[i] = f[rows, cols]; i += 1
            for s in FEATURE_SIGMAS[1:]:
                out[i] = _gauss(f, s)[rows, cols]; i += 1
        for nm in self.ext_names:
            out[i] = getattr(self.ext, nm)[rows, cols]; i += 1
        return out

    def band_rows(self, r0, r1):
        """(n_features, n_rows*W) float32 block for rows [r0, r1)."""
        W = self.foot.shape[1]
        out = np.empty((len(self.names), (r1 - r0) * W), np.float32)
        i = 0
        for nm in self.band_names:
            f = self._base(nm)
            out[i] = f[r0:r1].reshape(-1); i += 1
            for s in FEATURE_SIGMAS[1:]:
                out[i] = _gauss(f, s)[r0:r1].reshape(-1); i += 1
        for nm in self.ext_names:
            out[i] = np.asarray(getattr(self.ext, nm))[r0:r1].reshape(-1); i += 1
        return out


def block_labels(shape, n_blocks_y=4, n_blocks_x=6):
    """Deterministic rectangular spatial blocks, row-major id.

    DEFECT FIXED 2026-10-05: this used to default to a 3 x 4 grid while the
    evaluation in ``holdout40.blocks`` used 4 x 6.  The two grids do not nest, so
    every evaluation block lay partly inside H40-F's held-out training fold and
    partly inside folds the model *was* trained on -- roughly 75% of the SGMC-off
    positives in an evaluation block had been seen in training.  H40-F then
    "beat" the live-scored anchor by 2.05x on that instrument, which was
    leakage, not skill.  The block grid is now taken from ``holdout40.blocks``
    itself so the two cannot drift apart, and ``build_h40f`` returns the fold map
    so the alignment is asserted rather than assumed.
    """
    from .holdout40 import blocks as _blocks
    bid, _boxes, n = _blocks(shape, n_blocks_y, n_blocks_x)
    return bid.astype(np.int16), n


def build_h40f(bands, foot, ext, catalogue, log, cache=None, n_folds=4, max_pos=120_000,
               max_neg=480_000, seed=20261005, n_blocks_y=4, n_blocks_x=6,
               return_fold_map=False):
    """Supervised propensity for *fault pixels missing from the scored catalogue*.

    Positive class: USGS State Geologic Map Compilation fault pixels more than
    300 m from the scored catalogue (an official, independent realisation of the
    population the metric actually scores).  Negative class: catalogue-free
    background, stratified so that the model cannot win by learning distance to
    the catalogue.

    All predictions returned are OUT OF FOLD with respect to ``n_folds``
    rectangular spatial blocks, so evaluating this channel on the SGMC-off
    instrument is not circular.  The model is a scikit-learn
    HistGradientBoostingClassifier; no hidden or organizer label is read.

    Why this is different from prior work in this project's family: sibling
    rounds used SGMC-off only as an *evaluation* instrument.  Using it as the
    *training target* changes what is learned from "what a mapped fault looks
    like" to "what a fault the catalogue is missing looks like".
    """
    if ext.sgmc_off is None or ext.sgmc_n_px == 0:
        log.append("H40-F: SGMC off-catalogue layer unavailable -> zeros")
        z = np.zeros(foot.shape, np.float32)
        return (z, np.zeros(foot.shape, np.int8)) if return_fold_map else z
    fb = FeatureBuilder(bands, foot, ext, cache=cache)
    names = fb.names
    bid, nblocks = block_labels(foot.shape)
    active = foot & ~catalogue
    pos = active & ext.sgmc_off
    neg = active & ~ext.sgmc_off
    rng = np.random.default_rng(seed)
    # fold assignment: blocks are grouped into n_folds folds deterministically
    # fold = block id mod n_folds, so a block is predicted only by the model that
    # was trained WITHOUT that block.  With 24 blocks and 4 folds each fold holds
    # out 6 whole blocks.
    fold_of_block = np.arange(nblocks) % n_folds
    fold = fold_of_block[bid].astype(np.int8)

    oof = np.zeros(foot.shape, np.float32)
    from sklearn.ensemble import HistGradientBoostingClassifier
    for k in range(n_folds):
        tr = fold != k
        te = fold == k
        py, px = np.nonzero(pos & tr)
        ny, nx = np.nonzero(neg & tr)
        if py.size == 0:
            continue
        if py.size > max_pos:
            sel = rng.choice(py.size, max_pos, replace=False)
            py, px = py[sel], px[sel]
        # class-balance the negatives
        n_neg = min(ny.size, max_neg)
        sel = rng.choice(ny.size, n_neg, replace=False)
        ny, nx = ny[sel], nx[sel]
        Xtr = np.concatenate([fb.at_points(py, px), fb.at_points(ny, nx)], axis=1).T
        ytr = np.concatenate([np.ones(py.size, np.int8), np.zeros(ny.size, np.int8)])
        ty, tx = np.nonzero(active & te)
        if ty.size == 0:
            continue
        clf = HistGradientBoostingClassifier(
            max_iter=180, learning_rate=0.08, max_leaf_nodes=31,
            min_samples_leaf=60, l2_regularization=1.0, random_state=seed + k)
        clf.fit(Xtr, ytr)
        # Predict in ROW BANDS, not in flat point chunks.  ``at_points`` has to
        # recompute a full-array Gaussian for every chunk it is given, so
        # chunking 1.3 M points into 300 k pieces recomputed the smoothed bands
        # 5x per fold.  Row bands compute each Gaussian once per band.
        Hh, Ww = foot.shape
        band_rows = 256
        for r0 in range(0, Hh, band_rows):
            r1 = min(Hh, r0 + band_rows)
            sel = (ty >= r0) & (ty < r1)
            if not sel.any():
                continue
            yy, xx = ty[sel], tx[sel]
            oof[yy, xx] = clf.predict_proba(fb.at_points(yy, xx).T)[:, 1]
        del Xtr, ytr, clf
    log.append(f"H40-F: {int(pos.sum())} SGMC-off positives / {int(neg.sum())} background, "
               f"{len(names)} features, {n_folds} folds over {nblocks} spatial blocks "
               f"({n_blocks_y}x{n_blocks_x} -- the SAME grid holdout40 evaluates on), "
               f"out-of-fold predictions")
    out = norm(oof, foot)
    if return_fold_map:
        return out, fold
    return out


# ---------------------------------------------------------------- fusion / propensity
def fuse_play_fairway(channels, foot, log, weights=None):
    """Geometric (log-linear) fusion of the play-fairway channels.

    A geometric mean is used deliberately: a high score requires *several*
    independent factors to agree, so a single noisy channel cannot carry a pixel.
    Channels that are identically zero (unavailable external layer) are dropped
    and re-normalised rather than allowed to zero the product.
    """
    if weights is None:
        weights = {"H40-A": 1.30, "H40-B": 1.00, "H40-C": 1.10,
                   "H40-D": 0.90, "H40-E": 0.80, "H40-F": 1.60, "H40-G": 1.20}
    tot = 0.0
    acc = np.zeros(foot.shape, np.float64)
    used = []
    for k, v in channels.items():
        if v is None:
            continue
        vv = np.asarray(v, np.float64)
        if not np.any(vv[foot] > 0):
            continue
        w = weights.get(k, 1.0)
        # log of the unclipped, monotone-normalised channel: clipping at 1.0
        # would re-introduce exactly the tie mass that ``norm`` removed.
        acc += w * np.log(np.maximum(vv, 1e-6))
        tot += w
        used.append(k)
    if tot <= 0:
        return np.zeros(foot.shape, np.float32)
    fused = np.exp(acc / tot)
    log.append(f"fusion: geometric mean over {used} (weights sum {tot:.2f})")
    return norm(fused, foot)


def rank_field(a, foot):
    """Rank-normalise to [0,1] over the domain (ties broken by value)."""
    out = np.zeros(a.shape, np.float32)
    v = np.asarray(a, np.float64)[foot]
    if v.size == 0:
        return out
    order = np.argsort(v, kind="stable")
    r = np.empty(v.size, np.float64)
    r[order] = np.arange(v.size, dtype=np.float64) / max(v.size - 1, 1)
    out[foot] = r.astype(np.float32)
    return out


# ------------------------------------------ H40-G: tip-continuation / linkage corridor
def catalogue_tips_and_strikes(catalogue, foot, min_trace_px=25):
    """Skeleton tips of the supplied catalogue and the local strike at each tip.

    A tip is an endpoint of a 1-pixel-wide trace: it has exactly one 8-connected
    neighbour after thinning.  The local strike at the tip is measured from the
    trace within a 6 px window, which is what a continuation corridor has to
    follow.

    Returns (tip_mask, strike_deg_at_tip, tip_coords).
    """
    from scipy.ndimage import label as _label
    cat = catalogue & foot
    thin = _thin(cat)
    nb = np.zeros(thin.shape, np.int8)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            nb += np.roll(np.roll(thin.astype(np.int8), dy, 0), dx, 1)
    tips = thin & (nb == 1)
    comp, ncomp = _label(thin, structure=STRUCT8)
    sizes = np.bincount(comp.ravel(), minlength=ncomp + 1)
    tips &= np.isin(comp, np.nonzero(sizes >= min_trace_px)[0])
    ys, xs = np.nonzero(tips)
    strike = np.full(ys.size, np.nan, np.float32)
    for i in range(ys.size):
        y, x = int(ys[i]), int(xs[i])
        r = 6
        y0, y1 = max(0, y - r), min(thin.shape[0], y + r + 1)
        x0, x1 = max(0, x - r), min(thin.shape[1], x + r + 1)
        win = thin[y0:y1, x0:x1]
        c = comp[y0:y1, x0:x1] == comp[y, x]
        w = win & c
        if w.sum() < 4:
            continue
        yy, xx = np.nonzero(w)
        yy = yy.astype(float) - (y - y0)
        xx = xx.astype(float) - (x - x0)
        # principal direction by second moments (image space: x=east, y=south)
        sxx = float((xx * xx).mean()); syy = float((yy * yy).mean())
        sxy = float((xx * yy).mean())
        th = 0.5 * np.arctan2(2 * sxy, sxx - syy + 1e-12)
        gx, gy = np.cos(th), np.sin(th)
        strike[i] = image_to_geographic_strike(np.array([gx]), np.array([gy]))[0]
    return tips, strike, (ys, xs)


def _thin(mask, max_iter=60):
    """Morphological skeleton approximation (Zhang-Suen style via scipy).

    Implemented with binary erosion + hit-or-miss-free pruning so it stays
    deterministic and dependency-light: iteratively remove boundary pixels whose
    removal does not disconnect the local 3x3 neighbourhood.
    """
    from scipy.ndimage import binary_erosion
    skel = mask.copy()
    for _ in range(max_iter):
        er = binary_erosion(skel, structure=STRUCT8)
        boundary = skel & ~er
        if not boundary.any():
            break
        # keep a boundary pixel if it is an endpoint or if removing it would
        # disconnect its 3x3 neighbourhood
        keep = np.zeros(skel.shape, bool)
        ys, xs = np.nonzero(boundary)
        nbcount = np.zeros(skel.shape, np.int16)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy or dx:
                    nbcount += np.roll(np.roll(skel.astype(np.int16), dy, 0), dx, 1)
        endpoints = boundary & (nbcount <= 1)
        keep |= endpoints
        # simple connectivity proxy: a boundary pixel with >=3 neighbours is
        # redundant (interior of a thick trace); one with 2 neighbours on a
        # straight run is redundant too, but we keep corners
        straight = np.zeros(skel.shape, bool)
        n = skel.astype(np.int16)
        horiz = (np.roll(n, 1, 1) >= 1) & (np.roll(n, -1, 1) >= 1)
        vert = (np.roll(n, 1, 0) >= 1) & (np.roll(n, -1, 0) >= 1)
        straight = boundary & (horiz ^ vert) & (nbcount == 2)
        redundant = (boundary & ~straight) & (nbcount >= 4)
        new = skel & ~(redundant)
        if new.sum() == skel.sum():
            break
        skel = new
    return skel


def build_h40g(bands, foot, ext, catalogue, log, cache=None,
               corridor_len_px=60, corridor_half_width_px=3, min_dist_px=2.0):
    """Along-strike continuation corridor beyond catalogued fault tips.

    Signature targeted: a mapped trace that stops.  DrivenData staff confirmed
    that a "new fault" is "any fault pixel not already captured by
    USGS/INGENIOUS" and "can include newly mapped geometry of an existing fault
    system" (https://community.drivendata.org/t/11536/2), so the continuation of
    a mapped trace beyond its cartographic endpoint is squarely inside the
    scored population -- while the trace itself is masked out of it.

    Physical basis: Biasi & Wesnousky (2016, BSSA 106:1110) show that ground
    ruptures step and gap along strike with a characteristic distribution, and
    Faulds & Hinz (2015, https://www.osti.gov/servlets/purl/1724082) report that
    fault terminations and linkage/step-over settings host a large share of
    characterised Great Basin geothermal systems.  A trace that ends is either a
    real tip (where a relay ramp or splay may be unmapped) or a mapping
    limitation (where the fault simply continues).

    The corridor is emitted only BEYOND the tip, along the locally measured
    strike, at >= min_dist_px from the catalogue (so it never re-enters the
    pixel-exact mask), and it is weighted by the same structure strength and
    stress-favourability fields as H40-A so a corridor only fires where the
    geophysics actually continues the lineament.

    Difference from anything previously implemented in this repository: the
    earlier H39 family detected lineaments; none of them used the catalogue's
    own endpoint geometry to *direct* where to look, and none combined it with
    an empirically measured stress-favourability kernel.
    """
    tips, strike, (ys, xs) = catalogue_tips_and_strikes(catalogue, foot)
    if ys.size == 0:
        log.append("H40-G: no catalogue tips found -> zeros")
        return np.zeros(foot.shape, np.float32)
    coh, strength, _ = structure_field(bands, foot, ext, cache=cache)
    field = np.zeros(foot.shape, np.float32)
    H, W = foot.shape
    fav = None
    if ext.shmax_az is not None:
        p = None
        for cand in ext_dir_candidates():
            if cand.exists():
                p = cand
                break
        if p is not None:
            import csv as _csv
            sl, tl, al = [], [], []
            with p.open(newline="") as fh:
                for row in _csv.DictReader(fh):
                    try:
                        sl.append(float(row["Strike"])); tl.append(float(row["TS_norm"]))
                        al.append(float(row["SHmaxAz"]))
                    except (TypeError, ValueError, KeyError):
                        continue
            if len(sl) >= 500:
                centres, curve, _nb = ex.strike_favourability(
                    np.array(sl), np.array(tl), np.array(al), n_bins=18)
                fav = (centres, curve)
    dcat = distance_transform_edt(~(catalogue & foot))
    n_placed = 0
    for i in range(ys.size):
        y0, x0 = int(ys[i]), int(xs[i])
        s = strike[i]
        if not np.isfinite(s):
            continue
        # geographic strike -> image-space unit vector (x=east, y=south)
        az = np.deg2rad(s - 90.0)              # gradient azimuth of the trace
        ux, uy = np.cos(az), -np.sin(az)       # east, south components
        for sgn in (1.0, -1.0):
            for L in range(1, corridor_len_px + 1):
                y = int(round(y0 + sgn * uy * L))
                x = int(round(x0 + sgn * ux * L))
                if not (0 <= y < H and 0 <= x < W):
                    break
                if dcat[y, x] <= min_dist_px:
                    continue
                w = float(np.clip(1.0 - L / corridor_len_px, 0.0, 1.0)) ** 1.5
                base = 0.25 + 0.75 * float(strength[y, x]) * (0.3 + 0.7 * float(coh[y, x]))
                if fav is not None:
                    d = (float(s) - (float(ext.shmax_az[y, x]) + 90.0)) % 180.0
                    d = min(d, 180.0 - d)
                    fv = float(ex.favourability_field(np.array([d]), fav[0], fav[1])[0])
                    base *= (0.4 + 0.6 * fv)
                val = w * base
                for dy in range(-corridor_half_width_px, corridor_half_width_px + 1):
                    for dx in range(-corridor_half_width_px, corridor_half_width_px + 1):
                        rr = np.hypot(dy, dx)
                        if rr > corridor_half_width_px:
                            continue
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < H and 0 <= xx < W and foot[yy, xx] \
                                and dcat[yy, xx] > min_dist_px:
                            taper = 1.0 - rr / (corridor_half_width_px + 1.0)
                            if val * taper > field[yy, xx]:
                                field[yy, xx] = val * taper
                                n_placed += 1
    log.append(f"H40-G: {int(ys.size)} catalogue tips, {n_placed} corridor pixel writes, "
               f"corridor {corridor_len_px} px long x {corridor_half_width_px} px half-width, "
               f"all >= {min_dist_px} px from the catalogue")
    return norm(field, foot)


# ---------------------------------------------------------- H40 backbone (ridge consensus)
BACKBONE_BANDS = ("det_elev", "det_elev_slope", "tmi", "rtp", "tc",
                  "iso_grav_anom", "cond_surf", "geod_2ndinv")


def _rank_unit(v, foot):
    """Rank-normalise one surface over the domain, returning a flat vector.

    Kept separate so ``build_backbone`` can STREAM: each surface's ranks are
    accumulated into the harmonic-mean denominator and the surface is dropped,
    instead of holding 32 full-grid float32 arrays (1.6 GB) alive at once --
    which is what OOM-killed the first attempt at this function (exit 137).
    """
    vals = np.asarray(v, np.float64)[foot]
    if vals.size == 0:
        return np.zeros(0, np.float32)
    o = np.argsort(np.argsort(vals, kind="stable"), kind="stable")
    return (o / max(vals.size - 1, 1)).astype(np.float32)


def build_backbone(bands, foot, ext, log, cache=None, sigmas=(1.5, 3.0)):
    """Multi-scale lineament consensus -- the ranking backbone the H40 channels modulate.

    Rationale for having this at all: the empirically best-scoring artifacts in
    this project's family (0.2477 / 0.2600 / 0.2708 / 0.2778) are all *dotted
    ridge* emissions, i.e. their dots sit on a multi-scale lineament ranking.
    The H40 channels add physical favourability, but a favourability gate
    multiplied onto a weak backbone cannot beat a strong backbone -- measured,
    not assumed: the first H40 run scored 0.9-7.0% instrument hit rate against
    the anchor's 7.654%.  So the backbone is built first and the gates modulate
    it.

    Construction: for each of 8 physically distinct layers and 2 scales, take the
    signed-agnostic Hessian line response (ridge/valley detector) and the
    structure-tensor edge coherence times edge magnitude, then combine every
    surface by harmonic mean of ranks -- the fusion ``features.build_ensemble``
    uses, but streamed.  Harmonic mean of ranks is scale-free per surface and
    cannot be dominated by one heavy-tailed layer.

    Returns ``(backbone, n_surfaces)``.
    """
    from .features import lidar_scarp
    if cache is None:
        cache = {}
    nfoot = int(foot.sum())
    inv = np.zeros(nfoot, np.float64)
    wsum = 0.0
    nsurf = 0
    names = [nm for nm in BACKBONE_BANDS if nm in bands]
    for nm in names:
        f = band(bands, nm, foot, cache)
        for sig in sigmas:
            ridge = _line_response(f, sig, sign=0)
            rk = _rank_unit(ridge, foot)
            del ridge
            if rk.size:
                inv += 1.0 / (rk + 0.05); wsum += 1.0; nsurf += 1
            del rk
            gx, gy = _grad(f, sig)
            coh, _ = _structure_tensor(f, sig)
            edge = coh * np.hypot(gx, gy)
            del gx, gy, coh
            rk = _rank_unit(edge, foot)
            del edge
            if rk.size:
                inv += 1.0 / (rk + 0.05); wsum += 1.0; nsurf += 1
            del rk
    try:
        sc = lidar_scarp(Path(__file__).resolve().parents[2] / "data", foot)
        if sc.shape == foot.shape and float(sc.max()) > 0:
            rk = _rank_unit(_robust_unit(sc, foot), foot)
            if rk.size:
                inv += 1.0 / (rk + 0.05); wsum += 1.0; nsurf += 1
            del rk
        del sc
    except Exception as e:                       # optional corroboration, never fatal
        log.append(f"backbone: lidar_scarp unavailable ({e})")
    if wsum <= 0:
        log.append("backbone: no surfaces -> zeros")
        return np.zeros(foot.shape, np.float32), 0
    score = np.zeros(foot.shape, np.float32)
    score[foot] = (wsum / np.maximum(inv, 1e-9)).astype(np.float32)
    del inv
    log.append(f"backbone: harmonic-rank consensus of {nsurf} surfaces from "
               f"{len(names)} bands at sigmas {sigmas}")
    return norm(score, foot), nsurf
