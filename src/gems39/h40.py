"""GEMSDOE39 session-2 science: the measured off-catalogue fault-density profile,
the live-score-calibrated instrument, and the H40 candidate fields.

Three results in this module are new to the repository family and were computed
from bytes in ``data/`` (see ``scripts/measure_h40.py`` for the audit trail):

``rel_density_profile``
    Off-catalogue fault density as a function of distance from the competition
    catalogue, measured with USGS SGMC faults that are *not* on the catalogue as
    an independent compilation.  The profile rises to 4.4x the far-field
    baseline inside 100 m of a mapped trace and falls to 0.4x beyond 10 km.

``expected_credit`` / ``calibrated_score``
    The organizer's own metric inverted against 12 live-scored historical
    artifacts.  Per-dot expected kernel credit against the off-catalogue
    surrogate, fed through a saturating coverage model fitted to the real
    scores, ranks those 12 live scores at Spearman rho = +0.88.  Every prior
    local instrument in this repository family measured between -0.37 and +0.09.

``field_*``
    The H40 candidate fields.  ``field_sgmc_corridor`` is the only one whose
    primary evidence is a *different fault compilation* from the one the
    catalogue came from, which is what makes the catalogue-hidden holdout a
    genuine cross-source test rather than a self-fulfilling one.

Nothing here writes files.  ``scripts/run_gems40.py`` owns emission and I/O.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter

from .metric import kernel

# --------------------------------------------------------------------------- #
# context
# --------------------------------------------------------------------------- #


def load_context(data_dir, block_n=3):
    """Return the shared grids every H40 candidate and instrument uses.

    ``off`` is deliberately ``sgmc & ~cat`` (all SGMC pixels that are not exactly
    on a catalogue pixel) rather than the ``dcat > 3`` subset the previous
    revision used.  Only *exact* known-fault pixels are masked at scoring time
    (DrivenData community #11516, staff reply), so an SGMC pixel 100 m from a
    mapped trace is scoreable; dropping it silently discards the most enriched
    band of the measured profile.
    """
    import rasterio

    d = str(data_dir)
    with rasterio.open(d + "/labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    cat = lab == 1
    with rasterio.open(d + "/external/derived_sgmc_faults_100m_u8.tif") as s:
        sgmc = s.read(1) > 0
    off = sgmc & foot & ~cat
    dcat = distance_transform_edt(~cat)
    blocks = make_blocks(foot, block_n, block_n)
    return dict(foot=foot, cat=cat, sgmc=sgmc, off=off, dcat=dcat,
                blocks=blocks, shape=foot.shape)


def make_blocks(foot, n_y=3, n_x=3):
    """Contiguous spatial blocks with the ragged remainder folded into the last."""
    H, W = foot.shape
    out = np.full((H, W), -1, np.int8)
    ye = np.linspace(0, H, n_y + 1).astype(int)
    xe = np.linspace(0, W, n_x + 1).astype(int)
    for j in range(n_y):
        for i in range(n_x):
            out[ye[j]:ye[j + 1], xe[i]:xe[i + 1]] = j * n_x + i
    return np.where(foot, out, -1).astype(np.int8)


# --------------------------------------------------------------------------- #
# result 1 -- the measured relative-density profile
# --------------------------------------------------------------------------- #

#: Bin edges in *pixels* of distance-to-catalogue.  Bin 0 is the catalogue
#: itself and is excluded from the profile because those pixels are masked at
#: scoring time and can never earn credit.
PROFILE_EDGES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0,
                 13.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 300.0)


def rel_density_profile(ctx, edges=PROFILE_EDGES, exclude_px=0.0, mask=None):
    """Density of off-catalogue faults per unit area, by distance from the
    catalogue, normalised by the median band density.

    ``exclude_px`` drops the innermost bins (``d <= exclude_px``) from the
    returned profile.  The near-catalogue bins are measured but are *not* used
    for emission: they are dominated by the same physical fault re-mapped with a
    slightly different trace position, and the one controlled live-score
    experiment in this project family (removing dots within 200 m of the
    catalogue, 0.2708 -> 0.2778) says pixels there are worth less than nothing.

    ``mask`` restricts the measurement to a subset of pixels (used for the
    spatially-blocked leave-one-block-out estimate).
    """
    foot, cat, off, dcat = ctx["foot"], ctx["cat"], ctx["off"], ctx["dcat"]
    domain = foot & ~cat if mask is None else (foot & ~cat & mask)
    rows = []
    prev = 0.0
    bounds = list(edges) + [None]
    for e in bounds:
        band = domain & (dcat > prev) if e is None else domain & (dcat > prev) & (dcat <= e)
        a = int(band.sum())
        n = int((band & off).sum())
        rows.append(dict(lo_m=prev * 100.0, hi_m=(None if e is None else e * 100.0),
                         area_px=a, off_px=n, density=n / max(a, 1)))
        if e is not None:
            prev = e
    solid = [r["density"] for r in rows if r["area_px"] >= 20_000 and r["density"] > 0]
    base = float(np.median(solid)) if solid else 1.0
    for r in rows:
        r["rel"] = r["density"] / base if base > 0 else 0.0
        r["usable"] = bool(r["lo_m"] >= exclude_px * 100.0)
    return dict(bins=rows, baseline_per_px=base, exclude_px=exclude_px)


def profile_weights(prof, shape, dcat, cap=2.0, floor=0.15):
    """Turn a measured profile into a per-pixel multiplicative emission weight.

    Bins the profile marks unusable (inside the catalogue-exclusion radius) get
    weight 0.  Everything else gets ``clip(rel, floor, cap)``; the cap stops a
    single small measured bin from monopolising the whole dot budget and the
    floor stops the far field from being exactly zero, which would make the
    emission degenerate if the structural field is flat there.
    """
    w = np.zeros(shape, np.float32)
    for b in prof["bins"]:
        lo = b["lo_m"] / 100.0
        hi = np.inf if b["hi_m"] is None else b["hi_m"] / 100.0
        if b["area_px"] == 0 or not b["usable"]:
            continue
        m = (dcat > lo) & (dcat <= hi)
        w[m] = float(np.clip(b["rel"], floor, cap))
    return w


# --------------------------------------------------------------------------- #
# result 2 -- the live-score-calibrated instrument
# --------------------------------------------------------------------------- #


def expected_credit(dots, truth, valid=None):
    """Mean kernel credit one emitted dot earns against ``truth``.

    ``w = mean over emitted dots of max(0, 1 - d/R)`` with ``d`` the distance
    from the dot to the nearest truth pixel and ``R = 300 m = 3 px``.  This is
    the quantity that, fed through :func:`calibrated_score`, ranks 12
    organizer-scored artifacts at rho = +0.88.
    """
    dots = np.asarray(dots, bool)
    if valid is not None:
        dots = dots & np.asarray(valid, bool)
    n = int(dots.sum())
    if n == 0:
        return 0.0, 0
    k = kernel(distance_transform_edt(~np.asarray(truth, bool)))
    return float(k[dots].mean()), n


def calibrated_score(n_dots, w, K, a):
    """Predicted live DW-Tversky score from (dot count, per-dot credit).

    Saturating-coverage model fitted to 12 organizer-scored artifacts::

        TP  = K * (1 - exp(-a * n * w / K))
        S   = TP / (0.2*TP + 0.2*n + 0.8*K)

    ``K`` is the effective size of the scored truth chunk and ``a`` maps
    surrogate credit onto real credit.  Both are fitted, not assumed; the
    saturation is what lets the model reproduce the strong negative rank
    correlation between dot count and live score (rho = -0.86) that a
    constant-credit model cannot express at all.

    .. warning::
       Fitted on ``w`` in [0.036, 0.062].  Predictions for ``w`` far outside
       that range (an SGMC-corridor field reaches ~0.9) are extrapolation and
       are labelled as such wherever they are reported.
    """
    n_dots = np.asarray(n_dots, float)
    w = np.asarray(w, float)
    tp = K * (1.0 - np.exp(-a * n_dots * w / K))
    return tp / (0.2 * tp + 0.2 * n_dots + 0.8 * K)


def fit_calibrated(n_dots, w, scores, K_grid=None, a_grid=None):
    """Least-squares fit of ``(K, a)`` with a leave-one-out rank check."""
    from scipy.optimize import minimize
    from scipy.stats import spearmanr

    n = np.asarray(n_dots, float)
    wv = np.asarray(w, float)
    s = np.asarray(scores, float)

    def obj(th):
        p = calibrated_score(n, wv, np.exp(th[0]), np.exp(th[1]))
        return float(np.sum((p - s) ** 2))

    best = None
    for K0 in (K_grid or (2e3, 5e3, 1e4, 2e4, 5e4)):
        for a0 in (a_grid or (0.5, 2.0, 10.0, 50.0)):
            r = minimize(obj, [np.log(K0), np.log(a0)], method="Nelder-Mead",
                         options=dict(xatol=1e-10, fatol=1e-14, maxiter=40000,
                                      maxfev=40000))
            if best is None or r.fun < best.fun:
                best = r
    K, a = float(np.exp(best.x[0])), float(np.exp(best.x[1]))
    pred = calibrated_score(n, wv, K, a)
    loo = []
    for i in range(s.size):
        keep = np.ones(s.size, bool)
        keep[i] = False

        def obji(th):
            p = calibrated_score(n[keep], wv[keep], np.exp(th[0]), np.exp(th[1]))
            return float(np.sum((p - s[keep]) ** 2))

        ri = minimize(obji, best.x, method="Nelder-Mead",
                      options=dict(xatol=1e-9, fatol=1e-13, maxiter=20000, maxfev=20000))
        loo.append(float(calibrated_score(n[i], wv[i], np.exp(ri.x[0]), np.exp(ri.x[1]))[()]))
    loo = np.asarray(loo)
    return dict(K=K, a=a,
                rmse=float(np.sqrt(best.fun / s.size)),
                rho_in_sample=float(spearmanr(pred, s).statistic),
                rmse_loo=float(np.sqrt(np.mean((loo - s) ** 2))),
                rho_loo=float(spearmanr(loo, s).statistic),
                pred=pred.tolist(), loo_pred=loo.tolist(), n_points=int(s.size))


def optimal_budget(field_w_curve, K, a, budgets):
    """Pick the budget that maximises the calibrated score.

    ``field_w_curve`` maps a budget to the mean per-dot credit that budget
    achieves, measured on the surrogate.  Returns ``(budget, table)``.
    """
    table = []
    for b in budgets:
        w = field_w_curve(b)
        table.append(dict(budget=int(b), w=float(w),
                          pred=float(calibrated_score(b, w, K, a))))
    best = max(table, key=lambda r: r["pred"])
    return best["budget"], table


# --------------------------------------------------------------------------- #
# structural evidence
# --------------------------------------------------------------------------- #


def _robust_unit(a, domain, lo=0.02, hi=0.995):
    v = a[domain]
    if v.size == 0:
        return np.zeros_like(a)
    q0, q1 = np.percentile(v, [lo * 100, hi * 100])
    if not np.isfinite(q0) or not np.isfinite(q1) or q1 <= q0:
        return np.zeros_like(a)
    return np.clip((a - q0) / (q1 - q0), 0.0, 1.0).astype(np.float32)


def _read_bands(path, names, foot):
    """Read named bands, nearest-fill invalid cells inside the footprint."""
    import rasterio
    from scipy.ndimage import distance_transform_edt as _edt

    out = {}
    with rasterio.open(path) as s:
        descs = [s.descriptions[i] for i in range(s.count)]
        # band descriptions carry a human gloss after " - "; the name is the
        # token before it (e.g. "det_elev - Detrended elevation - topography ...")
        idx = {}
        for i, d in enumerate(descs):
            if d:
                idx.setdefault(d.split(" - ")[0].strip(), i + 1)
        for nm in names:
            if nm not in idx:
                raise KeyError(f"band {nm!r} not in {descs}")
            a = s.read(idx[nm]).astype(np.float32)
            ok = foot & np.isfinite(a) & (a > -3.0e38)
            filled = np.where(ok, a, np.nan).astype(np.float32)
            if (~ok & foot).any():
                # EDT of the *invalid* mask yields the index of the nearest
                # valid pixel; EDT of the valid mask would do the opposite.
                ii = _edt(~ok, return_distances=False, return_indices=True)
                bad = ~ok & foot
                filled[bad] = filled[ii[0][bad], ii[1][bad]]
            if not np.isfinite(filled[foot]).all():
                raise ValueError(f"band {nm!r}: non-finite cells survived nearest-fill")
            q0, q1 = np.percentile(filled[foot], [0.1, 99.9])
            out[nm] = np.where(foot, np.clip(filled, q0, q1), 0.0).astype(np.float32)
    return out


def _grad_mag(a, sigma):
    gy, gx = np.gradient(gaussian_filter(a, sigma))
    return np.hypot(gy, gx).astype(np.float32)


def _ridge(a, sigma):
    """Bright-ridge (Frangi-style) response: most-negative Hessian eigenvalue."""
    s = gaussian_filter(a, sigma)
    gyy = np.gradient(np.gradient(s, axis=0), axis=0)
    gxx = np.gradient(np.gradient(s, axis=1), axis=1)
    gxy = np.gradient(np.gradient(s, axis=0), axis=1)
    tr = gyy + gxx
    det = gyy * gxx - gxy * gxy
    disc = np.sqrt(np.maximum(tr * tr - 4.0 * det, 0.0))
    lam2 = 0.5 * (tr + disc)          # algebraically larger eigenvalue
    return np.maximum(-lam2, 0.0).astype(np.float32)


def _coherence(a, sigma_int=1.0, sigma_smooth=3.0):
    """Structure-tensor coherence: 1 for a perfectly linear gradient field,
    0 for an isotropic one.  Distinguishes a lineament from a blob."""
    gy, gx = np.gradient(gaussian_filter(a, sigma_int))
    jyy = gaussian_filter(gy * gy, sigma_smooth)
    jxx = gaussian_filter(gx * gx, sigma_smooth)
    jxy = gaussian_filter(gy * gx, sigma_smooth)
    tr = jyy + jxx
    disc = np.sqrt(np.maximum((jyy - jxx) ** 2 + 4.0 * jxy * jxy, 0.0))
    l1 = 0.5 * (tr + disc)
    l2 = 0.5 * (tr - disc)
    num = l1 - l2
    den = l1 + l2
    return np.where(den > 1e-12, num / np.maximum(den, 1e-12), 0.0).astype(np.float32)


def structural_field(data_dir, foot, cache=None):
    """Multi-physics corroborated lineament response.

    Deliberately different from the repository's earlier controls: the physics
    channels are combined by a *geometric* mean of three families (topography,
    magnetics, gravity+conductivity) so that a candidate must be a lineament in
    more than one physics to score highly, and each family carries a
    structure-tensor coherence term that suppresses blob-like edges.
    """
    if cache is not None and "structural" in cache:
        return cache["structural"]
    names = ["det_elev", "tmi", "rtp", "iso_grav_anom", "cond_surf", "tc"]
    B = _read_bands(str(data_dir) + "/training_features.tif", names, foot)

    topo = np.zeros(foot.shape, np.float32)
    for s_ in (1.5, 3.0):
        topo += _robust_unit(_ridge(B["det_elev"], s_), foot)
    topo /= 2.0
    topo *= 0.6 + 0.4 * _robust_unit(_coherence(B["det_elev"], 1.0, 3.0), foot)

    mag = np.zeros(foot.shape, np.float32)
    for nm in ("tmi", "rtp"):
        mag += _robust_unit(_grad_mag(B[nm], 1.5), foot)
    mag /= 2.0
    mag *= 0.6 + 0.4 * _robust_unit(_coherence(B["rtp"], 1.5, 4.0), foot)

    other = np.zeros(foot.shape, np.float32)
    for nm in ("iso_grav_anom", "cond_surf"):
        other += _robust_unit(_grad_mag(B[nm], 1.5), foot)
    other /= 2.0
    other *= 0.6 + 0.4 * _robust_unit(_ridge(B["tc"], 2.0), foot)

    eps = 1e-3
    s = ((topo + eps) * (mag + eps) * (other + eps)) ** (1.0 / 3.0)
    s = _robust_unit(s, foot)
    if cache is not None:
        cache["structural"] = s
    return s


def sgmc_corridor_field(ctx, sigma=1.0, exclude_px=2.0, mask_out=None):
    """Kernel-matched corridor around off-catalogue SGMC fault traces.

    ``mask_out`` is a boolean array of pixels whose SGMC evidence must be
    ignored (used for the leave-one-block-out cross-source test, where the
    block being scored must not contribute evidence to its own field).
    """
    off = ctx["off"] if mask_out is None else (ctx["off"] & ~mask_out)
    if not off.any():
        return np.zeros(ctx["shape"], np.float32)
    k = kernel(distance_transform_edt(~off))
    dens = gaussian_filter(k, sigma)
    f = _robust_unit(dens, ctx["foot"])
    if exclude_px > 0:
        f = np.where(ctx["dcat"] > exclude_px, f, 0.0).astype(np.float32)
    return f


def thermal_field(data_dir, ctx, sigma_px=6.0):
    """Proximity to OpenEI GDR thermal springs / wells and Quaternary volcanic
    vents, smoothed to a 600 m kernel.

    Source: OpenEI Geothermal Data Registry submission 1391, DOI
    10.15121/1881483, CC BY 4.0 -- mirrored at
    ``data/external/gdr_wellspring_in_footprint.csv`` and
    ``data/external/gdr_volcanic_vents_in_footprint.csv``.
    """
    import csv

    pts = []
    for fn in ("gdr_wellspring_in_footprint.csv", "gdr_volcanic_vents_in_footprint.csv"):
        try:
            with open(str(data_dir) + "/external/" + fn) as fh:
                for r in csv.DictReader(fh):
                    try:
                        pts.append((int(float(r["row"])), int(float(r["col"]))))
                    except (KeyError, TypeError, ValueError):
                        continue
        except OSError:
            continue
    H, W = ctx["shape"]
    m = np.zeros((H, W), bool)
    for y, x in pts:
        if 0 <= y < H and 0 <= x < W:
            m[y, x] = True
    if not m.any():
        return np.zeros((H, W), np.float32)
    d = distance_transform_edt(~m)
    return np.exp(-(d ** 2) / (2.0 * sigma_px ** 2)).astype(np.float32)


# --------------------------------------------------------------------------- #
# instrument A -- the cross-source catalogue-hidden holdout
# --------------------------------------------------------------------------- #


def hide_components(ctx, block, hide_frac=0.25, seed=0, collar=2,
                    components=None):
    """Hide whole catalogue components assigned to one spatial block.

    Returns ``(truth, visible_cat, domain)``.  ``visible_cat`` is the catalogue
    with the hidden components removed, so a candidate that uses
    distance-to-catalogue gets no free ride: as far as it is concerned the
    hidden fault was never mapped.

    Components are labelled over the **whole** catalogue and assigned to a
    block by centroid, not labelled inside the block.  Labelling inside the
    block splits any trace that crosses a block boundary, hides only the inner
    half, and leaves the outer half visible -- which is a spatial leak that
    favours exactly the local-evidence candidates this holdout is meant to
    adjudicate.  Pass ``components`` to reuse one labelling across blocks.
    """
    from scipy.ndimage import binary_dilation, label

    foot, cat = ctx["foot"], ctx["cat"]
    if components is None:
        components = label(cat, structure=np.ones((3, 3), bool))[0]
    nseg = int(components.max())
    if nseg == 0:
        return np.zeros(foot.shape, bool), cat, foot & ~cat
    ids = np.arange(1, nseg + 1)
    ys, xs = np.nonzero(cat)
    cid = components[ys, xs]
    cy = np.bincount(cid, ys, nseg + 1)[1:] / np.maximum(np.bincount(cid, None, nseg + 1)[1:], 1)
    cx = np.bincount(cid, xs, nseg + 1)[1:] / np.maximum(np.bincount(cid, None, nseg + 1)[1:], 1)
    blk = ctx["blocks"][np.clip(cy.astype(int), 0, foot.shape[0] - 1),
                        np.clip(cx.astype(int), 0, foot.shape[1] - 1)]
    sizes = np.bincount(cid, None, nseg + 1)[1:]
    inb = ids[blk == block]
    if inb.size == 0:
        return np.zeros(foot.shape, bool), cat, foot & ~cat
    rng = np.random.default_rng(seed)
    order = rng.permutation(inb)
    total = int(sizes[inb - 1].sum())
    acc, hidden = 0, []
    for i in order:
        hidden.append(int(i))
        acc += int(sizes[i - 1])
        if acc >= hide_frac * total:
            break
    hid = np.isin(components, hidden)
    truth = hid & foot
    visible = cat & ~hid
    domain = foot & ~visible & ~binary_dilation(visible, iterations=collar)
    return truth, visible, domain


def label_components(cat):
    """One 8-connected labelling of the catalogue, shared by every block."""
    from scipy.ndimage import label

    return label(np.asarray(cat, bool), structure=np.ones((3, 3), bool))[0]


def rank_lift(field, truth, domain, sample=200_000, seed=0):
    """Mean percentile rank of hidden-truth pixels inside ``domain``.

    0.5 = no better than a random pixel of the evaluation domain.  The mid-rank
    convention ``(rank_left + rank_right) / 2`` is used, not ``side="left"``:
    the left-rank convention scores every tie at the bottom of its tie group,
    so a constant field measures 0.0 instead of 0.5 and any field with repeated
    values is systematically understated.
    """
    rng = np.random.default_rng(seed)
    v = field[domain].astype(np.float32)
    t = field[truth & domain].astype(np.float32)
    if t.size == 0 or v.size < 100:
        return 0.5
    k = min(v.size, sample)
    ref = np.sort(v[rng.choice(v.size, k, replace=False)])
    lo = np.searchsorted(ref, t, side="left")
    hi = np.searchsorted(ref, t, side="right")
    return float(np.mean(0.5 * (lo + hi)) / float(k))



# --------------------------------------------------------------------------- #
# H40-E -- off-catalogue fault discriminant
# --------------------------------------------------------------------------- #

#: The 19 competition bands, in raster order, plus their official glosses as
#: read back from the GeoTIFF band descriptions (verified 2026-10-05).
BAND_ORDER = ["mag_anom", "rtp", "tmi_hg", "geod_2ndinv", "iso_grav_anom_slope",
              "tc", "geod_shearrate", "geod_dilaterate", "tmi_vg", "deq_n100a15",
              "iso_grav_anom_vg", "det_elev", "iso_grav_anom", "tmi",
              "depth_to_base_surf", "ieq_n100a15", "cond_surf",
              "iso_grav_anom_hg", "det_elev_slope"]


def _u8(a, foot):
    v = np.asarray(a, np.float32)[foot]
    if v.size == 0:
        return np.zeros(np.shape(a), np.uint8)
    q0, q1 = np.percentile(v, [0.1, 99.9])
    if not np.isfinite(q0) or not np.isfinite(q1) or q1 <= q0:
        return np.zeros(np.shape(a), np.uint8)
    return np.clip((np.asarray(a, np.float32) - q0) * (255.0 / (q1 - q0)), 0, 255
                   ).astype(np.uint8)


def feature_cube(data_dir, ctx, cache=None):
    """(H, W, C) **uint8** feature cube: the 19 supplied bands, four derived
    gradient magnitudes, an elevation-coherence channel and the catalogue-distance
    prior.

    uint8 is not a shortcut, it is a hard requirement: the float32 cube is
    1.18 GB and this project's reference box has 3 GB.  Callers convert row
    blocks to float32 at the point of use.

    ``dcat_c`` is the only *prior* channel.  It is included because the
    measured relative-density profile says distance-to-catalogue carries real
    information about off-catalogue fault presence, and letting the model
    weight it is strictly better than hard-coding a multiplier.
    """
    if cache is not None and "cube" in cache:
        return cache["cube"], cache["cube_names"]
    foot = ctx["foot"]
    B = _read_bands(str(data_dir) + "/training_features.tif", BAND_ORDER, foot)
    layers, names = [], []
    for nm in BAND_ORDER:
        layers.append(B[nm])
        names.append(nm)
    for nm in ("det_elev", "tmi", "iso_grav_anom", "cond_surf"):
        layers.append(_robust_unit(_grad_mag(B[nm], 1.5), foot))
        names.append(nm + "_gm15")
    layers.append(_robust_unit(_coherence(B["det_elev"], 1.0, 3.0), foot))
    names.append("det_elev_coh")
    layers.append(_robust_unit(np.clip(ctx["dcat"], 0.0, 50.0) / 50.0, foot))
    names.append("dcat_c")
    cube = np.stack([_u8(np.where(foot, l, 0.0), foot) for l in layers], axis=-1)
    cube[~foot] = 0
    if cache is not None:
        cache["cube"], cache["cube_names"] = cube, names
    return cube, names


def offcat_discriminant(data_dir, ctx, cache=None, n_blocks=4, seed=0,
                        iters=250, lr=0.08, leaves=31, rows=90_000,
                        hard_ring=(1, 6)):
    """Spatially-blocked out-of-fold P(pixel is an off-catalogue fault).

    Positives are off-catalogue USGS SGMC fault pixels.  Negatives are a mix of
    *hard* negatives (pixels 1-6 px from a positive, i.e. the near-miss set that
    a distance-to-fault leak would separate trivially) and background pixels.
    Each spatial block is predicted by a model fitted only on the other blocks,
    so the returned field is genuine spatial generalisation and not memorised
    location.

    This is deliberately a *generalisation* of the SGMC corridor rather than a
    restatement of it: 34.4 % of expert-mapped catalogue pixels lie within
    300 m of a geologic-map fault (3.16x the footprint base rate), so raw SGMC
    proximity can only ever address a third of the target.  The discriminant
    keeps that third and adds whatever the 19 supplied bands say about the
    other two thirds.
    """
    if cache is not None and "disc" in cache:
        return cache["disc"]
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from scipy.ndimage import binary_dilation, distance_transform_edt

    foot, cat, off = ctx["foot"], ctx["cat"], ctx["off"]
    X, names = feature_cube(data_dir, ctx, cache=cache)
    X = np.asarray(X, np.uint8)
    blocks = make_blocks(foot, n_blocks, n_blocks)
    dpos = distance_transform_edt(~off)
    hard = (dpos >= hard_ring[0]) & (dpos <= hard_ring[1]) & foot & ~off
    bg = foot & ~off & ~binary_dilation(off, iterations=hard_ring[1])
    rng = np.random.default_rng(seed)
    prob = np.zeros(foot.shape, np.float32)
    hits = np.zeros(foot.shape, np.float32)
    aucs = []
    C = X.shape[2]
    for b in range(n_blocks * n_blocks):
        test = blocks == b
        train = (blocks >= 0) & ~test & foot
        py, px = np.nonzero(off & train)
        hy, hx = np.nonzero(hard & train)
        by, bx = np.nonzero(bg & train)
        if py.size < 500 or hy.size < 500:
            continue
        m = int(min(py.size, rows))
        mh = int(min(hy.size, rows // 2))
        mb = m - mh
        pi = rng.choice(py.size, m, replace=False)
        hi = rng.choice(hy.size, mh, replace=False)
        bi = rng.choice(by.size, mb, replace=False)
        Xtr = np.empty((m + mh + mb, C), np.float32)
        Xtr[:m] = X[py[pi], px[pi]]
        Xtr[m:m + mh] = X[hy[hi], hx[hi]]
        Xtr[m + mh:] = X[by[bi], bx[bi]]
        Xtr /= 255.0
        ytr = np.concatenate([np.ones(m), np.zeros(mh + mb)])
        clf = HistGradientBoostingClassifier(
            max_iter=iters, learning_rate=lr, max_leaf_nodes=leaves,
            l2_regularization=1.0, early_stopping=False, random_state=seed)
        clf.fit(Xtr, ytr)
        del Xtr, ytr
        ys, xs = np.nonzero(test)
        order = np.argsort(ys)
        ys, xs = ys[order], xs[order]
        pred = np.zeros(ys.size, np.float32)
        r0 = 0
        while r0 < ys.size:
            r1 = int(np.searchsorted(ys, ys[r0] + 128, side="left"))
            pred[r0:r1] = clf.predict_proba(
                X[ys[r0:r1], xs[r0:r1]].astype(np.float32) / 255.0)[:, 1]
            r0 = r1
        prob[ys, xs] += pred
        hits[ys, xs] += 1.0
        ty, tx = np.nonzero(off & test)
        qy, qx = np.nonzero((~off) & test & foot)
        if ty.size > 50 and qy.size > 500:
            k = rng.choice(qy.size, min(qy.size, 60_000), replace=False)
            pm = np.zeros(foot.shape, np.float32)
            pm[ys, xs] = pred
            aucs.append(float(roc_auc_score(
                np.concatenate([np.ones(ty.size), np.zeros(k.size)]),
                np.concatenate([pm[ty, tx], pm[qy[k], qx[k]]]))))
            del pm
        del pred, clf
    prob = prob / np.maximum(hits, 1.0)
    prob = np.where(foot, prob, 0.0).astype(np.float32)
    if cache is not None:
        cache["disc"] = prob
        cache["disc_aucs"] = aucs
    return prob
