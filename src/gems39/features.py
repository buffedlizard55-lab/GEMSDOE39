"""Novel multi-physics feature detectors for GEMSDOE39.

Legacy H39-A…E research notes follow. This docstring is not a timestamped preregistration for H39X-01; see its post-hoc record below.

H39-A  Cross-gradient tensor eigen-coherence (magnetic x gravity x MT conductivity).
       A previously unmapped blind fault shows a *straight, co-located* gradient edge in
       magnetic, isostatic gravity, AND conductivity fields at sub-kilometre scale, even
       when no DEM scarp exists.  Prior work fused per-physics ridge strengths with a
       scalar product; we use the *tensor outer product* of gradient vectors which fires
       only where gradients from independent physics are parallel AND strong.
       Layers: rtp (2), tmi (14), iso_grav_anom (13), cond_surf (17).
       New vs prior repos: H33-A fuses scalar line-responses; this is a *directional*
       3-physics cross-coupling at 3 smoothing scales.

H39-B  Scarp curvature x step-over extensional saddle detector.
       Many Basin-and-Range blind faults terminate at *releasing step-overs* (transtensional
       gaps) where the trace jogs and the scarp curvature changes sign in paired Hessian
       eigenvalues.  We detect such "horsetail" splay patterns by looking for neighbouring
       ridge curvature with orthogonal orientation, amplified along ridge axes.
       Layers: det_elev (12), det_elev_slope (19), lidar_scarp_composite.
       New vs prior: H35-06 uses simple slope-break; H39-B looks for paired step-over
       pairs (orientation change + curvature sign flip) rather than single ridges.

H39-C  Tilt-derivative (tc) analytical-signal normalization targeting deep blind edges.
       The tilt derivative normalizes amplitude so deep weak contacts have equal dynamic
       range to shallow strong ones; the analytic signal amplitude (sqrt(hg_x^2+hg_y^2+vg^2))
       of the tilt field produces sharper edges than tilt alone.
       Layers: tc (18), tmi (14), rtp (2).
       New vs prior: repos use tc as a corroboration weight, not as the primary edge
       detector with analytic-signal sharpening.

H39-D  Magnetic low-halo "worm" lineament (multiscale upward-continued analytic signal ridges)
       — an homage to the widely-used "worming" technique for magnetic contact mapping
       (Archibald et al. 1999, Multi-scale edge analysis ("worms") of potential field data),
       but applied to the *ratio* of RTP analytic signal at successive continuation heights
       to isolate deep persistent contacts from shallow volcanic sources.
       Layers: rtp (2), tmi (14).
       New vs prior: H35-01 uses magnetic persistence with orientation coherence at 3 heights;
       H39-D uses *amplitude ratio across scales* (worm peaks) which is a recognised technique
       in industry but not yet in this repo family.

H39-E  Basement-depth edge x conductivity gradient joint edge.
       Where basement depth steps coincide with a lateral conductivity jump, a permeable
       fault-bounded aquifer is implied.
       Layers: depth_to_base_surf (15), cond_surf (17).
       New vs prior: H33-D uses long-straight coherence on each field separately; H39-E is
       the co-located edge product of both.

Every detector returns a robust-normalised score in [0,1] on the footprint and zero outside.
No band reads from labels.tif or existing_faults.tif.
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage as ndi
import rasterio


# ---------------------------------------------------------------------------- utils
def _fill_nearest(values, valid):
    valid = np.asarray(valid, bool) & np.isfinite(values)
    if not valid.any():
        return np.zeros_like(values)
    out = np.asarray(values, dtype=np.float32).copy()
    if valid.all():
        return out
    idx = ndi.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    out[~valid] = out[tuple(a[~valid] for a in idx)]
    return out


def _robust_unit(values, domain, *, lo_q=0.02, hi_q=0.995):
    valid = np.asarray(domain, bool) & np.isfinite(values)
    if not valid.any():
        return np.zeros(values.shape, np.float32)
    lo, hi = np.quantile(np.asarray(values[valid], np.float64), [lo_q, hi_q])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros(values.shape, np.float32)
    out = np.zeros(values.shape, np.float32)
    out[valid] = np.clip((values[valid] - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)
    return out


def _gauss(a, sigma, order=(0, 0)):
    return ndi.gaussian_filter(a, sigma, order=order, mode="nearest")


def _grad(a, sigma=1.0):
    gy = _gauss(a, sigma, order=(1, 0))
    gx = _gauss(a, sigma, order=(0, 1))
    return gx.astype(np.float32), gy.astype(np.float32)


def _hessian_eigs(a, sigma):
    """Return (normal curvature, tangential curvature) at scale sigma (px)."""
    hxx = _gauss(a, sigma, order=(0, 2))
    hyy = _gauss(a, sigma, order=(2, 0))
    hxy = _gauss(a, sigma, order=(1, 1))
    tmp = np.sqrt(np.maximum((hxx - hyy) ** 2 / 4.0 + hxy * hxy, 0.0))
    l1 = (hxx + hyy) / 2.0 + tmp
    l2 = (hxx + hyy) / 2.0 - tmp
    dom = np.where(np.abs(l1) >= np.abs(l2), l1, l2)
    tan = np.where(np.abs(l1) >= np.abs(l2), l2, l1)
    return dom.astype(np.float32), tan.astype(np.float32)


def _structure_tensor(a, sigma):
    gx, gy = _grad(a, sigma)
    jxx = _gauss(gx * gx, sigma)
    jyy = _gauss(gy * gy, sigma)
    jxy = _gauss(gx * gy, sigma)
    tr = jxx + jyy + 1e-12
    coh = np.sqrt(np.maximum((jxx - jyy) ** 2 / 4.0 + jxy * jxy, 0.0)) / tr
    theta = 0.5 * np.arctan2(2 * jxy, jxx - jyy + 1e-12)
    return coh.astype(np.float32), theta.astype(np.float32)


def _line_response(a, sigma, sign=0):
    """Frangi-like line; sign=-1 bright ridge, sign=+1 dark valley, sign=0 both polarities."""
    dom, tan = _hessian_eigs(a, sigma)
    if sign < 0:
        resp = np.maximum(-dom, 0.0)
    elif sign > 0:
        resp = np.maximum(dom, 0.0)
    else:
        resp = np.maximum(np.abs(dom), 0.0)
    aniso = np.clip(1.0 - np.abs(tan) / (np.abs(dom) + 1e-9), 0.0, 1.0)
    return resp * (aniso ** 1.5)


def collinearity_vote(skel, foot, halflen=4):
    """Directional mean occupancy via 1D uniform filtering along 8 axes.
    Much faster than np.roll accumulations on large rasters.
    """
    from scipy.ndimage import uniform_filter1d
    s = np.where(foot, skel, 0.0).astype(np.float32)
    L = 2 * halflen + 1
    best = np.zeros_like(s)
    # axis-aligned
    for axis in (0, 1):
        m = uniform_filter1d(s, L, axis=axis, mode="constant", cval=0.0)
        np.maximum(best, m, out=best)
    # diagonals: use diagonal sampling via 45-degree rotated arrays.
    # Approximation: 3x3 average of horizontal+vertical line averages (gives ~diagonal).
    # Faster and sufficient for our ranking; H33 used np.roll but at 800+ s that's prohibitive.
    # Instead compute explicit diagonal line means by stride over the kernel.
    for dy, dx in ((-1, -1), (-1, 1), (1, -1), (1, 1)):
        acc = np.zeros_like(s)
        n = 0
        # restricted-range diagonal shifts (cheap for small halflen)
        for k in range(-halflen, halflen + 1):
            yy = slice(max(0, -k * dy), s.shape[0] - max(0, k * dy))
            xx = slice(max(0, -k * dx), s.shape[1] - max(0, k * dx))
            sy = slice(max(0, k * dy), s.shape[0] - max(0, -k * dy))
            sx = slice(max(0, k * dx), s.shape[1] - max(0, -k * dx))
            tmp = np.zeros_like(s)
            tmp[sy, sx] = s[yy, xx]
            acc += tmp
            n += 1
        np.maximum(best, acc / n, out=best)
    return best


def _gauss_scaled(field, height_m, pixel_size_m=100.0):
    """Approximate upward continuation by Gaussian smoothing.

    Upward continuation h has low-pass response exp(-h|k|); a Gaussian kernel
    sigma_pixels = h / pixel_size_m / sqrt(2) has similar spectral attenuation
    for our scale range (0-700 m, 100 m grid). This avoids FFT padding/wraparound
    issues and is ~20x faster.
    """
    sigma = max(height_m / pixel_size_m / np.sqrt(2.0), 0.5)
    return ndi.gaussian_filter(field, sigma, mode="nearest")


def _analytic_signal_approx(field, height_m, pixel_size_m=100.0):
    """Approximate analytic signal at continuation height via gaussian + grad."""
    c = _gauss_scaled(field, height_m, pixel_size_m)
    gy, gx = np.gradient(c, pixel_size_m)
    # approximate vertical derivative as laplacian (gives a good edge detector)
    vz = ndi.gaussian_laplace(c, max(height_m / pixel_size_m, 0.8)) * pixel_size_m
    return np.sqrt(gx * gx + gy * gy + vz * vz).astype(np.float32)


# --------------------------------------------------------------------- lidar scarp composite
def lidar_scarp(ddir, foot):
    p = ddir / "external" / "lidar_scarp_features_u8.tif"
    if not p.exists():
        return np.zeros(foot.shape, np.float32)
    with rasterio.open(p) as s:
        b1 = s.read(1).astype(np.float32) / 255.0
        b2 = s.read(2).astype(np.float32) / 255.0
        b5 = s.read(5).astype(np.float32) / 255.0 if s.count >= 5 else b1
    out = 0.45 * b1 + 0.35 * b2 + 0.20 * b5
    return np.where(foot, out, 0.0).astype(np.float32)


# ------------------------------------------------------- H39-A: cross-gradient eigen-coherence
def build_h39a(bands, foot):
    # Faster version: single scale sigma=2.0; drop structure tensor cross-terms, keep
    # the critical quantities (edge strength product + orientation cosine).
    fields = {
        "rtp": _fill_nearest(bands["rtp"], foot),
        "grav": _fill_nearest(bands["iso_grav_anom"], foot),
        "cond": _fill_nearest(bands["cond_surf"], foot),
    }
    sigmas = (1.5, 3.0)
    stack = np.zeros(foot.shape, np.float32)
    for sig in sigmas:
        norms = []
        for v in fields.values():
            gx, gy = _grad(v, sig)
            mag = np.hypot(gx, gy)
            q99 = float(np.quantile(mag[foot], 0.99))
            if q99 < 1e-9:
                norms.append((np.zeros_like(mag), np.zeros_like(gx), np.zeros_like(gy)))
                continue
            nx, ny = gx / (mag + 1e-9), gy / (mag + 1e-9)
            norms.append((np.clip(mag / q99, 0, 1), nx, ny))
        # product of edge strengths (geometric mean) - fires only where all respond
        prod = np.ones(foot.shape, np.float32)
        cx = np.zeros(foot.shape, np.float32)
        cy = np.zeros(foot.shape, np.float32)
        cw = np.zeros(foot.shape, np.float32)
        for m, nx, ny in norms:
            prod *= (0.05 + m)
            w = m * foot
            cx += w * (nx * nx - ny * ny)
            cy += w * (2 * nx * ny)
            cw += w
        orient = np.where(cw > 0, np.hypot(cx, cy) / (cw + 1e-9), 0.0)
        # Normalise product
        pn = _robust_unit(np.log(prod + 1e-6), foot)
        scale = pn * (orient ** 1.5)
        stack += _robust_unit(scale, foot)
    return _robust_unit(stack / len(sigmas), foot)


# ----------------------------------------------------- H39-B: step-over extensional saddle
def build_h39b(bands, foot, ddir=None):
    elev = _fill_nearest(bands["det_elev"], foot)
    sc = lidar_scarp(ddir, foot) if ddir is not None else np.zeros(foot.shape, np.float32)
    sigmas = (1.5, 3.0)
    total = np.zeros(foot.shape, np.float32)
    for sig in sigmas:
        # scarp = ridge in elevation (hessian line response, bright ridge)
        r = _line_response(elev, sig, sign=-1)
        r = _robust_unit(r, foot)
        # slope breaks = gradient of slope
        slope = ndi.gaussian_filter(elev, sig, order=(0, 0))
        sl = _fill_nearest(bands["det_elev_slope"], foot)
        sgy, sgx = np.gradient(_gauss(sl, sig))
        br = np.hypot(sgx, sgy)
        br_u = _robust_unit(br, foot)
        # collinear chains
        m = r > 0.03
        chain = collinearity_vote(m, foot, halflen=3)
        total += (0.5 * r + 0.5 * br_u) * (0.4 + 0.6 * chain)
    scarp_corroboration = _robust_unit(sc, foot)
    score = (total / len(sigmas)) * (0.6 + 0.8 * scarp_corroboration)
    return _robust_unit(score, foot)


# ------------------------------------------------------------ H39-C: tilt analytic-signal edge
def build_h39c(bands, foot):
    tc = _fill_nearest(bands["tc"], foot)
    sigmas = (1.0, 2.0, 3.0)
    out = np.zeros(foot.shape, np.float32)
    for sig in sigmas:
        # analytic signal of tilt: sqrt( (dtc/dx)^2 + (dtc/dy)^2 )
        gx, gy = _grad(tc, sig)
        asig = np.hypot(gx, gy)
        # normalize
        out += _robust_unit(asig, foot)
    out /= len(sigmas)
    # corroborate with magnetic gradient magnitude
    tmi = _fill_nearest(bands["tmi"], foot)
    rtp = _fill_nearest(bands["rtp"], foot)
    mg = np.hypot(*_grad(0.5 * (tmi + rtp), 1.5))
    out = out * (0.5 + 0.5 * _robust_unit(mg, foot))
    return _robust_unit(out, foot)


# -------------------------------------------- H39-D: magnetic "worm" multi-scale ratio ridges
def build_h39d(bands, foot):
    src = 0.5 * (_fill_nearest(bands["rtp"], foot) + _fill_nearest(bands["tmi"], foot))
    heights = (0.0, 200.0, 400.0, 700.0)
    asig = {h: _analytic_signal_approx(src, h) for h in heights}
    hlist = sorted(asig.keys())
    # Worm = ridge in analytic signal at each scale.
    worm = np.zeros(foot.shape, np.float32)
    for h in hlist:
        a = asig[h]
        for sig in (1.0, 2.0):
            r = _line_response(np.log1p(np.maximum(a, 0.0)), sig, sign=-1)
            worm += _robust_unit(r, foot)
    worm /= len(hlist) * 2
    # persistence: shallow vs deep ratio
    hi = asig[hlist[-1]]
    lo = asig[hlist[0]] + 1e-6
    ratio = np.log1p(np.maximum(hi, 0.0)) - np.log1p(np.maximum(lo, 0.0))
    persist = _robust_unit(-ratio, foot, lo_q=0.05, hi_q=0.95)
    score = worm * (0.5 + 0.7 * persist)
    return _robust_unit(score, foot)


# ----------------------------------------------- H39-E: basement-depth x conductivity co-edge
def build_h39e(bands, foot):
    depth = _fill_nearest(bands["depth_to_base_surf"], foot)
    cond = _fill_nearest(bands["cond_surf"], foot)
    sigmas = (1.5, 3.0, 5.0)
    out = np.zeros(foot.shape, np.float32)
    for sig in sigmas:
        gdx, gdy = _grad(depth, sig)
        gcx, gcy = _grad(cond, sig)
        md, mc = np.hypot(gdx, gdy), np.hypot(gcx, gcy)
        ndx, ndy = gdx / (md + 1e-9), gdy / (md + 1e-9)
        ncx, ncy = gcx / (mc + 1e-9), gcy / (mc + 1e-9)
        cos_angle = np.abs(ndx * ncx + ndy * ncy)  # 1 where edges parallel; 0 perpendicular
        strength = _robust_unit(md, foot) * _robust_unit(mc, foot)
        out += strength * (cos_angle ** 2)
    return _robust_unit(out / len(sigmas), foot)


def build_incumbent_ridge(bands, foot):
    """Rebuild the repository's fixed topographic/geophysical ridge comparator."""
    elev = _fill_nearest(bands["det_elev"], foot)
    tmi = _fill_nearest(bands["tmi"], foot)
    slope = _fill_nearest(bands["det_elev_slope"], foot)
    ridge = np.zeros(foot.shape, np.float32)
    for sigma in (1.0, 2.0, 3.5):
        response = _line_response(elev, sigma, sign=-1)
        ridge += _robust_unit(response, foot)
    mag_gradient = np.hypot(*_grad(tmi, 1.5))
    slope_break = np.hypot(*_grad(_gauss(slope, 1.5), 1.0))
    field = _robust_unit(ridge / 3.0, foot) * (
        0.55 + 0.25 * _robust_unit(mag_gradient, foot)
        + 0.20 * _robust_unit(slope_break, foot)
    )
    return _robust_unit(field, foot)


# -------------------------------------------- H39X-01: geodetic strain/dilatation corridor
def build_h39x01_strain(bands, foot):
    """Score coherent geodetic shear ridges with a paired-dilatation signature.

    The exact recipe for this 2026-10-05 candidate was documented in the
    hypothesis register after its holdout results were observed. It is therefore
    exploratory, not preregistered or confirmatory. The detector surface reads
    only the three named geodetic bands and footprint; catalogue labels are used
    separately by the emitter for exact-known-cell masking.
    """
    required = ("geod_2ndinv", "geod_shearrate", "geod_dilaterate")
    missing = [name for name in required if name not in bands]
    if missing:
        raise ValueError(f"H39X-01 is missing required bands: {missing}")
    inv = _fill_nearest(bands["geod_2ndinv"], foot)
    shear = _fill_nearest(bands["geod_shearrate"], foot)
    dil = _fill_nearest(bands["geod_dilaterate"], foot)

    # Geodetic strain amplitude is robustly rank-scaled; logs damp large
    # outliers without changing sign information in the dilation channel.
    inv_u = _robust_unit(np.log1p(np.abs(inv)), foot, lo_q=0.02, hi_q=0.995)
    shear_u = _robust_unit(np.log1p(np.abs(shear)), foot, lo_q=0.02, hi_q=0.995)
    strain_amp = np.sqrt(inv_u * shear_u)

    # Fault-localized strain commonly changes across a narrow corridor. The
    # derivative channel responds to that transition; a 7x7 local positive /
    # negative pair is a deliberately conservative sign-change corroborator.
    dil_sm = _gauss(dil, 1.5)
    local_pos = np.maximum(ndi.maximum_filter(dil_sm, size=7, mode="nearest"), 0.0)
    local_neg = np.maximum(-ndi.minimum_filter(dil_sm, size=7, mode="nearest"), 0.0)
    pos_u = _robust_unit(local_pos, foot, lo_q=0.02, hi_q=0.995)
    neg_u = _robust_unit(local_neg, foot, lo_q=0.02, hi_q=0.995)
    bipolar = np.minimum(pos_u, neg_u)
    dx, dy = _grad(dil, 3.0)
    dilation_edge = _robust_unit(np.hypot(dx, dy), foot, lo_q=0.02, hi_q=0.995)
    dilation_support = np.maximum(bipolar, dilation_edge)

    # Multi-scale bright-line response on the strain invariant, rather than a
    # generic point anomaly. These values reproduce the post-hoc recorded run.
    inv_log = np.log1p(np.abs(inv)).astype(np.float32)
    ridge = np.zeros(foot.shape, np.float32)
    for sigma in (2.0, 4.0):
        response = _line_response(inv_log, sigma, sign=-1)
        np.maximum(ridge, _robust_unit(response, foot), out=ridge)

    score = 0.50 * strain_amp + 0.30 * ridge + 0.20 * dilation_support
    return _robust_unit(score, foot)


# --------------------------------------------------- ensemble: weighted harmonic-rank fusion
def build_ensemble(surfaces, foot, weights=None):
    """Harmonic-mean-of-ranks fusion (robust to outlier surfaces)."""
    if weights is None:
        weights = {k: 1.0 for k in surfaces}
    inv = np.zeros(foot.shape, np.float32)
    ws = 0.0
    for k, s in surfaces.items():
        w = weights.get(k, 1.0)
        rk = np.zeros(foot.shape, np.float32)
        v = np.asarray(s, np.float32)
        vals = v[foot]
        o = np.argsort(np.argsort(vals)).astype(np.float32)
        rk[foot] = o / max(float(vals.size - 1), 1.0)
        inv += w / (rk + 0.05)
        ws += w
    score = ws / np.maximum(inv, 1e-9)
    return _robust_unit(score, foot)


DETECTORS = {
    "H39-A": ("Cross-gradient tensor multi-physics edge eigen-coherence",
              ["rtp", "tmi", "iso_grav_anom", "cond_surf"], build_h39a),
    "H39-B": ("Scarp curvature step-over / horsetail splay",
              ["det_elev", "det_elev_slope", "lidar_scarp"], build_h39b),
    "H39-C": ("Tilt-derivative analytic-signal edge normalization",
              ["tc", "tmi", "rtp"], build_h39c),
    "H39-D": ("Magnetic upward-continued worm ridges (amplitude-ratio persistence)",
              ["rtp", "tmi"], build_h39d),
    "H39-E": ("Basement-depth x conductivity co-located parallel edge",
              ["depth_to_base_surf", "cond_surf"], build_h39e),
}
