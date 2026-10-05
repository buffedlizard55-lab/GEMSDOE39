"""Metric-priced emission: the budget is a *price*, not a copied constant.

Prior rounds in this project's family emitted a fixed budget (44 090, then
40 199, then 37 654) because that was what the previous artifact happened to
contain.  The calibrated inverse-DTI model in ``calibrate.py`` makes the budget
a derived quantity instead:

  * a dot is worth adding iff its marginal hit probability pi exceeds
    ``calibrate.break_even_pi`` (about 3.2% at DTI 0.2778);
  * the hit probability decays with rank down the propensity list, and that
    decay is *measured* on the prevalence-matched holdout instrument rather
    than assumed;
  * emission therefore stops at the rank where the measured marginal hit rate
    crosses the break-even price, transferred to the live scale by the anchor
    artifact's instrument-to-live ratio.

Geometry is unchanged from the empirically validated winning family:
best-first Poisson-disk thinning at ~2.8 px (280 m) minimum spacing, with a
hard exclusion of any pixel within ``cat_buffer_px`` of the supplied catalogue.
The exclusion is justified by an organizer-scored differential, not by theory:
deleting the 2 545 dots at 1.41-2.00 px from the catalogue moved 0.2708 ->
0.2778 (see ``calibrate.verify_nested_pair``).  Staff confirmed the mask is
pixel-exact, so those dots were not masked -- they were ordinary false
positives, and the near-catalogue ring is where new truth is empirically
sparsest.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from .calibrate import break_even_pi, forward, price


def allowed_domain(foot, catalogue=None, cat_buffer_px=2):
    ok = foot.copy()
    if catalogue is not None and cat_buffer_px > 0:
        d = distance_transform_edt(~(catalogue & foot))
        ok &= d > cat_buffer_px
    return ok


def emit_ranked(prop, allowed, min_dist=2.8, max_n=200_000, pool_mult=14):
    """Best-first Poisson-disk sampling.

    Returns (mask, rank_of_chosen) where ``rank_of_chosen[i]`` is the position
    (0-based) of the i-th accepted dot in the descending-propensity order.  The
    rank array is what makes the stopping rule auditable: the marginal hit rate
    can be plotted against rank instead of being asserted.
    """
    prop = np.asarray(prop, np.float32)
    H, W = prop.shape
    flat = np.where(allowed, prop, np.float32(-np.inf)).astype(np.float32).reshape(-1)
    pool = int(min(pool_mult * max_n, int(allowed.sum())))
    pool = max(pool, 1)
    idx = np.argpartition(flat, -pool)[-pool:]
    vals = flat[idx]
    order = np.argsort(-vals, kind="stable")
    idx_sorted = idx[order]
    r2 = min_dist * min_dist
    cell = max(1.0, float(min_dist))
    grid = {}
    mask = np.zeros((H, W), bool)
    ranks = []
    for i in range(idx_sorted.size):
        fi = int(idx_sorted[i])
        if not np.isfinite(flat[fi]) or flat[fi] <= 0:
            continue
        y, x = fi // W, fi % W
        if not allowed[y, x]:
            continue
        cy, cx = int(y // cell), int(x // cell)
        ok = True
        for gy in range(cy - 2, cy + 3):
            for gx in range(cx - 2, cx + 3):
                for (ky, kx) in grid.get((gy, gx), ()):
                    if (y - ky) ** 2 + (x - kx) ** 2 < r2:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            mask[y, x] = True
            ranks.append(i)
            grid.setdefault((cy, cx), []).append((y, x))
            if len(ranks) >= max_n:
                break
    return mask, np.asarray(ranks, np.int64)


def dot_hits(mask, truth, r_px=3.0, active=None):
    """Per-dot hit indicator: a dot 'hits' if a truth pixel lies within r_px."""
    ys, xs = np.nonzero(mask)
    if ys.size == 0 or not np.any(truth):
        return ys, xs, np.zeros(ys.size, bool)
    d = distance_transform_edt(~truth)
    hit = d[ys, xs] <= r_px
    if active is not None:
        hit &= active[ys, xs]
    return ys, xs, hit


def marginal_hit_rate(hit_sorted, window=1500):
    """Centred rolling hit rate over the emission order (marginal, not average)."""
    n = hit_sorted.size
    if n == 0:
        return np.zeros(0)
    c = np.concatenate([[0.0], np.cumsum(hit_sorted.astype(np.float64))])
    half = window // 2
    lo = np.maximum(np.arange(n) - half, 0)
    hi = np.minimum(np.arange(n) + half + 1, n)
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)


def choose_budget(hit_sorted, pi_star, window=1500, min_n=5000, transfer=1.0,
                  max_n=None):
    """Largest n whose *marginal* (windowed) hit rate is still above the price.

    ``transfer`` converts an instrument hit rate into the predicted live hit
    rate (live = instrument / transfer).  ``transfer`` is estimated from the
    live-scored anchor artifact, so it is measured, not assumed -- but it is
    estimated from a single anchor and that limitation is reported.
    """
    if hit_sorted.size == 0:
        return min_n, dict(reason="no dots")
    marg = marginal_hit_rate(hit_sorted, window) / max(transfer, 1e-9)
    n = int(hit_sorted.size if max_n is None else min(max_n, hit_sorted.size))
    cut = n
    for i in range(min_n, n):
        if marg[i] < pi_star:
            cut = i
            break
    else:
        cut = n
    return int(cut), dict(pi_star=pi_star, window=window, transfer=transfer,
                          marginal_at_cut=float(marg[max(cut - 1, 0)]),
                          mean_hit_rate_over_cut=float(hit_sorted[:cut].mean()) if cut else 0.0,
                          max_available=n)


def trim_to_budget(mask, ranks, n):
    """Keep only the first ``n`` accepted dots in emission order."""
    if n >= ranks.size:
        return mask.copy()
    keep = np.zeros(mask.shape, bool)
    ys, xs = np.nonzero(mask)
    # recompute acceptance order deterministically: sort chosen dots by rank
    order = np.argsort(ranks)
    ys, xs = ys[order][:n], xs[order][:n]
    keep[ys, xs] = True
    return keep


def price_candidate(n, hit_rate, G, l0=3.0, kb=0.5):
    return price(n, G, hit_rate, l0=l0, kb=kb)
