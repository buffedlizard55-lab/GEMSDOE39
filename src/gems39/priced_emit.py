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

    Returns ``(mask, ys, xs)`` where ``ys[i], xs[i]`` is the i-th ACCEPTED dot, in
    descending-propensity (best-first) order.  Returning the coordinates in
    acceptance order is what makes the stopping rule and any budget trim correct.

    DEFECT FIXED 2026-10-05: this used to return ``(mask, rank_of_chosen)`` where
    ``rank_of_chosen[i]`` was the position of the i-th accepted dot in the
    *candidate pool*.  Callers then did ``ys, xs = np.nonzero(mask)`` -- which is
    RASTER order -- and reordered with ``ys[np.argsort(ranks)]``, indexing a
    raster-ordered coordinate array with pool-order positions.  The result was a
    scrambled subset: ``trim_to_budget`` did not return the best n dots, it
    returned n dots chosen by an arbitrary permutation, and the "marginal hit rate
    curve" the stopping rule reads was noise around the mean.  Both callers and
    the emitted artifact were affected.
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
    ys_out = np.empty(max_n, np.int32)
    xs_out = np.empty(max_n, np.int32)
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
            k = len(ranks)
            ys_out[k] = y
            xs_out[k] = x
            ranks.append(i)
            grid.setdefault((cy, cx), []).append((y, x))
            if k + 1 >= max_n:
                break
    k = len(ranks)
    return mask, ys_out[:k], xs_out[:k]


def dot_hits(mask, truth, r_px=3.0, active=None, ys=None, xs=None):
    """Per-dot hit indicator: a dot 'hits' if a truth pixel lies within r_px.

    Pass ``ys``/``xs`` (as returned by ``emit_ranked``) to keep the result in
    ACCEPTANCE order; omit them to get raster order from the mask.
    """
    if ys is None or xs is None:
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


def trim_to_budget(shape, ys, xs, n):
    """Keep only the first ``n`` accepted dots in best-first emission order.

    ``ys``/``xs`` must be the acceptance-ordered coordinates returned by
    ``emit_ranked``.  See the defect note there: an earlier signature took the
    mask plus a pool-rank array and silently returned a permuted subset.
    """
    n = int(min(max(n, 0), ys.size))
    keep = np.zeros(shape, bool)
    if n:
        keep[ys[:n], xs[:n]] = True
    return keep


def price_candidate(n, hit_rate, G, l0=3.0, kb=0.5):
    return price(n, G, hit_rate, l0=l0, kb=kb)


def choose_budget_priced(hit_sorted, G, transfer=1.0, s_start=0.2778, window=1500,
                         min_n=1000, l0=3.0, kb=0.5, iters=12, tol=50):
    """Self-consistent emission budget: iterate the break-even price to a fixed point.

    DEFECT FIXED 2026-10-05: ``choose_budget`` compared the marginal hit rate
    against ``pi*`` evaluated at the ANCHOR's score (DTI 0.2778 -> pi* = 3.149%).
    But ``pi*`` rises steeply with the operating score -- at DTI 0.4774 it is
    10.5%, not 3.15% -- because the saturating coverage term ``l_marg`` shrinks as
    more of |G| is already credited.  Using the anchor's price therefore told a
    strong candidate to keep dotting long after its own marginal dot had become
    score-negative, and the "priced" arm ran to the pool maximum at a LOWER price
    than the anchor-budget arm.  That contradiction is what exposed the bug.

    The budget and the price determine each other, so they are solved together:

        s_0     = anchor live score
        repeat: pi* = break_even_pi(s_k, G, hits(n_k))
                n_{k+1} = first rank whose windowed marginal live hit rate < pi*
                s_{k+1} = price(n_{k+1}, G, mean_hit_rate(n_{k+1}) / transfer)
        until |n_{k+1} - n_k| <= tol

    ``hit_sorted`` is the per-dot hit indicator in best-first acceptance order on
    the local instrument; ``transfer`` converts an instrument hit rate to the live
    scale.  Returns (n, receipt) where the receipt carries the whole iteration
    path so the stopping point can be audited.
    """
    from .calibrate import break_even_pi as _bep, price as _price
    hs = np.asarray(hit_sorted, bool)
    n_avail = hs.size
    if n_avail == 0:
        return int(min_n), dict(reason="no dots emitted", iterations=[])
    cum = np.concatenate([[0], np.cumsum(hs.astype(np.float64))])
    marg = marginal_hit_rate(hs, window) / max(transfer, 1e-9)
    s = float(s_start)
    n_cut = None            # must start unset: initialising it to n_avail made the
                            # first iteration "converge" against itself and return
                            # the whole pool without ever re-pricing at the new score
    path = []
    new_n = None
    for it in range(iters):
        # price the CURRENT iterate so pi* is evaluated at the score this budget
        # actually achieves, not at the anchor's
        if n_cut is not None:
            h_now = float(cum[n_cut] / max(n_cut, 1))
            s = _price(n_cut, G, h_now / max(transfer, 1e-9), l0=l0, kb=kb)["score"]
        hits_live = (cum[n_cut] / max(transfer, 1e-9)) if n_cut else 0.0
        bep = _bep(s, G, hits_live, l0=l0, kb=kb)
        pi = bep["pi_star"]
        below = np.nonzero(marg[min_n:n_avail] < pi)[0]
        new_n = int(min_n + below[0]) if below.size else int(n_avail)
        hits_at_n = cum[new_n] / max(transfer, 1e-9)
        h_instr = float(cum[new_n] / new_n)
        pr = _price(new_n, G, h_instr / max(transfer, 1e-9), l0=l0, kb=kb)
        path.append(dict(iteration=it, score_in=s, pi_star=pi,
                         l_marginal_TP=bep["l_marginal_TP"], n_out=new_n,
                         mean_instrument_hit_rate=h_instr,
                         mean_live_hit_rate=h_instr / max(transfer, 1e-9),
                         score_out=pr["score"]))
        converged = n_cut is not None and abs(new_n - n_cut) <= tol
        n_cut, s = new_n, pr["score"]
        if converged:
            break
    return n_cut, dict(n=n_cut, score=s, transfer=transfer, G=G, window=window,
                       min_n=min_n, l0=l0, kb=kb, iterations=len(path), path=path,
                       converged=bool(new_n is not None and abs(new_n - n_cut) <= tol),
                       marginal_at_cut=float(marg[max(n_cut - 1, 0)]),
                       pi_star_at_cut=float(path[-1]["pi_star"]) if path else None,
                       n_available=int(n_avail))


def optimal_budget_from_curve(hit_sorted, G, transfer=1.0, l0=3.0, kb=0.5,
                              n_grid=400, min_n=500):
    """The budget as the argmax of the calibrated price curve.

    ``choose_budget_priced`` iterates "budget -> price -> break-even price ->
    budget" to a fixed point.  That map has no stable fixed point when coverage
    saturates: as more of |G| is credited, ``l_marginal`` shrinks and
    ``pi* = 0.2 s / [l_marg (1-0.2 s) + 0.2 s kb]`` rises past 1.0, which correctly
    means "no further dot is worth adding" but makes the iteration oscillate
    between the whole pool and the floor.  Rather than damp it, this function does
    the thing the iteration was trying to approximate: evaluate the price at every
    budget and take the maximum.

        h_instr(n) = (hits among the first n accepted dots) / n
        s(n)       = price(n, G, h_instr(n) / transfer)
        n*         = argmax_n s(n)

    The break-even price at ``n*`` is reported next to the measured marginal hit
    rate at ``n*`` so the two can be compared: at an interior optimum they should
    be close, and if they are not, that disagreement is printed rather than hidden.

    Returns ``(n_star, receipt)``.
    """
    from .calibrate import break_even_pi as _bep, price as _price
    hs = np.asarray(hit_sorted, bool)
    n_avail = hs.size
    if n_avail <= min_n:
        return int(n_avail), dict(reason="pool smaller than min_n", n_available=n_avail)
    cum = np.concatenate([[0], np.cumsum(hs.astype(np.float64))])
    ns = np.unique(np.linspace(min_n, n_avail, n_grid).astype(int))
    best = None
    curve = []
    for n in ns:
        h_instr = float(cum[n] / n)
        pr = _price(int(n), G, h_instr / max(transfer, 1e-9), l0=l0, kb=kb)
        curve.append((int(n), h_instr, pr["score"]))
        if best is None or pr["score"] > best[2]["score"]:
            best = (int(n), h_instr, pr)
    n_star, h_star, pr_star = best
    marg = marginal_hit_rate(hs, 1500)
    i = max(n_star - 1, 0)
    marg_live = float(marg[i] / max(transfer, 1e-9))
    bep = _bep(pr_star["score"], G, cum[n_star] / max(transfer, 1e-9), l0=l0, kb=kb)
    return n_star, dict(n=n_star, score=pr_star["score"], price=pr_star,
                        mean_instrument_hit_rate=h_star,
                        mean_live_hit_rate=h_star / max(transfer, 1e-9),
                        marginal_live_hit_rate_at_optimum=marg_live,
                        pi_star_at_optimum=bep["pi_star"],
                        l_marginal_TP=bep["l_marginal_TP"],
                        marginal_agrees_with_price=bool(
                            abs(marg_live - bep["pi_star"]) < 0.5 * max(bep["pi_star"], 1e-9)),
                        transfer=transfer, G=G, n_available=int(n_avail),
                        curve=[dict(n=a, mean_instrument_hit_rate=b, score=c)
                               for a, b, c in curve])
