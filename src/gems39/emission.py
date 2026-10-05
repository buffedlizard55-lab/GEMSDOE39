"""Emission: convert a [0,1] detector field into a binary prediction set.

Two emitters are provided:

  * ``fast_poisson``: top-q support then raster-order Poisson-disk thinning.
    This is O(N) and gives the canonical "dotted ridge" layout of prior winning
    submissions (d2.8, 0.2600-0.2708 scores) with a minimum spacing ~2.8 px.
  * ``greedy_maxcover``: metric-native lazy greedy (more principled but ~N*29
    shifts per iteration; fine for budgets up to ~50k and is used here only
    for the final primary after ranking with the fast emitter).

Catalogue-buffer exclusion is a hard gate (no pixels within ``cat_buffer_px``
(>=2 px) of a known fault, matching the H33-2-B2 0.2778 rule, because known
faults are masked at scoring time so any pixel placed on them is wasted FP).
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import distance_transform_edt
from .metric import R_PX, kernel


def _offsets():
    r = int(np.ceil(R_PX))
    dy, dx, k = [], [], []
    for yy in range(-r, r + 1):
        for xx in range(-r, r + 1):
            kk = float(kernel(np.hypot(yy, xx)))
            if kk > 0:
                dy.append(yy); dx.append(xx); k.append(kk)
    return np.array(dy, int), np.array(dx, int), np.array(k, np.float32)


_DY, _DX, _KW = _offsets()


def _shift(a, dy, dx):
    out = np.zeros_like(a)
    H, W = a.shape
    y0, y1 = max(0, dy), min(H, H + dy)
    x0, x1 = max(0, dx), min(W, W + dx)
    out[y0:y1, x0:x1] = a[max(0, -dy):min(H, H - dy), max(0, -dx):min(W, W - dx)]
    return out


def dot_thin_poisson(support, min_dist=2.8, order=None):
    """Poisson-disk thinning over a binary support in field-rank order.

    If ``order`` is provided, pixels are visited in descending order of
    ``order`` (higher values first) — the classic "best-first" disk sample.
    Otherwise raster order is used.
    """
    support = np.asarray(support, bool)
    if support.sum() == 0:
        return np.zeros(support.shape, bool)
    if order is None:
        ys, xs = np.nonzero(support)
    else:
        ys, xs = np.nonzero(support)
        o = np.asarray(order, np.float32)[ys, xs]
        perm = np.argsort(-o)
        ys, xs = ys[perm], xs[perm]
    r2 = min_dist * min_dist
    cell = max(1.0, float(min_dist))
    grid = {}
    out = np.zeros(support.shape, bool)
    r = int(np.ceil(min_dist))
    for i in range(ys.size):
        y, x = int(ys[i]), int(xs[i])
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
            out[y, x] = True
            grid.setdefault((cy, cx), []).append((y, x))
    return out


def _catalogue_block(catalogue, foot, cat_buffer_px):
    if catalogue is None or cat_buffer_px <= 0:
        return np.zeros(foot.shape, bool)
    d = distance_transform_edt(~(catalogue & foot))
    return (d <= cat_buffer_px)


def emit_fast(field, foot, target_n, min_dist=2.7, catalogue=None, cat_buffer_px=2):
    """Match-budget fast best-first Poisson-disk emission.

    Take the top `candidate_mult * target_n` pixels by field value and iterate
    in descending order, adding each pixel if it is farther than `min_dist`
    from any already-added dot. This is best-first thinning over a pre-screened
    candidate set (O(cand_n) with spatial hash grid), much faster than
    quantile-binary-search raster-order variants.
    """
    candidate_mult = 6
    f = np.asarray(field, np.float32)
    f = np.where(foot, f, -np.inf).astype(np.float32)
    blocked = _catalogue_block(catalogue, foot, cat_buffer_px)
    allowed = foot & ~blocked
    flat = f.reshape(-1)
    cand_n = int(min(candidate_mult * target_n, int(allowed.sum())))
    # partial sort for top-k
    idx = np.argpartition(flat, -cand_n)[-cand_n:]
    vals = flat[idx]
    # sort descending
    order = np.argsort(-vals)
    idx_sorted = idx[order]
    r2 = min_dist * min_dist
    cell = max(1.0, float(min_dist))
    grid_d = {}
    out = np.zeros(f.shape, bool)
    H, W = f.shape
    chosen = 0
    r = int(np.ceil(min_dist))
    for flat_i in idx_sorted:
        y, x = int(flat_i // W), int(flat_i % W)
        if not allowed[y, x] or f[y, x] <= 0:
            continue
        cy, cx = int(y // cell), int(x // cell)
        ok = True
        for gy in range(cy - 2, cy + 3):
            for gx in range(cx - 2, cx + 3):
                for (ky, kx) in grid_d.get((gy, gx), ()):
                    if (y - ky) ** 2 + (x - kx) ** 2 < r2:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            out[y, x] = True
            chosen += 1
            grid_d.setdefault((cy, cx), []).append((y, x))
            if chosen >= target_n:
                break
    return out


def greedy_maxcover(field, foot, max_n, catalogue=None, cat_buffer_px=2, verbose=False):
    """Lazy-greedy max-coverage. Use only for small budgets (≤ 5k) or primaries."""
    f = np.asarray(field, np.float32) * foot
    H, W = f.shape
    C = np.zeros(f.shape, np.float32)
    chosen = np.zeros(f.shape, bool)
    blocked = _catalogue_block(catalogue, foot, cat_buffer_px)
    gain = np.zeros(f.shape, np.float32)
    for dy, dx, kk in zip(_DY, _DX, _KW):
        gain += kk * _shift(f, dy, dx)
    gain *= foot
    n_taken = 0
    for _ in range(max_n):
        g = np.where(chosen | blocked | ~foot, -np.inf, gain)
        i = int(np.argmax(g))
        v = float(g.reshape(-1)[i])
        if not np.isfinite(v) or v <= 0:
            break
        y, x = divmod(i, W)
        chosen[y, x] = True
        n_taken += 1
        for dy, dx, kk in zip(_DY, _DX, _KW):
            qy, qx = y + dy, x + dx
            if 0 <= qy < H and 0 <= qx < W:
                if kk > C[qy, qx]:
                    old = C[qy, qx]
                    C[qy, qx] = kk
                    wq = f[qy, qx]
                    if wq > 0:
                        for d2, e2, k2 in zip(_DY, _DX, _KW):
                            xy, xx = qy + d2, qx + e2
                            if 0 <= xy < H and 0 <= xx < W:
                                dgw = wq * (max(0.0, k2 - kk) - max(0.0, k2 - old))
                                if dgw:
                                    gain[xy, xx] += dgw
        gain[y, x] = -np.inf
        if verbose and n_taken % 1000 == 0:
            print(f"  greedy taken={n_taken} last_marg={v:.4f}")
    return chosen
