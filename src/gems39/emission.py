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
    """Return the exact/near-catalogue cells excluded from emission.

    A zero-pixel buffer means exclude the catalogue cells themselves; it does
    not disable masking. The competition organizer clarified that known-fault
    pixels are excluded from scoring, but did not say that a 200 m halo is
    masked. Buffering neighboring pixels can therefore throw away useful DTI
    credit for a distinct new fault.
    """
    if catalogue is None:
        return np.zeros(foot.shape, bool)
    if cat_buffer_px < 0:
        raise ValueError("cat_buffer_px must be non-negative")
    known = np.asarray(catalogue, bool) & np.asarray(foot, bool)
    if known.shape != foot.shape:
        raise ValueError("catalogue and footprint shape mismatch")
    if not known.any():
        return np.zeros(foot.shape, bool)
    d = distance_transform_edt(~known)
    return d <= float(cat_buffer_px)


def emit_fast(field, foot, target_n, min_dist=2.7, catalogue=None, cat_buffer_px=0):
    """Best-first Poisson-disk emission with a matched positive-pixel budget.

    Pixels are ordered by detector score and accepted if they are at least
    ``min_dist`` pixels from previously accepted points. The top-k pool expands
    until the requested budget is met or every allowed pixel has been examined;
    unlike the earlier fixed 6x pool, this does not silently underfill merely
    because high-ranked points clustered too tightly.
    """
    foot = np.asarray(foot, bool)
    f = np.asarray(field, np.float32)
    if f.shape != foot.shape:
        raise ValueError("field and footprint shape mismatch")
    if target_n < 0:
        raise ValueError("target_n must be non-negative")
    if min_dist < 0:
        raise ValueError("min_dist must be non-negative")
    if target_n == 0:
        return np.zeros(foot.shape, bool)
    blocked = _catalogue_block(catalogue, foot, cat_buffer_px)
    allowed = foot & ~blocked & np.isfinite(f)
    if not allowed.any():
        return np.zeros(foot.shape, bool)
    score = np.where(allowed, f, -np.inf).astype(np.float32, copy=False)
    flat = score.reshape(-1)
    allowed_n = int(allowed.sum())
    target_n = min(int(target_n), allowed_n)
    H, W = f.shape
    r2 = float(min_dist) * float(min_dist)
    cell = max(1.0, float(min_dist))
    r = int(np.ceil(min_dist))
    pool_n = min(allowed_n, max(target_n, 6 * target_n))

    while True:
        if pool_n == allowed_n:
            idx = np.flatnonzero(np.isfinite(flat))
        else:
            idx = np.argpartition(flat, -pool_n)[-pool_n:]
        idx = idx[np.argsort(-flat[idx], kind="stable")]
        grid_d = {}
        out = np.zeros(f.shape, bool)
        chosen = 0
        for flat_i in idx:
            y, x = divmod(int(flat_i), W)
            if flat[y * W + x] <= 0:
                continue
            cy, cx = int(y // cell), int(x // cell)
            ok = True
            for gy in range(cy - 2, cy + 3):
                for gx in range(cx - 2, cx + 3):
                    for ky, kx in grid_d.get((gy, gx), ()):
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
                    return out
        if pool_n == allowed_n:
            return out
        pool_n = min(allowed_n, max(pool_n + 1, pool_n * 2))


def emit_sequence(field, foot, max_n, min_dist=2.7, catalogue=None, cat_buffer_px=0):
    """Best-first Poisson-disk emission that returns the *acceptance order*.

    Same greedy rule as :func:`emit_fast` (pixels visited in descending field
    order, accepted when at least ``min_dist`` from every already-accepted
    point), but instead of stopping at a budget it records the ordered list of
    accepted pixels up to ``max_n``.  Because the rule is a prefix-stable
    greedy, ``emit_prefix(seq, k)`` is exactly what ``emit_fast`` would return
    for budget ``k`` -- so a whole budget sweep costs one emission instead of
    one per budget.  The candidate pool is always the full allowed set, which
    is what makes the prefixes identical across budgets.
    """
    foot = np.asarray(foot, bool)
    f = np.asarray(field, np.float32)
    if f.shape != foot.shape:
        raise ValueError("field and footprint shape mismatch")
    if max_n <= 0:
        return np.zeros((0, 2), np.int32)
    blocked = _catalogue_block(catalogue, foot, cat_buffer_px)
    allowed = foot & ~blocked & np.isfinite(f) & (f > 0)
    ys, xs = np.nonzero(allowed)
    if ys.size == 0:
        return np.zeros((0, 2), np.int32)
    order = np.argsort(-f[ys, xs], kind="stable")
    ys, xs = ys[order], xs[order]
    H, W = f.shape
    r2 = float(min_dist) * float(min_dist)
    cell = max(1.0, float(min_dist))
    grid_d: dict = {}
    out = []
    for i in range(ys.size):
        y, x = int(ys[i]), int(xs[i])
        cy, cx = int(y // cell), int(x // cell)
        ok = True
        for gy in range(cy - 2, cy + 3):
            for gx in range(cx - 2, cx + 3):
                for ky, kx in grid_d.get((gy, gx), ()):
                    if (y - ky) ** 2 + (x - kx) ** 2 < r2:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            out.append((y, x))
            grid_d.setdefault((cy, cx), []).append((y, x))
            if len(out) >= max_n:
                break
    return np.asarray(out, np.int32).reshape(-1, 2)


def emit_prefix(seq, k, shape):
    """Rasterise the first ``k`` points of an :func:`emit_sequence` result."""
    m = np.zeros(shape, bool)
    if len(seq):
        s = seq[: int(k)]
        if s.size:
            m[s[:, 0], s[:, 1]] = True
    return m


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
