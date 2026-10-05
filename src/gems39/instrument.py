"""Local proxy instruments, calibrated against owner-reported live scores.

Why this module exists
----------------------
The only scored truth is DrivenData's private holdout, which is not accessible.
`scripts/forensics.py` measured, on 30 re-readable historical artifacts with
owner-reported scores, that:

  * emitted-pixel count        Spearman -0.686  (p<1e-4)  -> budget dominates
  * fraction of dots on the known catalogue      -0.598  (p=5e-4)
  * DTI against the *visible* catalogue          -0.344  (p=0.06)
  * DTI against off-catalogue SGMC faults        -0.080  (p=0.67)  -> NO signal

So the SGMC instrument used by sibling repositories is not supported by this
calibration set, and the catalogue-DTI is confounded by budget.  This module
implements the two instruments that survive and the partial-correlation test
that removes the budget confound.

Instrument H (catalogue-hidden, spatially blocked) - the promotion instrument
    Hide a random 25 % of the catalogue's connected components; evaluate DTI
    against the hidden components inside spatial blocks, with the visible
    catalogue (plus a small collar) removed from the evaluation domain.  This
    measures *generalisation to faults not seen*, which is the closest local
    analogue of "find a fault the catalogue does not contain".

Instrument B (budget-matched catalogue DTI)
    Same as catalogue DTI but only comparable between candidates emitted at
    the same budget; the historical correlation above is computed as a partial
    Spearman against log(n_pos).

All numbers this module produces are LOCAL PROXIES.  They are never reported as
leaderboard scores.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt, label

from .metric import ALPHA, BETA, EPS, kernel


@dataclass
class Cell:
    key: str
    sl: tuple
    active: np.ndarray      # evaluation domain (subset)
    truth: np.ndarray       # hidden truth inside `active`
    k_dg: np.ndarray        # kernel(distance to truth) inside `active`
    n_truth: int
    ty: np.ndarray = None   # truth row coords (local)
    tx: np.ndarray = None   # truth col coords (local)


@dataclass
class Holdout:
    foot: np.ndarray
    cat: np.ndarray
    blocks: np.ndarray
    cells: list = field(default_factory=list)


def make_blocks(foot, n_y=2, n_x=2):
    """Contiguous rectangular spatial blocks over the footprint."""
    yy, xx = np.nonzero(foot)
    H, W = foot.shape
    yq = np.quantile(yy, np.linspace(0, 1, n_y + 1)[1:-1])
    xq = np.quantile(xx, np.linspace(0, 1, n_x + 1)[1:-1])
    gy, gx = np.mgrid[0:H, 0:W]
    by = np.digitize(gy, yq)
    bx = np.digitize(gx, xq)
    b = (by * (n_x) + bx).astype(np.int32)
    b[~foot] = -1
    return b


def build_holdout(foot, cat, hide_frac=0.25, n_y=2, n_x=2, seeds=(11, 12, 13),
                  collar=2, domain_erode=0):
    """Hide `hide_frac` of catalogue components; return evaluation cells."""
    comp, nc = label(cat, structure=np.ones((3, 3), int))
    ids = np.arange(1, nc + 1)
    sizes = np.bincount(comp.ravel(), minlength=nc + 1)  # indexed by component id 1..nc
    blocks = make_blocks(foot, n_y, n_x)
    cells = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(ids)
        cum = np.cumsum(sizes[perm])
        k = int(np.searchsorted(cum, hide_frac * float(sizes[ids].sum()))) + 1
        hidden_ids = perm[:min(k, perm.size)]
        hidden = np.isin(comp, hidden_ids)
        visible = cat & ~hidden
        vis_buf = binary_dilation(visible, iterations=collar)
        for b in np.unique(blocks[foot]):
            q = blocks == b
            rows = np.flatnonzero(q.any(axis=1))
            cols = np.flatnonzero(q.any(axis=0))
            sl = (slice(int(rows[0]), int(rows[-1]) + 1),
                  slice(int(cols[0]), int(cols[-1]) + 1))
            qsub = q[sl]
            active = (qsub & foot[sl] & ~vis_buf[sl]).astype(bool)
            truth = (hidden[sl] & active).astype(bool)
            n_t = int(truth.sum())
            if n_t < 50:
                continue
            dg = distance_transform_edt(~truth)
            ty, tx = (np.nonzero(truth)[0].astype(np.int32),
                       np.nonzero(truth)[1].astype(np.int32))
            cells.append(Cell(key=f"s{seed}_b{int(b)}", sl=sl, active=active,
                              truth=truth, k_dg=kernel(dg).astype(np.float32), n_truth=n_t,
                              ty=ty, tx=tx))
    return Holdout(foot=foot, cat=cat, blocks=blocks, cells=cells)


def evaluate(pred, ctx, soft=False):
    """Evaluate a prediction (bool mask or soft float field) on the holdout.

    Returns per-cell DTI plus the pooled (TP/FP/FN summed over cells) DTI.
    Pooling is what the organizer does over the evaluation domain, so the
    pooled value is the headline; per-cell values feed the Wald SPRT.
    """
    p = np.asarray(pred, np.float64)
    p = np.where(np.isfinite(p), p, 0.0)
    if not soft:
        p = (p > 0.5).astype(np.float64)
    p = np.clip(p, 0.0, 1.0)
    TP = FP = FN = 0.0
    per_cell = {}
    for c in ctx.cells:
        ps = p[c.sl] * c.active
        if ps.sum() <= 0:
            per_cell[c.key] = dict(tp=0.0, fp=0.0, fn=float(c.n_truth),
                                   n_truth=c.n_truth, dti=0.0)
            FN += c.n_truth
            continue
        # max over the 3-px neighbourhood weighted by the triangular kernel
        best = np.zeros(c.n_truth, np.float64)
        ty, tx = c.ty, c.tx
        r = 3
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                kk = float(kernel(np.hypot(dy, dx)))
                if kk <= 0:
                    continue
                ny = np.clip(ty + dy, 0, ps.shape[0] - 1)
                nx = np.clip(tx + dx, 0, ps.shape[1] - 1)
                np.maximum(best, ps[ny, nx] * kk, out=best)
        tp = float(best.sum())
        fn = float(c.n_truth) - tp
        py, px = np.nonzero(ps > 0)
        fp = float(((1.0 - c.k_dg[py, px]) * ps[py, px]).sum())
        per_cell[c.key] = dict(tp=tp, fp=fp, fn=fn, n_truth=c.n_truth,
                               dti=float(tp / (tp + ALPHA * fp + BETA * fn + EPS)))
        TP += tp
        FP += fp
        FN += fn
    pooled = float(TP / (TP + ALPHA * FP + BETA * FN + EPS))
    return dict(pooled_dti=pooled, tp=TP, fp=FP, fn=FN, cells=per_cell)


def partial_spearman(x, y, z):
    """Spearman correlation between rank(x) and rank(y) with rank(z) removed."""
    from scipy.stats import spearmanr, linregress
    rx = _rank(x)
    ry = _rank(y)
    rz = _rank(z)
    rz = (rz - rz.mean()) / (rz.std() + 1e-12)
    rx_res = rx - (rx.mean() + _beta(rz, rx) * rz)
    ry_res = ry - (ry.mean() + _beta(rz, ry) * rz)
    rho, pv = spearmanr(rx_res, ry_res)
    return float(rho), float(pv)


def _beta(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / ((a * a).sum() + 1e-12))


def _rank(v):
    v = np.asarray(v, float)
    order = np.argsort(np.argsort(v))
    return order.astype(float)
