"""Distance-Weighted Tversky Index (DTI) for DOE GEMS #306.

Parameters (verified against DrivenData page 967 and organizer forum posts):
  alpha = 0.2 (FP weight), beta = 0.8 (FN weight), R = 300 m = 3 px (100 m grid),
  triangular kernel k(d) = max(1 - d/R, 0).

References (manual review):
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric
  https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
R_PX = 3.0
EPS = 1e-12


def kernel(d, radius=R_PX):
    return np.maximum(1.0 - np.asarray(d, np.float64) / radius, 0.0)


def _offsets(radius=R_PX):
    r = int(np.ceil(radius))
    out = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            d = np.hypot(dy, dx)
            k = float(kernel(d, radius))
            if k > 0:
                out.append((dy, dx, k))
    return out


_OFFS = _offsets()


def dti_components(pred, truth, valid=None, known=None):
    """Return dict with TP_w, FP_w, FN_w, dti.  Masks are boolean; pred may be soft [0,1]."""
    pred = np.asarray(pred, np.float64)
    truth = np.asarray(truth, bool)
    if pred.shape != truth.shape:
        raise ValueError("shape mismatch")
    valid = np.ones(pred.shape, bool) if valid is None else np.asarray(valid, bool)
    known = np.zeros(pred.shape, bool) if known is None else np.asarray(known, bool)
    active = valid & ~known
    p = np.where(active & np.isfinite(pred), pred, 0.0).astype(np.float64)
    g = active & truth
    H, W = p.shape
    yy, xx = np.nonzero(g)
    n = int(yy.size)
    if n == 0:
        return dict(tp=0.0, fp=float(p.sum()), fn=0.0, n_truth=0, dti=0.0)
    credit = np.zeros(n, np.float64)
    for dy, dx, k in _OFFS:
        ny, nx = yy + dy, xx + dx
        ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W)
        credit[ok] = np.maximum(credit[ok], p[ny[ok], nx[ok]] * k)
    tp = float(credit.sum())
    fn = float(n) - tp
    d = distance_transform_edt(~g)
    fp = float((p * (1.0 - kernel(d))).sum())
    dti = tp / (tp + ALPHA * fp + BETA * fn + EPS)
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=float(dti))


def dti_binary(pred_bool, truth, valid=None, known=None):
    pred_bool = np.asarray(pred_bool, bool)
    truth = np.asarray(truth, bool)
    valid = np.ones(pred_bool.shape, bool) if valid is None else np.asarray(valid, bool)
    known = np.zeros(pred_bool.shape, bool) if known is None else np.asarray(known, bool)
    active = valid & ~known
    p = pred_bool & active
    g = truth & active
    n = int(g.sum())
    if n == 0:
        return dict(tp=0.0, fp=float(p.sum()), fn=0.0, n_truth=0, dti=0.0)
    if not p.any():
        return dict(tp=0.0, fp=0.0, fn=float(n), n_truth=n, dti=0.0)
    dp = distance_transform_edt(~p)
    tp = float(kernel(dp[g]).sum())
    fn = float(n) - tp
    dg = distance_transform_edt(~g)
    fp = float((1.0 - kernel(dg[p])).sum())
    dti = tp / (tp + ALPHA * fp + BETA * fn + EPS)
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=float(dti))


def dti(tp, fp, fn):
    return float(tp / (tp + ALPHA * fp + BETA * fn + EPS))
