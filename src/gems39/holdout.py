"""Spatially-blocked holdout for candidate evaluation.

Splits the footprint into 4 quadrants (NW/NE/SW/SE) at the median valid row/col,
erodes the domain by a collar to prevent bleed, and runs 2 draws of 20% label
component hiding per quadrant. Returns per-fold DTI for SPRT consumption.

No live-score peeking is possible here: the holdout uses only `labels.tif`
(visible catalogue) and is run before any submission-slot decision.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import rasterio
from scipy.ndimage import (binary_dilation, binary_erosion, distance_transform_edt,
                           label)
from .metric import ALPHA, BETA, EPS, kernel, dti_binary

FOLD_NAMES = ("NW", "NE", "SW", "SE")
COLLAR_PX = 15
DOMAIN_ERODE = 12
NEAR_LABEL_PX = 3


def quadrant_ids(footprint):
    footprint = np.asarray(footprint, bool)
    yy, xx = np.nonzero(footprint)
    ym, xm = int(np.median(yy)), int(np.median(xx))
    H, W = footprint.shape
    gy, gx = np.ogrid[:H, :W]
    q = np.full((H, W), -1, np.int8)
    q[(gy < ym) & (gx < xm) & footprint] = 0
    q[(gy < ym) & (gx >= xm) & footprint] = 1
    q[(gy >= ym) & (gx < xm) & footprint] = 2
    q[(gy >= ym) & (gx >= xm) & footprint] = 3
    return q


@dataclass
class Cell:
    key: str
    bbox: tuple[slice, slice]
    active_sub: np.ndarray
    truth_sub: np.ndarray
    k_dg_sub: np.ndarray
    n_truth: int


@dataclass
class Holdout:
    foot: np.ndarray
    labels: np.ndarray
    near_known: np.ndarray
    sgmc_off: np.ndarray
    quad: np.ndarray
    quad_bboxes: dict
    cells: list
    # SGMC off-catalogue (independent from the visible labels)
    sgmc_truth_coords: tuple
    sgmc_k_dg: np.ndarray
    sgmc_n_truth: int


def load_holdout(ddir: Path, hide_frac=0.20, seeds=(20, 21)):
    with rasterio.open(ddir / "sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1))
    with rasterio.open(ddir / "labels.tif") as ds:
        lbl = ds.read(1)
    labels = (lbl >= 1) & foot
    with rasterio.open(ddir / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgm = ds.read(1) > 0
    near = binary_dilation(labels, iterations=NEAR_LABEL_PX)
    sgmc_off = sgm & ~labels & ~near & foot
    quad = quadrant_ids(foot)
    comp, n_comp = label(labels, structure=np.ones((3, 3), int))
    comp_size = np.bincount(comp.ravel(), minlength=n_comp + 1)
    all_ids = np.arange(1, n_comp + 1)
    quad_bboxes = {}
    collars = {}
    domains = {}
    for fold in range(4):
        q = quad == fold
        rows = np.flatnonzero(q.any(axis=1))
        cols = np.flatnonzero(q.any(axis=0))
        r0 = max(0, int(rows[0]) - 6)
        r1 = min(foot.shape[0], int(rows[-1]) + 7)
        c0 = max(0, int(cols[0]) - 6)
        c1 = min(foot.shape[1], int(cols[-1]) + 7)
        quad_bboxes[fold] = (slice(r0, r1), slice(c0, c1))
        collars[fold] = binary_dilation(q, iterations=COLLAR_PX) & foot
        domains[fold] = binary_erosion(q, iterations=DOMAIN_ERODE)
    cells = []
    for seed in seeds:
        for fold in range(4):
            rng = np.random.default_rng(10_000 * (seed + 1) + fold)
            q = quad == fold
            collar = collars[fold]
            domain = domains[fold]
            touch_collar = np.isin(all_ids, np.unique(comp[collar & labels]))
            in_test = np.isin(all_ids, np.unique(comp[q & labels]))
            test_ids = all_ids[in_test]
            train_ids = all_ids[~touch_collar]
            # pick components to hide
            def _pick(ids, target):
                if ids.size == 0:
                    return ids
                perm = rng.permutation(ids)
                cum = np.cumsum(comp_size[perm])
                k = int(np.searchsorted(cum, target)) + 1
                return perm[:min(k, perm.size)]
            hid_test = _pick(test_ids, hide_frac * float((labels & q).sum()))
            hid_train = _pick(train_ids, hide_frac * float(comp_size[train_ids].sum()))
            hidden_full = np.isin(comp, hid_test)
            hidden_train = np.isin(comp, hid_train)
            visible = labels & ~hidden_full & ~hidden_train
            sl = quad_bboxes[fold]
            active_sub = domain[sl] & ~visible[sl]
            truth_sub = (hidden_full[sl] & domain[sl]) & active_sub
            t_coords = np.nonzero(truth_sub)
            n_t = int(t_coords[0].size)
            dg_sub = distance_transform_edt(~truth_sub)
            k_dg_sub = kernel(dg_sub)
            cells.append(Cell(key=f"s{seed}_f{fold}", bbox=sl, active_sub=active_sub,
                              truth_sub=truth_sub, k_dg_sub=k_dg_sub, n_truth=n_t))
    # SGMC off-catalogue ground (independent test)
    sgmc_active = foot & ~labels
    sgmc_truth = sgmc_off & sgmc_active
    sgc = np.nonzero(sgmc_truth)
    n_s = int(sgc[0].size)
    k_s = kernel(distance_transform_edt(~sgmc_truth))
    return Holdout(foot=foot, labels=labels, near_known=near, sgmc_off=sgmc_off,
                   quad=quad, quad_bboxes=quad_bboxes, cells=cells,
                   sgmc_truth_coords=sgc, sgmc_k_dg=k_s, sgmc_n_truth=n_s)


def evaluate(mask, ctx):
    """Binary-mask evaluation across all cells and SGMC off-catalogue."""
    mask = np.asarray(mask, bool) & ctx.foot
    n_emit = int(mask.sum())
    on_cat = int((mask & ctx.labels).sum())
    off_mask = mask & ~ctx.labels
    dp_off = distance_transform_edt(~off_mask) if off_mask.any() else np.full(mask.shape, 999.0, np.float32)
    folds = {}
    for cell in ctx.cells:
        sl = cell.bbox
        p_sub = mask[sl] & cell.active_sub
        n_p = int(p_sub.sum())
        n_t = cell.n_truth
        if n_t == 0 or n_p == 0:
            folds[cell.key] = dict(dti=0.0, tp=0.0, fp=float(n_p), fn=float(n_t), n_truth=n_t)
            continue
        dp_sub = distance_transform_edt(~p_sub)
        d_at = dp_sub[cell.truth_sub]
        tp = float(kernel(d_at).sum())
        fn = float(n_t) - tp
        fp = float((1.0 - cell.k_dg_sub[p_sub]).sum())
        folds[cell.key] = dict(tp=tp, fp=fp, fn=fn, n_truth=n_t,
                               dti=float(tp / (tp + ALPHA * fp + BETA * fn + EPS)))
    if off_mask.any() and ctx.sgmc_n_truth > 0:
        tp_s = float(kernel(dp_off[ctx.sgmc_truth_coords]).sum())
        fn_s = float(ctx.sgmc_n_truth) - tp_s
        fp_s = float((1.0 - ctx.sgmc_k_dg[off_mask]).sum())
        dti_sg = float(tp_s / (tp_s + ALPHA * fp_s + BETA * fn_s + EPS))
    else:
        dti_sg = 0.0
    fold_dtis = np.array([v["dti"] for v in folds.values()])
    per_quad = {FOLD_NAMES[f]: float(np.mean([folds[f"s{s}_f{f}"]["dti"] for s in (20, 21)]))
                for f in range(4)}
    return dict(
        emitted_pixels=n_emit, on_catalogue_pixels=on_cat,
        catalogue_hidden_mean=float(fold_dtis.mean()),
        catalogue_hidden_per_quadrant=per_quad,
        catalogue_hidden_folds=folds,
        sgmc_off_catalogue_dti=dti_sg,
    )
