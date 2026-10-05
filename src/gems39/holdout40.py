"""Spatially-blocked holdout with two independent, prevalence-matched instruments.

Instrument I1 -- hidden catalogue components
    Whole connected components of the supplied USGS/INGENIOUS catalogue are
    hidden inside a block.  The emulation of the organizer's protocol is exact:
    the *visible* catalogue is masked (``known``), the hidden component is
    truth, and predictions are scored on ``active & ~known``.  DrivenData staff:
    "Pixels corresponding to known USGS/INGENIOUS faults are masked / excluded
    from evaluation, so they do not count towards penalty terms" and
    "Re-evaluation will also mask/exclude the existing USGS/INGENIOUS faults"
    (https://community.drivendata.org/t/11516/2).  This measures "can the
    detector find a mapped fault blind?".

Instrument I2 -- SGMC off-catalogue faults, prevalence matched
    USGS State Geologic Map Compilation fault pixels more than 300 m from the
    scored catalogue.  Staff: a "new fault" is "any fault pixel not already
    captured by USGS/INGENIOUS" and "can include newly mapped geometry of an
    existing fault system" (https://community.drivendata.org/t/11536/2).  I2 is
    an official, independent realisation of exactly that population.

    I2 is *prevalence matched* to the calibrated |G| from ``calibrate.py``.
    Unmatched it carries 62 122 truth pixels where the scored population has
    ~14 100, so a dot field's hit rate on I2 would be inflated roughly 4x and
    could not be compared with the calibrated live numbers.  Matching drops
    whole connected components -- never partial traces.

Memory discipline: this runs on a 3.9 GB sandbox against a 12.3 Mpx grid, so
every per-cell array is a *crop* of the block's bounding box, and the truth
kernel field is computed on demand inside ``dti_cell`` rather than stored.  An
earlier draft stored full-grid block masks and was killed by the OOM reaper
(exit 137); that is why this file is written the way it is.

Leakage rule: H40-F (the supervised channel) contributes only OUT-OF-FOLD
predictions with respect to the same block grid, so no channel is trained on a
block it is scored in.  ``leakage_probe`` measures the residual leak of a
candidate against the evaluation truth and the build gate refuses to promote
anything above a pre-declared threshold.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt, label

from .metric import ALPHA, BETA, EPS, kernel

STRUCT8 = np.ones((3, 3), int)


def blocks(shape, n_y=4, n_x=6):
    """Deterministic rectangular block ids (row-major) and their slices."""
    H, W = shape
    ys = np.linspace(0, H, n_y + 1).astype(int)
    xs = np.linspace(0, W, n_x + 1).astype(int)
    bid = np.full(shape, -1, np.int16)
    boxes = {}
    k = 0
    for i in range(n_y):
        for j in range(n_x):
            sl = (slice(ys[i], ys[i + 1]), slice(xs[j], xs[j + 1]))
            bid[sl] = k
            boxes[k] = sl
            k += 1
    return bid, boxes, k


def prevalence_match(mask, target_px, seed=20261005, min_comp_px=6):
    """Drop whole connected components until the pixel count is ~target_px."""
    comp, n = label(mask, structure=STRUCT8)
    if n == 0:
        return mask.copy(), 0
    sizes = np.bincount(comp.ravel(), minlength=n + 1)
    ids = np.arange(1, n + 1)
    ids = ids[sizes[ids] >= min_comp_px]
    if ids.size == 0:
        return mask.copy(), int(mask.sum())
    rng = np.random.default_rng(seed)
    perm = rng.permutation(ids)
    cum = np.cumsum(sizes[perm])
    k = min(int(np.searchsorted(cum, target_px)) + 1, perm.size)
    lut = np.zeros(n + 1, bool)
    lut[perm[:k]] = True
    out = lut[comp]
    return out, int(out.sum())


@dataclass
class Cell:
    key: str
    instrument: str
    bbox: tuple                 # (slice_row, slice_col) into the full grid
    active_sub: np.ndarray      # bool crop: pixels that are scored in this cell
    truth_sub: np.ndarray       # bool crop: truth pixels in this cell
    known_sub: np.ndarray       # bool crop: pixels masked out of scoring here
    n_truth: int = 0
    n_active: int = 0


@dataclass
class Holdout40:
    foot: np.ndarray
    catalogue: np.ndarray
    active: np.ndarray
    sgmc_off: np.ndarray
    sgmc_matched: np.ndarray
    bid: np.ndarray
    boxes: dict
    cells: list = field(default_factory=list)
    n_blocks: int = 0
    collar_px: int = 10
    erode_px: int = 6
    info: dict = field(default_factory=dict)


def load_holdout40(ddir: Path, target_truth_px=None, n_y=4, n_x=6, seeds=(31, 32),
                   collar_px=10, erode_px=6, sgmc_min_dist_px=3.0,
                   hidden_fraction=0.25):
    from .external import sgmc_off_catalogue

    with rasterio.open(ddir / "sample_submission.tif") as s:
        foot = np.isfinite(s.read(1))
    with rasterio.open(ddir / "labels.tif") as s:
        cat = (s.read(1) >= 1) & foot
    active = foot & ~cat
    sgmc_off, n_off = sgmc_off_catalogue(ddir, cat, foot, sgmc_min_dist_px)
    if target_truth_px:
        matched, n_m = prevalence_match(sgmc_off, target_truth_px)
    else:
        matched, n_m = sgmc_off, int(sgmc_off.sum())

    bid, boxes, nb = blocks(foot.shape, n_y, n_x)
    comp, ncomp = label(cat, structure=STRUCT8)
    comp = comp.astype(np.int32)
    comp_size = np.bincount(comp.ravel(), minlength=ncomp + 1)

    cells = []
    for seed in seeds:
        for b in range(nb):
            sl = boxes[b]
            q = np.zeros(foot.shape, bool)
            q[sl] = True
            interior = binary_erosion(q, iterations=erode_px)
            dom_full = interior & foot
            # DEFECT FIXED 2026-10-05: the collar was dilate(q, 10) & foot, which
            # CONTAINS q, so every component inside the block also "touched the
            # collar" and the candidate set was empty for all 24 blocks -- I1 had
            # zero truth pixels and silently reported DTI 0.0000.  The collar must
            # be the ANNULUS around the scored interior.
            col_full = binary_dilation(q, iterations=collar_px) & ~interior & foot
            if int(dom_full.sum()) == 0:
                del q, interior, dom_full, col_full
                continue
            # ---- I1: hide whole catalogue components that live inside the block
            cat_b = comp[(bid == b) & cat]
            in_blk = np.unique(cat_b)
            in_blk = in_blk[in_blk > 0]
            touch = np.unique(comp[col_full & cat])
            touch = touch[touch > 0]
            cand = np.setdiff1d(in_blk, touch, assume_unique=False)
            if cand.size:
                rng = np.random.default_rng(100_000 * (seed + 1) + b)
                perm = rng.permutation(cand)
                cum = np.cumsum(comp_size[perm])
                k = min(int(np.searchsorted(cum, hidden_fraction *
                                            float(comp_size[cand].sum()))) + 1, perm.size)
                lut = np.zeros(ncomp + 1, bool)
                lut[perm[:k]] = True
                hidden = lut[comp]
            else:
                hidden = np.zeros(foot.shape, bool)
            known_i1 = cat & ~hidden
            truth_i1 = hidden & dom_full
            active_i1 = (foot & ~known_i1) & dom_full
            cells.append(_crop_cell(f"I1_s{seed}_b{b:02d}", "I1", sl,
                                    active_i1, truth_i1, known_i1))
            del hidden, known_i1, truth_i1, active_i1
            # ---- I2: prevalence-matched SGMC off-catalogue ----
            truth_i2 = matched & dom_full
            cells.append(_crop_cell(f"I2_s{seed}_b{b:02d}", "I2", sl,
                                    active & dom_full, truth_i2, cat))
            del truth_i2
            del q, interior, dom_full, col_full

    del comp
    info = dict(n_blocks=nb, seeds=list(seeds), sgmc_off_px=int(n_off),
                sgmc_matched_px=int(n_m), target_truth_px=target_truth_px,
                collar_px=collar_px, erode_px=erode_px, hidden_fraction=hidden_fraction,
                catalogue_px=int(cat.sum()), active_px=int(active.sum()),
                footprint_px=int(foot.sum()), n_cells=len(cells),
                n_i1=sum(1 for c in cells if c.instrument == "I1"),
                n_i2=sum(1 for c in cells if c.instrument == "I2"))
    return Holdout40(foot=foot, catalogue=cat, active=active, sgmc_off=sgmc_off,
                     sgmc_matched=matched, bid=bid, boxes=boxes, cells=cells,
                     n_blocks=nb, collar_px=collar_px, erode_px=erode_px, info=info)


def _crop_cell(key, instr, sl, active_full, truth_full, known_full):
    a = np.ascontiguousarray(active_full[sl])
    t = np.ascontiguousarray(truth_full[sl])
    k = np.ascontiguousarray(known_full[sl])
    return Cell(key=key, instrument=instr, bbox=sl, active_sub=a, truth_sub=t,
                known_sub=k, n_truth=int(t.sum()), n_active=int(a.sum()))


def dti_cell(mask, cell):
    """Official DTI restricted to one cell, with that cell's own known mask.

    The truth kernel field is computed here rather than stored, which is what
    keeps 48 cells inside the memory budget.
    """
    p = np.asarray(mask[cell.bbox], bool) & cell.active_sub
    n_p = int(p.sum())
    n_t = cell.n_truth
    if n_t == 0:
        return dict(dti=0.0, tp=0.0, fp=float(n_p), fn=0.0, n_truth=0, n_pred=n_p,
                    truth_covered=0, truth_coverage=0.0, dot_hit_rate=0.0)
    if n_p == 0:
        return dict(dti=0.0, tp=0.0, fp=0.0, fn=float(n_t), n_truth=n_t, n_pred=0,
                    truth_covered=0, truth_coverage=0.0, dot_hit_rate=0.0)
    dp = distance_transform_edt(~p)
    d_at = dp[cell.truth_sub]
    tp = float(kernel(d_at).sum())
    covered = int((d_at <= 3.0).sum())
    fn = float(n_t) - tp
    dg = distance_transform_edt(~cell.truth_sub)
    ys, xs = np.nonzero(p)
    d_hits = int((dg[ys, xs] <= 3.0).sum())
    fp = float((1.0 - kernel(dg[p])).sum())
    return dict(dti=float(tp / (tp + ALPHA * fp + BETA * fn + EPS)), tp=tp, fp=fp,
                fn=fn, n_truth=n_t, n_pred=n_p,
                truth_covered=covered,
                truth_coverage=float(covered / max(n_t, 1)),
                dot_hit_rate=float(d_hits / max(n_p, 1)))


def instrument_truth(ctx, instr):
    """Full-grid truth for one instrument (I1 = union of every hidden component)."""
    if instr == "I2":
        return ctx.sgmc_matched
    t = np.zeros(ctx.foot.shape, bool)
    for cell in ctx.cells:
        if cell.instrument == "I1":
            t[cell.bbox] |= cell.truth_sub
    return t


def evaluate(mask, ctx, pooled=True):
    """Per-cell evaluation on both instruments, plus optional pooled DTI."""
    mask = np.asarray(mask, bool) & ctx.foot
    out = {}
    for instr in ("I1", "I2"):
        per = {c.key: dti_cell(mask, c) for c in ctx.cells if c.instrument == instr}
        dt = np.array([v["dti"] for v in per.values()]) if per else np.array([0.0])
        hr = np.array([v["dot_hit_rate"] for v in per.values()]) if per else np.array([0.0])
        rec = dict(cells=per, cell_mean=float(dt.mean()) if dt.size else 0.0,
                   n_cells=len(per),
                   cell_mean_dot_hit_rate=float(hr.mean()) if hr.size else 0.0)
        if pooled:
            truth = instrument_truth(ctx, instr)
            # I2: the whole catalogue is masked.  I1: only the *visible* part is
            # masked -- the hidden components ARE the truth, so masking the full
            # catalogue here would empty the truth set and report DTI 0.
            known = ctx.catalogue if instr == "I2" else (ctx.catalogue & ~truth)
            domain = ctx.foot & ~known
            from .metric import dti_binary
            r = dti_binary(mask, truth, valid=ctx.foot, known=known)
            dp = distance_transform_edt(~mask) if mask.any() else None
            hits = int((dp[truth & domain] <= 3.0).sum()) if dp is not None else 0
            # DEFECT FIXED 2026-10-05: this used to report
            #     hit_rate = (truth pixels within R of a dot) / (number of dots)
            # which divides a truth-side count by a prediction-side count.  It is
            # neither a hit rate nor a coverage, and it was the quantity the
            # instrument->live transfer factor was computed from.  Both real
            # quantities are now reported separately and named for what they are:
            #   truth_covered   = truth pixels within R of some dot  (0..|G|)
            #   truth_coverage  = truth_covered / |G|                (0..1)
            #   dot_hit_rate    = dots within R of some truth / n_dots (0..1)
            n_dots = int(mask.sum())
            ys, xs = np.nonzero(mask) if n_dots else (np.zeros(0, int),) * 2
            if n_dots:
                dg = distance_transform_edt(~(truth & domain))
                d_hits = int((dg[ys, xs] <= 3.0).sum())
                del dg
            else:
                d_hits = 0
            rec.update(pooled_dti=r["dti"], pooled_tp=r["tp"], pooled_fp=r["fp"],
                       pooled_fn=r["fn"], pooled_n_truth=r["n_truth"],
                       truth_covered=hits,
                       truth_coverage=float(hits / max(int((truth & domain).sum()), 1)),
                       dot_hit_rate=float(d_hits / max(n_dots, 1)),
                       n_pred=n_dots)
            del dp, truth, domain
        out[instr] = rec
    out["emitted_pixels"] = int(mask.sum())
    out["on_catalogue_pixels"] = int((mask & ctx.catalogue).sum())
    return out


def fold_score_delta(cell, cand, anchor, s_anchor):
    """First-order change in the POOLED score contributed by one fold.

    The pooled metric is ``s = TP / D`` with ``D = 0.2 TP + 0.2 FP + 0.8 |G|``.
    Its first-order variation is

        ds = (dTP - s * dD) / D,   dD = 0.2 dTP + 0.2 dFP
           = [ dTP (1 - 0.2 s) - 0.2 s dFP ] / D

    ``D`` is a positive constant shared by every fold, so the sign of ``ds`` is
    the sign of ``dTP (1 - 0.2 s) - 0.2 s dFP`` -- and that expression is
    ADDITIVE over folds.  Summing it over all folds therefore reproduces the
    first-order change in the pooled score exactly.

    This is why the SPRT uses this statistic rather than "did the fold's own DTI
    go up".  Per-fold DTI is a ratio: a candidate can win the pooled metric while
    losing a majority of folds (by concentrating gains where truth is dense), and
    the fold-level DTI test would then reject a genuine pooled improvement.  The
    first run of the selection stage hit exactly that, so the statistic was
    changed -- and the change is declared here rather than buried.

    Returns (delta, dTP, dFP, dti_cand, dti_anchor, fold_wins).
    """
    a = dti_cell(anchor, cell)
    c = dti_cell(cand, cell)
    dTP = c["tp"] - a["tp"]
    dFP = c["fp"] - a["fp"]
    delta = dTP * (1.0 - ALPHA * s_anchor) - ALPHA * s_anchor * dFP
    return dict(delta=float(delta), dTP=float(dTP), dFP=float(dFP),
                dti_cand=c["dti"], dti_anchor=a["dti"],
                fold_dti_win=bool(c["dti"] > a["dti"]),
                n_truth=cell.n_truth)


def leakage_probe(mask, ctx, truth):
    """Fraction of emitted pixels that sit on the evaluation truth itself.

    A detector *built from* the truth returns ~1.0 here.  Anything above the
    pre-declared gate means contamination and the candidate must not be
    promoted.
    """
    m = np.asarray(mask, bool) & ctx.foot
    if m.sum() == 0:
        return 0.0
    return float((m & truth).sum() / m.sum())
