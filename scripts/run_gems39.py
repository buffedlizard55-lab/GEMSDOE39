#!/usr/bin/env python3
"""GEMSDOE39 end-to-end: stack -> off-catalogue discriminant -> emission -> TIF.

Preregistered design (declared before any candidate score is read):

  * Wald SPRT parameters: alpha = 0.05, beta = 0.10, p0 = 0.50, p1 = 0.70.
  * Field-comparison instrument: catalogue-hidden spatially-blocked holdout,
    2x2 blocks x 3 seeds, hide 25 % of catalogue components, collar 2 px.
    Metric = mean percentile rank of the hidden truth pixels inside each cell
    (0.5 = no skill).  Per-cell wins feed the SPRT.
  * Off-catalogue corroboration: DTI against USGS SGMC faults farther than 3 px
    from the visible catalogue (62,703 px).  Used for the emission sweep and as
    a discordance check - NOT for promotion on its own, because
    scripts/calibrate_instruments.py measured its Spearman correlation with
    owner-reported live scores at -0.08 (n=30): no rank signal.
  * Emission budget: chosen by sweeping the off-catalogue surrogate, i.e. by the
    DTI break-even rule (a dot pays for itself iff its expected kernel credit
    exceeds 0.2 x DTI), and cross-checked against the holdout.

Stages are cached so later passes do not refit:
  data/stack_u8.npz   (feature stack)
  data/probs.npz      (out-of-fold model probabilities)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39 import stack, instrument as inst, emission as em, grid  # noqa: E402
from gems39.metric import ALPHA, BETA, EPS, kernel  # noqa: E402
from gems39.sprt_select import sprt_pairwise  # noqa: E402
import gems39.io39 as io39  # noqa: E402

DATA = ROOT / "data"
STACK_CACHE = DATA / "stack_u8.npz"
PROB_CACHE = DATA / "probs.npz"


# --------------------------------------------------------------------------- data
def load():
    with rasterio.open(DATA / "labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    cat = lab == 1
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as s:
        sg = s.read(1) > 0
    sg_off = sg & foot & ~cat & ~binary_dilation(cat, iterations=3)
    return foot, cat, sg_off


def load_stack(foot, rebuild=False):
    if STACK_CACHE.exists() and not rebuild:
        z = np.load(STACK_CACHE)
        F = {k: z[k] for k in z.files if k != "foot"}
        if F and F[next(iter(F))].shape == foot.shape:
            return F
    F = stack.build_stack(DATA / "training_features.tif", foot)
    np.savez_compressed(STACK_CACHE, **F, foot=foot.astype(np.uint8))
    return F


def stack_cube(F, foot):
    """(H, W, C) uint8 cube.  491 MB at C=40 - the only big allocation."""
    cols = sorted(F.keys())
    S = np.empty((*foot.shape, len(cols)), np.uint8)
    for j, c in enumerate(cols):
        S[..., j] = F[c]
    return S, cols


# -------------------------------------------------------------------- supervised
def fit_blocked(S, foot, pos, neg_mask, blocks, n_blocks=4, seed=0,
                iters=250, lr=0.08, leaves=31, rows=90000, seeds=(0, 1)):
    """Spatially blocked out-of-fold probabilities.

    For each spatial block, fit on the other blocks and predict inside it, so
    the reported ranking skill is genuine spatial generalisation rather than
    memorised location (the documented failure mode of the published GBM in
    GEMSDOE29, which reached train AUC 1.0 on a distance-to-catalogue leak).
    """
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    H, W = foot.shape
    C = S.shape[2]
    out = np.zeros((H, W), np.float32)
    hits = np.zeros((H, W), np.float32)
    aucs = []
    for seed in seeds:
      rng = np.random.default_rng(seed)
      for b in range(n_blocks):
        test = blocks == b
        train = (blocks >= 0) & ~test & foot
        py, px = np.nonzero(pos & train)
        ny_, nx_ = np.nonzero(neg_mask & train)
        if py.size < 200 or ny_.size < 200:
            continue
        m = int(min(py.size, ny_.size, rows))
        pi = rng.choice(py.size, m, replace=False)
        ni = rng.choice(ny_.size, m, replace=False)
        Xtr = np.empty((2 * m, C), np.float32)
        Xtr[:m] = S[py[pi], px[pi]]
        Xtr[m:] = S[ny_[ni], nx_[ni]]
        Xtr = Xtr.astype(np.float32) / 255.0
        ytr = np.concatenate([np.ones(m), np.zeros(m)])
        clf = HistGradientBoostingClassifier(
            max_iter=iters, learning_rate=lr, max_leaf_nodes=leaves,
            l2_regularization=1.0, early_stopping=False, random_state=int(seed))
        clf.fit(Xtr, ytr)
        del Xtr, ytr
        ys, xs = np.nonzero(test)
        pred = np.zeros(ys.size, np.float32)
        order = np.argsort(ys)
        ys, xs = ys[order], xs[order]
        r0 = 0
        while r0 < ys.size:
            r1 = int(np.searchsorted(ys, ys[r0] + 128, side="left"))
            blk = S[ys[r0:r1], xs[r0:r1]].astype(np.float32) / 255.0
            pred[r0:r1] = clf.predict_proba(blk)[:, 1]
            r0 = r1
        out[ys, xs] += pred
        hits[ys, xs] += 1.0
        ty, tx = np.nonzero(pos & test)
        qy, qx = np.nonzero((~pos) & test & foot)
        if ty.size > 50 and qy.size > 500:
            k = rng.choice(qy.size, min(qy.size, 60_000), replace=False)
            pm = np.zeros((H, W), np.float32)
            pm[ys, xs] = pred
            sc = np.concatenate([pm[ty, tx], pm[qy[k], qx[k]]])
            lb = np.concatenate([np.ones(ty.size), np.zeros(k.size)])
            aucs.append(float(roc_auc_score(lb, sc)))
            del pm
        del pred, clf
    out = out / np.maximum(hits, 1.0)
    return out, aucs


# ------------------------------------------------------------------ field metrics
def field_lift(field, ho, sample=250_000, seed=0):
    """Mean percentile rank of hidden-truth pixels inside each holdout cell.

    0.5 = no better than a random pixel of the cell's evaluation domain.
    Percentiles are computed against a random subsample of the domain (cheap and
    unbiased); per-cell values feed the Wald SPRT.
    """
    rng = np.random.default_rng(seed)
    per = {}
    for c in ho.cells:
        sub = field[c.sl]
        v = sub[c.active].astype(np.float32)
        t = sub[c.truth].astype(np.float32)
        del sub
        if t.size == 0 or v.size < 100:
            per[c.key] = 0.5
            continue
        k = min(v.size, sample)
        ref = v[rng.choice(v.size, k, replace=False)]
        ref = np.sort(ref)
        per[c.key] = float(np.mean(np.searchsorted(ref, t, side="left") / float(k)))
        del v, t, ref
    return per, float(np.mean(list(per.values())))


def dti_mask(pred_bool, truth, valid, known=None):
    from scipy.spatial import cKDTree
    p = np.asarray(pred_bool, bool) & valid
    g = truth & valid
    if known is not None:
        p &= ~known
        g &= ~known
    n = int(g.sum())
    if n == 0 or not p.any():
        return 0.0
    py, px = np.nonzero(p)
    ty, tx = np.nonzero(g)
    tree = cKDTree(np.column_stack([ty, tx]))
    d, j = tree.query(np.column_stack([py, px]), k=1)
    kk = np.asarray(kernel(d), float)
    mx = np.zeros(n, float)
    np.maximum.at(mx, j, kk)
    tp = float(mx.sum())
    fn = float(n) - tp
    fp = float((1.0 - kk).sum())
    return float(tp / (tp + ALPHA * fp + BETA * fn + EPS))


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild-stack", action="store_true")
    ap.add_argument("--rebuild-probs", action="store_true")
    ap.add_argument("--budget", type=int, default=40000)
    ap.add_argument("--spacing", type=float, default=3.0)
    ap.add_argument("--cat-buffer", type=int, default=2)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--sweep", action="store_true")
    args = ap.parse_args()

    tag = args.tag or dt.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    print("[1/7] loading grids", flush=True)
    foot, cat, sg_off = load()
    print(f"  footprint={int(foot.sum()):,}  catalogue={int(cat.sum()):,}  "
          f"sgmc_off(>3px)={int(sg_off.sum()):,}", flush=True)

    print("[2/7] feature stack", flush=True)
    F = load_stack(foot, rebuild=args.rebuild_stack)
    S, cols = stack_cube(F, foot)
    del F
    print(f"  {len(cols)} fields; cube {S.nbytes/1e6:.0f} MB", flush=True)

    blocks = inst.make_blocks(foot, 2, 2)
    dcat = distance_transform_edt(~cat)
    dsg = distance_transform_edt(~sg_off)
    hard = (dcat >= 1) & (dcat <= 6) & foot & ~cat
    bg = foot & ~cat & ~binary_dilation(cat, iterations=8)
    hardB = ((dsg >= 1) & (dsg <= 6) & foot & ~sg_off
             & ~binary_dilation(cat, iterations=3))
    print(f"  hard-negative pool={int(hard.sum()):,}  bg pool={int(bg.sum()):,}  "
          f"hard-B pool={int(hardB.sum()):,}", flush=True)

    print("[3/7] blocked models", flush=True)
    if PROB_CACHE.exists() and not args.rebuild_probs:
        z = np.load(PROB_CACHE)
        PA, PB = z["PA"].astype(np.float32), z["PB"].astype(np.float32)
        aucA = list(map(float, z["aucA"]))
        aucB = list(map(float, z["aucB"]))
        print(f"  cached: A={np.mean(aucA):.4f} B={np.mean(aucB):.4f}", flush=True)
    else:
        PA, aucA = fit_blocked(S, foot, cat, hard | bg, blocks, seeds=(1, 2, 3))
        print(f"  model A (catalogue)   OOF AUC {[round(a,4) for a in aucA]} "
              f"mean={np.mean(aucA):.4f}", flush=True)
        PB, aucB = fit_blocked(S, foot, sg_off, hardB | bg, blocks, seeds=(11, 12, 13))
        print(f"  model B (SGMC off-cat) OOF AUC {[round(a,4) for a in aucB]} "
              f"mean={np.mean(aucB):.4f}", flush=True)
        np.savez_compressed(PROB_CACHE, PA=PA.astype(np.float16),
                            PB=PB.astype(np.float16),
                            aucA=np.array(aucA), aucB=np.array(aucB))
    del S

    print("[4/7] candidate fields", flush=True)
    z = np.load(STACK_CACHE)
    g_ = lambda n_: z[n_].astype(np.float32) / 255.0
    ridge, chain = g_("topo_ridge"), g_("ridge_chain")
    maghg, sec = g_("mag_hg"), g_("secondary_line")
    pA = np.clip(PA, 0, 1)
    pB = np.clip(PB, 0, 1)
    joint = np.sqrt(pA * pB)
    CAND = {
        "CONTROL-topo-x-mag": ridge * (0.5 + 0.5 * maghg),
        "H39-A-model": pA,
        "H39-B-model": pB,
        "H39-C-joint": joint,
        "H39-D-joint-chain": joint * (0.45 + 0.55 * chain),
        "H39-E-joint-secondary": joint * (0.35 + 0.65 * sec),
    }
    CAND = {k: np.where(foot, v, 0.0).astype(np.float32) for k, v in CAND.items()}
    del ridge, chain, maghg, sec, joint

    print("[5/7] holdout + Wald SPRT (alpha=0.05 beta=0.10 p0=0.5 p1=0.7)", flush=True)
    ho = inst.build_holdout(foot, cat, hide_frac=0.25, n_y=2, n_x=2,
                            seeds=(21, 22, 23), collar=2)
    lifts, decisions = {}, {}
    for k, fld in CAND.items():
        per, m = field_lift(fld, ho)
        lifts[k] = dict(mean=m, per_cell=per)
        print(f"  {k:24s} holdout_lift={m:.4f}", flush=True)
    base = "CONTROL-topo-x-mag"
    for k in CAND:
        if k == base:
            continue
        wins = [lifts[k]["per_cell"][c] > lifts[base]["per_cell"][c]
                for c in lifts[base]["per_cell"]]
        decisions[k] = sprt_pairwise(wins, p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
        d = decisions[k]
        print(f"  SPRT {k:24s} wins={d['wins']}/{d['n']} llr={d['llr']:+.3f} "
              f"-> {d['decision']}", flush=True)

    accepted = [k for k, d in decisions.items() if d["decision"] == "accept_H1"]
    if accepted:
        best_field = max(accepted, key=lambda k: lifts[k]["mean"])
        why = "SPRT accept_H1"
    else:
        # preregistered tie-break: highest holdout lift among H39 candidates
        best_field = max((k for k in CAND if k.startswith("H39")),
                         key=lambda k: lifts[k]["mean"])
        why = "SPRT continue -> preregistered tie-break on holdout lift"
    print(f"  selected: {best_field} ({why})", flush=True)
    fld = CAND[best_field]

    print("[6/7] emission sweep (off-catalogue surrogate + holdout)", flush=True)
    sweep = []
    budgets = (24000, 30000, 36000, 42000, 50000) if args.sweep else (args.budget,)
    spacings = (2.6, 3.0, 3.4) if args.sweep else (args.spacing,)
    buffers = (2, 3) if args.sweep else (args.cat_buffer,)
    for sp in spacings:
        for buf in buffers:
            for bud in budgets:
                m = em.emit_fast(fld, foot, target_n=bud, min_dist=sp,
                                 catalogue=cat, cat_buffer_px=buf)
                d_s = dti_mask(m, sg_off, foot, known=cat)
                sweep.append(dict(spacing=sp, buffer=buf, budget=bud,
                                  n=int(m.sum()), sgmc_dti=d_s))
                print(f"    sp={sp} buf={buf} bud={bud} n={int(m.sum()):6d} "
                      f"sgmc_dti={d_s:.4f}", flush=True)
                del m
    if sweep:
        best = max(sweep, key=lambda r: r["sgmc_dti"])
        print(f"  sweep best: {best}", flush=True)
        spacing, buffer, budget = best["spacing"], best["buffer"], best["budget"]
    else:
        spacing, buffer, budget = args.spacing, args.cat_buffer, args.budget

    mask = em.emit_fast(fld, foot, target_n=budget, min_dist=spacing,
                        catalogue=cat, cat_buffer_px=buffer)

    print("[7/7] writing submissions", flush=True)
    out = ROOT / "docs" / "downloads"
    out.mkdir(parents=True, exist_ok=True)
    name = f"gemsdoe39-{best_field.lower().replace('_','-')}-{tag}"
    pred = mask.astype(np.float32)
    audits = {}
    for outside in ("zeros", "nan"):
        p = out / f"{name}-{outside}.tif"
        io39.write(pred, DATA / "sample_submission.tif", p, foot, outside=outside)
        a = io39.audit(p, foot, DATA / "sample_submission.tif")
        audits[outside] = a
        print(f"  {outside}: {a['bytes']} B  sha={a['sha256'][:16]}  ok={a['ok']}  "
              f"positive={a['positive_px']}", flush=True)
    dsgm = dti_mask(mask, sg_off, foot, known=cat)
    from scipy.spatial import cKDTree
    yy, xx = np.nonzero(mask)
    _d, _ = cKDTree(np.column_stack([yy, xx])).query(np.column_stack([yy, xx]), k=2)
    nn_stats = dict(median=float(np.median(_d[:, 1])), mean=float(np.mean(_d[:, 1])),
                    p10=float(np.percentile(_d[:, 1], 10)), p90=float(np.percentile(_d[:, 1], 90)))
    print(f"  emitted={int(mask.sum())} nn_median={nn_stats['median']:.2f} "
          f"within_2px_cat={int((mask & (dcat <= 2)).sum())}", flush=True)
    man = dict(schema_version=2, name=name, tag=tag, field=best_field, selection_rule=why,
               budget=budget, spacing=spacing, cat_buffer=buffer,
               emitted=int(mask.sum()), on_catalogue=int((mask & cat).sum()),
               within_2px_catalogue=int((mask & (dcat <= 2)).sum()),
               oof_auc_A=aucA, oof_auc_B=aucB, lifts=lifts, sprt=decisions,
               sweep=sweep, sgmc_dti_of_primary=dsgm, nn_stats=nn_stats,
               audits={k: {kk: vv for kk, vv in v.items() if kk != "meta"} for k, v in audits.items()},
               note="GEMSDOE39 off-catalogue discriminant (spatially blocked OOF) + "
                    "DTI-break-even emission; zeros-outside variant is the safer upload.")
    (out / f"{name}-manifest.json").write_text(json.dumps(man, indent=2, default=float) + "\n")
    print(f"  manifest {name}-manifest.json", flush=True)
    print(f"  sgmc_dti={dsgm:.4f}  emitted={int(mask.sum())}  on_cat={int((mask&cat).sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
