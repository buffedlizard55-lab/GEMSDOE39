#!/usr/bin/env python3
"""Calibrate candidate local proxy instruments against owner-reported live scores.

Reads registry/score_ledger.csv (owner-reported scores; NOT organizer receipts),
re-opens every artifact, and computes:

  H  catalogue-hidden spatially-blocked holdout DTI (pooled + per-cell)
  C  catalogue DTI (masked and unmasked)
  S  off-catalogue SGMC DTI

then reports raw and budget-adjusted (partial on log n) Spearman correlations
with the reported live score.

Output: registry/instrument_calibration.json
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39 import instrument as inst  # noqa: E402
from gems39.metric import ALPHA, BETA, EPS, kernel  # noqa: E402
from gems39.instrument import partial_spearman  # noqa: E402

DATA = ROOT / "data"
REF_TR = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def dti(pred, truth, valid, known=None):
    p = np.asarray(pred, float) * valid
    if known is not None:
        p = p * ~known
    g = truth & valid
    if known is not None:
        g = g & ~known
    n = int(g.sum())
    if n == 0:
        return 0.0
    if p.sum() <= 0:
        return 0.0
    from scipy.ndimage import distance_transform_edt
    dg = distance_transform_edt(~g)
    py, px = np.nonzero(p > 0)
    # weighted TP via max-kernel over the 3px neighbourhood
    best = np.zeros(n, float)
    ty, tx = np.nonzero(g)
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([ty, tx]))
    tpw = 0.0
    # per prediction pixel: assign to nearest truth; approximate max-credit by
    # grouping (exact for well separated traces, conservative otherwise)
    d, j = tree.query(np.column_stack([py, px]), k=1)
    k = kernel(d)
    # exact metric: TP = sum_g max_x p(x)k(d); compute via scatter-max
    mx = np.zeros(n, float)
    np.maximum.at(mx, j, p[py, px] * k)
    tp = float(mx.sum())
    fn = float(n) - tp
    fp = float(((1.0 - kernel(d)) * p[py, px]).sum())
    return float(tp / (tp + ALPHA * fp + BETA * fn + EPS))


def main():
    with rasterio.open(DATA / "labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    cat = lab == 1
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as s:
        sg = s.read(1) > 0
    sg_off = sg & foot & ~binary_dilation(cat, iterations=3)

    print("building catalogue-hidden holdout (2x2 blocks x 3 seeds, hide 25% of components)...")
    ho = inst.build_holdout(foot, cat, hide_frac=0.25, n_y=2, n_x=2,
                            seeds=(11, 12, 13), collar=2)
    print(f"  {len(ho.cells)} cells; hidden truth px per cell: "
          f"{[c.n_truth for c in ho.cells]}")

    rows = []
    for name, score_s, url in csv.reader((ROOT / "registry" / "score_ledger.csv").open()):
        p = DATA / "scored" / name
        if not p.exists():
            continue
        with rasterio.open(p) as s:
            if (s.height, s.width) != (3730, 3292) or \
               tuple(round(v, 3) for v in s.transform)[:6] != REF_TR:
                continue
            a = s.read(1)
        a = np.where(foot, np.clip(np.nan_to_num(a, nan=0.0), 0.0, 1.0), 0.0)
        ev = inst.evaluate(a, ho)
        cells = ev["cells"]
        cv = np.array([v["dti"] for v in cells.values()])
        rows.append(dict(
            name=name,
            score=None if score_s in ("NA", "") else float(score_s),
            n_pos=int((a > 0).sum()),
            H_pooled=ev["pooled_dti"],
            H_mean=float(cv.mean()),
            H_min=float(cv.min()),
            H_cells={k: v["dti"] for k, v in cells.items()},
            C_cat=dti(a > 0.5, cat, foot),
            C_cat_masked=dti(a > 0.5, cat, foot, known=cat),
            S_sgmc=dti(a > 0.5, sg_off, foot, known=cat),
        ))
        r = rows[-1]
        print(f"  {name[:50]:52s} n={r['n_pos']:7d} H={r['H_pooled']:.5f} "
              f"C={r['C_cat']:.4f} S={r['S_sgmc']:.4f} score={r['score']}")

    sc = [r for r in rows if r["score"] is not None]
    y = np.array([r["score"] for r in sc])
    logn = np.log(np.array([r["n_pos"] for r in sc], float))
    res = {}
    for k in ("H_pooled", "H_mean", "H_min", "C_cat", "C_cat_masked", "S_sgmc", "n_pos"):
        x = np.array([r[k] for r in sc], float)
        rho, pv = spearmanr(x, y)
        prho, ppv = partial_spearman(x, y, logn)
        res[k] = dict(spearman=float(rho), p=float(pv),
                      partial_given_logn=float(prho), partial_p=float(ppv))
    print(f"\n--- instrument calibration (n={len(sc)} artifacts with owner-reported scores) ---")
    print(f"{'instrument':16s} {'rho':>8s} {'p':>9s} {'rho|log n':>10s} {'p':>9s}")
    for k, v in res.items():
        print(f"{k:16s} {v['spearman']:+8.3f} {v['p']:9.4f} "
              f"{v['partial_given_logn']:+10.3f} {v['partial_p']:9.4f}")

    (ROOT / "registry" / "instrument_calibration.json").write_text(
        json.dumps(dict(n_artifacts=len(sc), n_cells=len(ho.cells),
                        hide_frac=0.25, blocks="2x2", seeds=[11, 12, 13],
                        collar=2, instruments=res,
                        artifacts=[{kk: vv for kk, vv in r.items() if kk != "H_cells"}
                                   for r in rows]), indent=2, default=float) + "\n")
    print("\nwrote registry/instrument_calibration.json")


if __name__ == "__main__":
    main()
