#!/usr/bin/env python3
"""Forensic analysis of historical scored artifacts.

Goal: establish, from bytes we can re-read, WHY the highest-scoring family wins,
and calibrate a local proxy instrument against owner-reported live scores.

Outputs (written to registry/):
  forensics_artifacts.csv   per-artifact geometry / emission statistics
  forensics_instruments.json  Spearman correlation of each candidate proxy with
                              the owner-reported live score

Evidence classes are preserved: every "score" here is OWNER-REPORTED, not
organizer-authenticated (no public receipt ties a hash to a score).
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt, binary_dilation
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39.metric import ALPHA, BETA, EPS, kernel  # noqa: E402

DATA = ROOT / "data"
REF_TR = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def load_grid():
    with rasterio.open(DATA / "labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    cat = lab == 1
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as s:
        sg = s.read(1) > 0
    sg_off = sg & foot & ~cat & ~binary_dilation(cat, iterations=3)
    return foot, cat, sg_off


def dti_of(pred_bool, truth, valid, known=None):
    """Distance-weighted Tversky with an optional masked (known) set."""
    p = np.asarray(pred_bool, bool) & valid
    if known is not None:
        p = p & ~known
    g = truth & valid
    if known is not None:
        g = g & ~known
    n = int(g.sum())
    if n == 0 or not p.any():
        return dict(tp=0.0, fp=float(p.sum()), fn=float(n), n_truth=n, dti=0.0)
    dp = distance_transform_edt(~p)
    tp = float(kernel(dp[g]).sum())
    fn = float(n) - tp
    dg = distance_transform_edt(~g)
    fp = float((1.0 - kernel(dg[p])).sum())
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n,
                dti=float(tp / (tp + ALPHA * fp + BETA * fn + EPS)))


def main():
    foot, cat, sg_off = load_grid()
    cat_near2 = binary_dilation(cat, iterations=2)
    dcat = distance_transform_edt(~cat)
    print(f"footprint={int(foot.sum()):,}  catalogue={int(cat.sum()):,}  "
          f"sgmc_off_catalogue(d>3px)={int(sg_off.sum()):,}")

    rows = []
    ledger = list(csv.reader((ROOT / "registry" / "score_ledger.csv").open()))
    for name, score_s, url in ledger:
        p = DATA / "scored" / name
        if not p.exists():
            print(f"  MISSING {name}")
            continue
        with rasterio.open(p) as s:
            if (s.height, s.width) != (3730, 3292) or \
               tuple(round(v, 3) for v in s.transform)[:6] != REF_TR:
                print(f"  GRID-MISMATCH (skipped) {name}")
                continue
            a = s.read(1)
        a = np.where(foot, a, 0.0)
        a = np.where(np.isfinite(a), a, 0.0)
        a = np.clip(a, 0.0, 1.0)
        pos = a > 0
        n = int(pos.sum())
        mass = float(a.sum())
        # nearest-neighbour spacing between positive pixels
        if 0 < n <= 400000:
            yy, xx = np.nonzero(pos)
            from scipy.spatial import cKDTree
            samp = min(n, 60000)
            idx = np.random.default_rng(0).choice(n, samp, replace=False)
            tree = cKDTree(np.column_stack([yy[idx], xx[idx]]))
            dd, _ = tree.query(np.column_stack([yy[idx], xx[idx]]), k=2)
            nn_med = float(np.median(dd[:, 1]))
            nn_mean = float(np.mean(dd[:, 1]))
        else:
            nn_med = nn_mean = float("nan")
        rows.append(dict(
            name=name,
            score=None if score_s in ("NA", "") else float(score_s),
            url=url,
            n_pos=n,
            mass=mass,
            binaryness=float(mass / max(n, 1)),
            nn_median=nn_med,
            nn_mean=nn_mean,
            frac_on_cat=float((pos & cat).sum()) / max(n, 1),
            frac_within_2px_cat=float((pos & cat_near2).sum()) / max(n, 1),
            frac_within_3px_cat=float((pos & (dcat <= 3)).sum()) / max(n, 1),
            median_dist_cat=float(np.median(dcat[pos])) if n else float("nan"),
            dti_catalogue=dti_of(a > 0.5, cat, foot)["dti"],
            dti_cat_masked=dti_of(a > 0.5, cat, foot, known=cat)["dti"],
            dti_sgmc_off=dti_of(a > 0.5, sg_off, foot, known=cat)["dti"],
            dti_sgmc_off_soft=dti_of(a, sg_off, foot, known=cat)["dti"],
        ))
        r = rows[-1]
        print(f"  {name[:52]:54s} n={n:7d} nn={nn_med:5.2f} oncat={r['frac_on_cat']:.3f} "
              f"w2={r['frac_within_2px_cat']:.3f} dti_cat={r['dti_catalogue']:.4f} "
              f"dti_sgmc={r['dti_sgmc_off']:.4f} score={r['score']}")

    # --- instrument calibration -------------------------------------------
    scored = [r for r in rows if r["score"] is not None]
    keys = ["n_pos", "nn_median", "nn_mean", "frac_on_cat", "frac_within_2px_cat",
            "frac_within_3px_cat", "median_dist_cat", "dti_catalogue",
            "dti_cat_masked", "dti_sgmc_off", "dti_sgmc_off_soft"]
    out = {}
    y = np.array([r["score"] for r in scored])
    for k in keys:
        x = np.array([r[k] for r in scored], float)
        ok = np.isfinite(x)
        rho, pv = spearmanr(x[ok], y[ok])
        out[k] = dict(spearman=float(rho), p=float(pv), n=int(ok.sum()))
    print("\n--- instrument vs owner-reported live score (n=%d artifacts) ---" % len(scored))
    for k, v in sorted(out.items(), key=lambda kv: -abs(kv[1]["spearman"])):
        print(f"  {k:24s} rho={v['spearman']:+.3f}  p={v['p']:.4f}")

    with (ROOT / "registry" / "forensics_artifacts.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    (ROOT / "registry" / "forensics_instruments.json").write_text(
        json.dumps(dict(n_scored=len(scored), instruments=out,
                        artifacts=[{k: v for k, v in r.items() if k != "url"}
                                   for r in rows]), indent=2, default=float) + "\n")
    print("\nwrote registry/forensics_artifacts.csv and registry/forensics_instruments.json")


if __name__ == "__main__":
    main()
