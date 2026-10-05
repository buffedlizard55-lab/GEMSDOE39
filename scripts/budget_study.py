#!/usr/bin/env python3
"""Choose the emission budget from the DTI break-even rule, not from a surrogate.

Derivation (from the organizer's own equations, page 967)
--------------------------------------------------------
DTI = TP / (TP + 0.2*FP + 0.8*FN + eps),  FN = |G| - TP,  and for a binary
emission every emitted pixel x contributes kernel weight k(d(x,G)) to TP (via the
nearest truth pixel) and 1 - k(d(x,G)) to FP.  Adding one pixel with expected
kernel credit w therefore changes the denominator by

    dTP + 0.2*dFP + 0.8*dFN = w + 0.2(1-w) + 0.8(-w) = 0.2

so the added pixel pays for itself iff

    w > 0.2 * DTI                                       (break-even)

independent of everything else.  At the historical best DTI of 0.2778 that is
w* = 0.0556; GEMSDOE32 measured the same bar empirically (0.0548) and derived it
as 0.2 x 0.26 = 0.0520.  Two independent routes, one number.

What we need, then, is the marginal expected credit w(n) of the n-th ranked dot.
We measure its shape on the catalogue-hidden holdout and rescale it to the
prevalence of the real (unobserved) truth set, which is declared in advance.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gems39 import instrument as inst, emission as em  # noqa: E402
from gems39.metric import kernel  # noqa: E402

DATA = ROOT / "data"


def load():
    with rasterio.open(DATA / "labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    cat = lab == 1
    return foot, cat


def marginal_credit(field, ho, n_grid, cat_for_buffer, buffer_px=2, spacing=3.0):
    """For each n in n_grid, emit n dots on the holdout domain and return
    (n, TP_w summed over cells, per-dot mean credit)."""
    from scipy.ndimage import distance_transform_edt as edt
    # pooled emission over the whole footprint, then evaluate per cell
    out = []
    prev = None
    for n in n_grid:
        m = em.emit_fast(field, ho.foot, target_n=n, min_dist=spacing,
                         catalogue=cat_for_buffer, cat_buffer_px=buffer_px)
        TP = FP = FN = 0.0
        for c in ho.cells:
            p = m[c.sl] & c.active
            if not p.any():
                FN += c.n_truth
                continue
            dp = edt(~p)
            d = dp[c.ty, c.tx]
            k = kernel(d)
            tp = float(k.sum())
            py, px = np.nonzero(p)
            fp = float((1.0 - c.k_dg[py, px]).sum())
            TP += tp
            FP += fp
            FN += float(c.n_truth) - tp
        dots = int(m.sum())
        out.append(dict(target=n, dots=dots, tp=TP, fp=FP, fn=FN,
                        dti=float(TP / (TP + 0.2 * FP + 0.8 * FN + 1e-12)),
                        mean_credit=float(TP / max(dots, 1))))
        del m
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default="A", choices=["A", "B", "C"])
    ap.add_argument("--k-real", type=float, default=30000.0,
                    help="declared prior size (px) of the scored truth set")
    ap.add_argument("--spacing", type=float, default=3.0)
    ap.add_argument("--buffer", type=int, default=2)
    ap.add_argument("--dti", type=float, default=0.27)
    args = ap.parse_args()

    foot, cat = load()
    z = np.load(DATA / "probs.npz")
    PA, PB = z["PA"].astype(np.float32), z["PB"].astype(np.float32)
    pA, pB = np.clip(PA, 0, 1), np.clip(PB, 0, 1)
    field = {"A": pA, "B": pB, "C": np.sqrt(pA * pB)}[args.field]
    field = np.where(foot, field, 0.0).astype(np.float32)

    ho = inst.build_holdout(foot, cat, hide_frac=0.25, n_y=2, n_x=2,
                            seeds=(21, 22, 23), collar=2)
    k_ho = float(sum(c.n_truth for c in ho.cells))
    print(f"holdout hidden-truth pixels (summed over cells): {k_ho:,.0f}")
    print(f"declared real truth size K_real: {args.k_real:,.0f}")
    scale = args.k_real / max(k_ho, 1.0)
    print(f"prevalence scale w_real = w_holdout x {scale:.3f}")
    bar = 0.2 * args.dti
    print(f"break-even credit bar at DTI={args.dti}: w* = {bar:.4f}\n")

    grid = [8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000,
            44000, 50000, 56000, 64000]
    res = marginal_credit(field, ho, grid, cat, args.buffer, args.spacing)
    print(f"{'n':>7s} {'dots':>7s} {'TP_w':>9s} {'FP_w':>9s} {'hoDTI':>8s} "
          f"{'mean_w':>7s} {'marg_w':>7s} {'marg_scaled':>11s} {'>bar':>5s}")
    prev = None
    rec = []
    for r in res:
        if prev is None:
            marg = r["tp"] / max(r["dots"], 1)
        else:
            marg = (r["tp"] - prev["tp"]) / max(r["dots"] - prev["dots"], 1)
        ms = marg * scale
        rec.append(dict(**r, marg_w=marg, marg_scaled=ms))
        print(f"{r['target']:7d} {r['dots']:7d} {r['tp']:9.1f} {r['fp']:9.1f} "
              f"{r['dti']:8.5f} {r['mean_credit']:7.4f} {marg:7.4f} {ms:11.4f} "
              f"{'YES' if ms > bar else 'no':>5s}")
        prev = r
    # largest n whose marginal scaled credit still clears the bar
    ok = [r for r in rec if r["marg_scaled"] > bar]
    pick = max(ok, key=lambda r: r["dots"])["target"] if ok else min(rec, key=lambda r: r["dots"])["target"]
    print(f"\nbreak-even budget (field {args.field}, K_real={args.k_real:.0f}, "
          f"spacing={args.spacing}, buffer={args.buffer}): n = {pick}")
    (ROOT / "registry" / f"budget_study_{args.field}.json").write_text(
        json.dumps(dict(args=vars(args), k_holdout=k_ho, scale=scale, bar=bar,
                        pick=pick, table=rec), indent=2, default=float) + "\n")


if __name__ == "__main__":
    main()
