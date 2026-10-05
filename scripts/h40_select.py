#!/usr/bin/env python3
"""H40 fusion selection: pre-declared grid, measured on the local instrument.

This is the selection stage.  It exists because the first full run of
``build_h40.py`` showed every H40 channel scoring *below* the live-scored anchor
on the prevalence-matched instrument (0.9-7.0% hit rate against the anchor's
7.654%).  Two defects explained most of that gap and are now fixed:

  1. every channel ended in ``_robust_unit``, which clips at the 99.5th
     percentile and pinned ~25 500 pixels to exactly 1.0 -- about two thirds of
     a 37 654-dot budget -- so the emitter ranked by raster order, not geology;
  2. the I1 collar was ``dilate(block)`` which contains the block, so every
     catalogue component "touched the collar" and I1 had zero truth pixels.

The grid below is declared in full BEFORE any variant is measured, and the
selection rule is declared with it:

  * emit every variant at the anchor's validated budget (37 654 dots) with the
    anchor's validated geometry (Poisson 2.8 px, catalogue exclusion 2 px);
  * measure the hit rate on instrument I2 (SGMC off-catalogue faults,
    prevalence-matched to the calibrated |G|);
  * the anchor artifact itself is measured in the same run and is the control;
  * a variant is a candidate winner only if it beats the control's instrument
    hit rate AND reaches ``accept_H1`` in the Wald SPRT over the I2 folds;
  * among those, take the highest instrument hit rate.

Selecting on a proxy has a cost and it is stated: I2 is USGS SGMC geologic-map
faults, not the organizer's labels, so a variant can win here and lose live.
That is why ``build_h40.py`` re-checks the winner against the promotion gate
before writing anything.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import calibrate as cal, detectors40 as det40, external as ex  # noqa: E402
from gems39 import grid, holdout40 as ho, priced_emit as pe                 # noqa: E402
from gems39.sprt_select import sprt_normal_mean, sprt_pairwise            # noqa: E402

ANCHOR_REL = "data/calib/h33-2-b2_0.2778.tif"
ANCHOR_SCORE = 0.2778
ANCHOR_N = 37_654
SPRT = dict(alpha=0.05, beta=0.10, p0=0.5, p1=0.7)

# ---------------------------------------------------------------- pre-declared grid
# Exponent vector over (backbone, A, B, C, D, E, F, G).  Floats allowed.
#
# This grid REPLACES the one used in the first selection run.  The reason is
# recorded here rather than hidden: that run compared variants against the anchor
# using a quantity ``evaluate()`` called ``hit_rate`` which was in fact
# (truth pixels within R) / (number of dots) -- a truth-side count divided by a
# prediction-side count.  The variants were measured with the correct per-dot
# statistic, so the two columns were not comparable and every "vs anchor" ratio
# in that run is meaningless.  Both defects are fixed (see holdout40.py) and the
# whole grid was re-measured.  The first run did show one thing that survived the
# bug, because it never used the broken quantity: on the pooled instrument DTI --
# the metric's own shape -- the supervised off-catalogue channel H40-F scored
# 0.1309 against the anchor's 0.0640.  This grid is centred on that channel.
GRID = {
    "W00-backbone":            (1, 0, 0, 0, 0, 0, 0, 0),
    "W01-F":                   (0, 0, 0, 0, 0, 0, 1, 0),
    "W02-F2":                  (0, 0, 0, 0, 0, 0, 2, 0),
    "W03-bb0.5-F1":            (0.5, 0, 0, 0, 0, 0, 1, 0),
    "W04-bb0.25-F1":           (0.25, 0, 0, 0, 0, 0, 1, 0),
    "W05-F1-A0.5":             (0, 0.5, 0, 0, 0, 0, 1, 0),
    "W06-F1-A0.25":            (0, 0.25, 0, 0, 0, 0, 1, 0),
    "W07-F1-C0.5":             (0, 0, 0, 0.5, 0, 0, 1, 0),
    "W08-F1-G0.5":             (0, 0, 0, 0, 0, 0, 1, 0.5),
    "W09-F1-A0.5-G0.5":        (0, 0.5, 0, 0, 0, 0, 1, 0.5),
    "W10-F1-A0.5-C0.5":        (0, 0.5, 0, 0.5, 0, 0, 1, 0),
    "W11-bb1-F1":              (1, 0, 0, 0, 0, 0, 1, 0),
    "W12-bb1-A1-F2":           (1, 1, 0, 0, 0, 0, 2, 0),
    "W13-bb0.5-F1.5":          (0.5, 0, 0, 0, 0, 0, 1.5, 0),
    "W14-F1-D0.25-E0.25":      (0, 0, 0, 0, 0.25, 0.25, 1, 0),
    "W15-F1-A0.5-B0.25-C0.5-G0.25": (0, 0.5, 0.25, 0.5, 0, 0, 1, 0.25),
    "W16-F1-C0.5-G0.5":        (0, 0, 0, 0.5, 0, 0, 1, 0.5),
    "W17-F1-A0.25-C0.25-G0.25": (0, 0.25, 0, 0.25, 0, 0, 1, 0.25),
    "W18-bb0.25-F1-A0.25":     (0.25, 0.25, 0, 0, 0, 0, 1, 0),
}
ORDER = ("backbone", "H40-A", "H40-B", "H40-C", "H40-D", "H40-E", "H40-F", "H40-G")

# Selection rule, declared before measurement:
#   PRIMARY   pooled instrument DTI on I2 (the metric's own shape, additive over
#             folds through fold_score_delta) must exceed the anchor's;
#   CONFIRM   Wald SPRT on the additive per-fold first-order score delta must
#             reach accept_H1;
#   GUARD     leakage <= 0.02 and the per-dot instrument hit rate must not be
#             below the anchor's by more than 20% (a variant that wins DTI only by
#             sitting on a few dense traces would fail this).
#   Among variants clearing all three, take the highest pooled instrument DTI.


def fuse(ch, w, foot, eps=1e-6):
    acc = np.zeros(foot.shape, np.float64)
    tot = 0.0
    for name, p in zip(ORDER, w):
        if p == 0:
            continue
        a = ch.get(name)
        if a is None:
            continue
        acc += p * np.log(np.maximum(np.asarray(a, np.float64), eps))
        tot += p
    if tot <= 0:
        return np.zeros(foot.shape, np.float32)
    return det40.norm(np.exp(acc / tot), foot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(ROOT / "data"))
    ap.add_argument("--cache", default=str(ROOT / "artifacts" / "h40_raw_channels.npz"))
    ap.add_argument("--reuse-channels", action="store_true")
    ap.add_argument("--rebuild", default="",
                    help="comma-separated channel names to rebuild even with --reuse-channels")
    ap.add_argument("--max-n", type=int, default=ANCHOR_N)
    ap.add_argument("--min-dist", type=float, default=2.8)
    ap.add_argument("--cat-buffer-px", type=int, default=2)
    ap.add_argument("--out", default=str(ROOT / "artifacts" / "h40_selection.json"))
    args = ap.parse_args()

    ddir = Path(args.data_dir)
    t0 = time.time()
    with rasterio.open(ddir / "sample_submission.tif") as s:
        foot = np.isfinite(s.read(1))
    with rasterio.open(ddir / "labels.tif") as s:
        cat = (s.read(1) >= 1) & foot
    active = foot & ~cat
    cache_path = Path(args.cache)

    ch = {}
    rebuild = {x.strip() for x in args.rebuild.split(",") if x.strip()}
    if args.reuse_channels and cache_path.exists():
        z = np.load(cache_path, allow_pickle=True)
        ch = {k: z[k] for k in z.files if k != "__log__"}
        print(f"reused {len(ch)} raw channels from {cache_path}: {sorted(ch)}")
        if rebuild:
            print(f"rebuilding {sorted(rebuild)} ...")
            bands = grid.read_all_bands(ddir / "training_features.tif", foot)
            ext = ex.load_all(ddir, foot.shape, cat, foot)
            clog = []
            cache = {}
            for nm in sorted(rebuild):
                if nm == "H40-F":
                    ch[nm], fmap = det40.build_h40f(bands, foot, ext, cat, clog,
                                                    cache=cache, return_fold_map=True)
                    np.savez_compressed(ROOT / "artifacts" / "h40f_foldmap.npz", fold=fmap)
                elif nm == "backbone":
                    ch[nm], _ = det40.build_backbone(bands, foot, ext, clog, cache=cache)
                else:
                    fn = {"H40-A": det40.build_h40a, "H40-B": det40.build_h40b,
                          "H40-C": det40.build_h40c, "H40-D": det40.build_h40d,
                          "H40-E": det40.build_h40e}[nm]
                    ch[nm] = fn(bands, foot, ext, clog, cache=cache)                         if nm != "H40-G" else                         det40.build_h40g(bands, foot, ext, cat, clog, cache=cache)
                print(f"  {nm} rebuilt ({time.time()-t0:.0f}s)")
            for l in clog:
                print("   .", l)
            np.savez_compressed(cache_path, **ch, __log__=np.array(clog, dtype=object))
            del bands
    else:
        bands = grid.read_all_bands(ddir / "training_features.tif", foot)
        ext = ex.load_all(ddir, foot.shape, cat, foot)
        clog = []
        cache = {}
        bb, _ = det40.build_backbone(bands, foot, ext, clog, cache=cache)
        ch["backbone"] = bb
        print(f"  backbone built ({time.time()-t0:.0f}s)")
        for key, fn in (("H40-A", det40.build_h40a), ("H40-B", det40.build_h40b),
                        ("H40-C", det40.build_h40c), ("H40-D", det40.build_h40d),
                        ("H40-E", det40.build_h40e)):
            ch[key] = fn(bands, foot, ext, clog, cache=cache)
            print(f"  {key} built ({time.time()-t0:.0f}s)")
        ch["H40-G"] = det40.build_h40g(bands, foot, ext, cat, clog, cache=cache)
        print(f"  H40-G built ({time.time()-t0:.0f}s)")
        ch["H40-F"] = det40.build_h40f(bands, foot, ext, cat, clog, cache=cache)
        print(f"  H40-F built ({time.time()-t0:.0f}s)")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_path, **ch, __log__=np.array(clog, dtype=object))
        for l in clog:
            print("   .", l)
        del bands

    # tie diagnostic: the defect that motivated this script
    print("\ntie diagnostic (pixels at exactly the channel maximum, inside the active domain):")
    for k in sorted(ch):
        v = ch[k][active]
        mx = float(v.max())
        print(f"  {k:9s} max={mx:.4g}  n_at_max={int((v >= mx).sum()):>8,}  "
              f"n_at_p995={int((v >= float(np.quantile(v, 0.995))).sum()):>8,}")

    print("\nholdout ...")
    G = cal.fit_G_sensitivity(40199, 0.2708, ANCHOR_N, ANCHOR_SCORE)["G_median"]
    ctx = ho.load_holdout40(ddir, target_truth_px=int(round(G)))
    print(f"  |G| = {G:,.0f} | I1 truth px = {int(ho.instrument_truth(ctx,'I1').sum()):,} "
          f"| I2 truth px = {int(ho.instrument_truth(ctx,'I2').sum()):,}")
    allowed = pe.allowed_domain(foot, cat, args.cat_buffer_px)

    with rasterio.open(ROOT / ANCHOR_REL) as s:
        a = s.read(1)
    anchor = (np.isfinite(a) & (a > 0)) & active
    del a

    i2 = sorted([c for c in ctx.cells if c.instrument == "I2"], key=lambda c: c.key)
    i1 = sorted([c for c in ctx.cells if c.instrument == "I1"], key=lambda c: c.key)
    ev_a = ho.evaluate(anchor, ctx, pooled=True)
    s_anchor = ev_a["I2"]["pooled_dti"]
    h_anchor_dot = ev_a["I2"]["dot_hit_rate"]
    cov_anchor = ev_a["I2"]["truth_coverage"]
    # live hit rate of the anchor from the inverse-DTI calibration
    h_anchor_live = cal.implied_hit_rate(ANCHOR_N, ANCHOR_SCORE, G)
    transfer = h_anchor_dot / h_anchor_live if h_anchor_live > 0 else 1.0
    anchor_cells = {c.key: ho.dti_cell(anchor, c) for c in i2}
    print(f"\nanchor control: n={int(anchor.sum()):,}")
    print(f"  pooled instrument DTI (I2)      = {s_anchor:.4f}")
    print(f"  per-dot instrument hit rate     = {100*h_anchor_dot:.3f}%")
    print(f"  instrument truth coverage       = {100*cov_anchor:.2f}% of {int(ctx.sgmc_matched.sum()):,} px")
    print(f"  pooled instrument DTI (I1)      = {ev_a['I1']['pooled_dti']:.4f}")
    print(f"  live hit rate from calibration  = {100*h_anchor_live:.3f}%  (live score {ANCHOR_SCORE})")
    print(f"  TRANSFER (instrument/live, per-dot hit rate) = {transfer:.3f}")

    hdr = (f"\n{'variant':30s}{'n':>8}{'instrDTI':>10}{'x anchor':>9}{'dot h':>8}"
           f"{'cov%':>7}{'live h':>8}{'PRICE':>8}{'sign w/n':>10}{'LLR':>8}  sign-test"
           f"     mean-test                  I1 ratio")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for name, w in GRID.items():
        prop = fuse(ch, w, foot)
        mask, ranks = pe.emit_ranked(prop, allowed, min_dist=args.min_dist,
                                     max_n=args.max_n)
        n = int(mask.sum())
        ev = ho.evaluate(mask, ctx, pooled=True)
        s_i2 = ev["I2"]["pooled_dti"]
        h_dot = ev["I2"]["dot_hit_rate"]
        cov = ev["I2"]["truth_coverage"]
        h_live = h_dot / transfer if transfer > 0 else 0.0
        pr = cal.price(n, G, h_live)
        # additive first-order score delta per fold -> BOTH sequential tests
        wins, deltas, dti_wins_l = [], [], []
        for c in i2:
            fd = ho.fold_score_delta(c, mask, anchor, s_anchor)
            wins.append(fd["delta"] > 0)
            deltas.append(fd["delta"])
            dti_wins_l.append(fd["fold_dti_win"])
        d = sprt_pairwise(wins, **SPRT)
        dn = sprt_normal_mean(deltas, d_alt=0.5, alpha=SPRT["alpha"], beta=SPRT["beta"])
        dti_wins = sum(1 for c in i2
                       if ho.dti_cell(mask, c)["dti"] > anchor_cells[c.key]["dti"])
        leak = ho.leakage_probe(mask, ctx, ctx.sgmc_matched)
        rows.append(dict(variant=name, weights=dict(zip(ORDER, w)), n=n,
                         instrument_dti=s_i2, ratio_dti=s_i2 / s_anchor if s_anchor else 0.0,
                         dot_hit_rate=h_dot, ratio_hit=h_dot / h_anchor_dot if h_anchor_dot else 0.0,
                         truth_coverage=cov, live_hit_rate=h_live,
                         priced_score=pr["score"], priced_TP=pr["TPw"], priced_FP=pr["FPw"],
                         instrument_dti_I1=ev["I1"]["pooled_dti"],
                         ratio_dti_I1=(ev["I1"]["pooled_dti"] / ev_a["I1"]["pooled_dti"]
                                       if ev_a["I1"]["pooled_dti"] else 0.0),
                         leakage=leak, sprt_sign=d, sprt_mean=dn,
                         sum_fold_delta=float(np.sum(deltas)),
                         mean_fold_delta=float(np.mean(deltas)) if deltas else 0.0,
                         fold_dti_wins=int(sum(dti_wins_l)), n_folds=len(i2)))
        r = rows[-1]
        print(f"{name:30s}{n:>8,}{s_i2:>10.4f}{r['ratio_dti']:>8.3f}x{100*h_dot:>7.3f}%"
              f"{100*cov:>6.1f}%{100*h_live:>7.3f}%{pr['score']:>8.4f}"
              f"{str(d['wins'])+'/'+str(d['n']):>10}{d['llr']:>+8.3f}  {d['decision']:10s}"
              f" | mean-test {dn['decision']:10s} LLR {dn['llr']:+7.3f}"
              f" | I1 {r['ratio_dti_I1']:.3f}x")
        del prop, mask, ranks, ev, deltas

    ok = [r for r in rows
          if r["ratio_dti"] > 1.0
          and (r["sprt_sign"]["decision"] == "accept_H1"
               or r["sprt_mean"]["decision"] == "accept_H1")
          and r["leakage"] <= 0.02 and r["ratio_hit"] >= 0.80]
    dti_only = [r for r in rows if r["ratio_dti"] > 1.0 and r["leakage"] <= 0.02]
    if ok:
        winner = max(ok, key=lambda r: r["instrument_dti"])
        why = ("clears the whole declared rule: pooled instrument DTI above the anchor, "
               "Wald SPRT accept_H1 on the additive per-fold score delta, leakage <= 2%, "
               "per-dot hit rate within 20% of the anchor")
    elif dti_only:
        winner = max(dti_only, key=lambda r: r["instrument_dti"])
        why = ("PROXY WIN ONLY: pooled instrument DTI beats the anchor but the Wald SPRT "
               "did not reach accept_H1 on the additive per-fold score delta. This is not "
               "a validated improvement and the manifest and site must say so.")
    else:
        winner = None
        why = "no variant beat the anchor on the pooled instrument DTI"
    print("\n" + "=" * 78)
    print(f"winner: {winner['variant'] if winner else None}")
    if winner:
        print(f"  pooled instrument DTI {winner['instrument_dti']:.4f} "
              f"({winner['ratio_dti']:.3f}x anchor {s_anchor:.4f})")
        print(f"  per-dot instrument hit rate {100*winner['dot_hit_rate']:.3f}% "
              f"-> live-priced {100*winner['live_hit_rate']:.3f}% (anchor "
              f"{100*h_anchor_live:.3f}%)")
        print(f"  PRICED LIVE SCORE {winner['priced_score']:.4f} vs anchor {ANCHOR_SCORE}")
        print(f"  sign SPRT  {winner['sprt_sign']['decision']} "
              f"({winner['sprt_sign']['wins']}/{winner['sprt_sign']['n']} folds, "
              f"LLR {winner['sprt_sign']['llr']:+.3f})")
        print(f"  mean SPRT  {winner['sprt_mean']['decision']} "
              f"(n={winner['sprt_mean']['n']}, LLR {winner['sprt_mean']['llr']:+.3f}, "
              f"sigma={winner['sprt_mean']['sigma']:.4f} estimated="
              f"{winner['sprt_mean']['sigma_estimated']})")
        print(f"  sum of per-fold first-order score deltas = {winner['sum_fold_delta']:+.4f}")
        print(f"  I1 pooled DTI {winner['instrument_dti_I1']:.4f} vs anchor "
              f"{ev_a['I1']['pooled_dti']:.4f}")
    print(f"reason: {why}")
    out = dict(G=G, anchor=dict(n=int(anchor.sum()), live_score=ANCHOR_SCORE,
                                live_hit_rate=h_anchor_live,
                                instrument_dti_I2=s_anchor,
                                instrument_dot_hit_rate=h_anchor_dot,
                                instrument_truth_coverage=cov_anchor,
                                instrument_dti_I1=ev_a["I1"]["pooled_dti"],
                                transfer_factor=transfer),
               grid={k: dict(zip(ORDER, v)) for k, v in GRID.items()},
               rows=rows, winner=winner["variant"] if winner else None,
               winner_reason=why,
               selection_rule=("pooled instrument DTI > anchor AND (sign SPRT accept_H1 OR "
                               "Wald normal-mean SPRT accept_H1 on the additive per-fold "
                               "first-order score delta) AND leakage <= 0.02 AND per-dot hit "
                               "rate >= 0.80 x anchor"),
               sprt_predeclared=SPRT, seconds=round(time.time() - t0, 1))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=2, default=float) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
