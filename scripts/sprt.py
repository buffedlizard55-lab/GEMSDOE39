#!/usr/bin/env python3
"""Wald Sequential Probability Ratio Test CLI.

Usage:
  python scripts/sprt.py 1 1 0 1 1 1 0 1
  (each positional arg is 1 = candidate wins that fold, 0 = baseline wins)

Defaults (pre-declared for GEMSDOE39):
  p0 = 0.5   (H0: candidate doesn't beat baseline — win prob 0.5 = chance)
  p1 = 0.7   (H1: candidate beats baseline with 70% win prob per fold)
  alpha = 0.05  (Type-I error: falsely declaring an improvement)
  beta  = 0.10  (Type-II error: missing a real improvement)
"""
import argparse, math, sys, json
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'src'))
from gems39.sprt_select import sprt_pairwise  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('outcomes', nargs='*', help='Per-fold outcomes: 1 = candidate wins, 0 = baseline wins.')
    p.add_argument('--p0', type=float, default=0.5)
    p.add_argument('--p1', type=float, default=0.7)
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--beta', type=float, default=0.10)
    args = p.parse_args()
    wins = [bool(int(x)) for x in args.outcomes] if args.outcomes else []
    r = sprt_pairwise(wins, p0=args.p0, p1=args.p1, alpha=args.alpha, beta=args.beta)
    print(json.dumps(r, indent=2))
    if r['decision'] == 'accept_H1':
        print('\n>>> STOP: accept H1 (candidate beats baseline at declared error rates).')
    elif r['decision'] == 'accept_H0':
        print('\n>>> STOP: accept H0 (candidate does NOT beat baseline at declared error rates).')
    else:
        print(f'\n>>> CONTINUE sampling (n={r["n"]}, llr={r["llr"]:.3f} in [{r["lower"]:.3f}, {r["upper"]:.3f}]).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
