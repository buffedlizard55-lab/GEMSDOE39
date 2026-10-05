#!/usr/bin/env python3
"""Wald SPRT for paired fold-level holdout outcomes.
H1: candidate wins with probability p1; H0: win probability p0.
One observation per independent spatial block; do not use per-pixel rows as independent.
"""
import argparse, math

def sprt(outcomes, p0=.5, p1=.7, alpha=.05, beta=.10):
    if not (0 < p0 < p1 < 1 and 0 < alpha < 1 and 0 < beta < 1): raise ValueError('require 0<p0<p1<1 and valid errors')
    upper=math.log((1-beta)/alpha); lower=math.log(beta/(1-alpha)); llr=0.0
    for i, win in enumerate(outcomes, 1):
        llr += math.log((p1 if win else 1-p1)/(p0 if win else 1-p0))
        if llr >= upper: return {'decision':'accept_H1','n':i,'llr':llr,'upper':upper,'lower':lower}
        if llr <= lower: return {'decision':'accept_H0','n':i,'llr':llr,'upper':upper,'lower':lower}
    return {'decision':'continue','n':len(outcomes),'llr':llr,'upper':upper,'lower':lower}

if __name__ == '__main__':
 p=argparse.ArgumentParser(); p.add_argument('outcomes', nargs='+', choices=['1','0']); p.add_argument('--p0',type=float,default=.5); p.add_argument('--p1',type=float,default=.7); p.add_argument('--alpha',type=float,default=.05); p.add_argument('--beta',type=float,default=.1); a=p.parse_args(); print(sprt([x=='1' for x in a.outcomes],a.p0,a.p1,a.alpha,a.beta))
