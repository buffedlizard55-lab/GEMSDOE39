"""Wald's Sequential Probability Ratio Test for paired fold-level outcomes.

Pre-declared parameters (alpha, beta, p0, p1) are set before peeking at any fold
score.  After each independent spatial-fold comparison we accumulate the
log-likelihood ratio and compare against the Wald boundaries:
    A = log(beta / (1 - alpha))   (accept H0: candidate doesn't beat baseline)
    B = log((1 - beta) / alpha)   (accept H1: candidate beats baseline)
Stopping is principled: no visual peeking, no "run a bit longer" bias.
"""
import math


def sprt_pairwise(wins, p0=0.5, p1=0.7, alpha=0.05, beta=0.10):
    """Return dict with decision ('accept_H0', 'accept_H1', 'continue'),
    n, llr, upper, lower after processing the supplied wins (bool iterable)."""
    if not (0 < p0 < p1 < 1 and 0 < alpha < 1 and 0 < beta < 1):
        raise ValueError("invalid parameters")
    upper = math.log((1 - beta) / alpha)
    lower = math.log(beta / (1 - alpha))
    llr = 0.0
    decision = "continue"
    n = 0
    for w in wins:
        n += 1
        if w:
            llr += math.log(p1 / p0)
        else:
            llr += math.log((1 - p1) / (1 - p0))
        if llr >= upper:
            decision = "accept_H1"
            break
        if llr <= lower:
            decision = "accept_H0"
            break
    return dict(decision=decision, n=n, llr=llr, upper=upper, lower=lower,
                p0=p0, p1=p1, alpha=alpha, beta=beta,
                wins=int(sum(1 for w in wins if w)), total=len(list(wins)))
