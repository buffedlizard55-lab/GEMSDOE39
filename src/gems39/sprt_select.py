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
    n, llr, upper, lower after processing the supplied wins (bool iterable).

    Boundaries are Wald's (1945) with alpha = P(accept H1 | H0) and
    beta = P(accept H0 | H1):

        B = log((1 - beta) / alpha)   accept H1 when llr >= B
        A = log(beta / (1 - alpha))   accept H0 when llr <= A

    DEFECT FIXED 2026-10-05: the previous revision computed ``wins`` and
    ``total`` by re-iterating ``wins`` *after* the accumulation loop.  When the
    caller passed a generator the loop had already exhausted it, so both fields
    were reported as 0 while ``n`` was correct -- an internally inconsistent
    audit record.  The input is now materialised once, and ``wins`` counts only
    the folds actually consumed before the stopping boundary (``wins_total``
    reports the whole vector for transparency).
    """
    if not (0 < p0 < p1 < 1 and 0 < alpha < 1 and 0 < beta < 1):
        raise ValueError("invalid parameters")
    seq = [bool(w) for w in wins]
    upper = math.log((1 - beta) / alpha)
    lower = math.log(beta / (1 - alpha))
    llr = 0.0
    decision = "continue"
    n = 0
    wins_used = 0
    path = []
    for w in seq:
        n += 1
        if w:
            llr += math.log(p1 / p0)
            wins_used += 1
        else:
            llr += math.log((1 - p1) / (1 - p0))
        path.append(round(llr, 6))
        if llr >= upper:
            decision = "accept_H1"
            break
        if llr <= lower:
            decision = "accept_H0"
            break
    return dict(decision=decision, n=n, llr=llr, upper=upper, lower=lower,
                p0=p0, p1=p1, alpha=alpha, beta=beta,
                wins=wins_used, wins_total=int(sum(seq)), total=len(seq),
                consumed=n, llr_path=path)


def folds_needed(p0=0.5, p1=0.7, alpha=0.05, beta=0.10):
    """Minimum all-win / all-loss fold counts that can reach a Wald boundary.

    Reported so the design of the holdout is checked *before* any fold is
    evaluated: an SPRT that cannot reach either boundary at the planned n is a
    design defect, not a "run it longer" invitation.
    """
    upper = math.log((1 - beta) / alpha)
    lower = math.log(beta / (1 - alpha))
    w = math.log(p1 / p0)
    l = math.log((1 - p1) / (1 - p0))
    return dict(upper=upper, lower=lower, llr_per_win=w, llr_per_loss=l,
                min_all_wins=math.ceil(upper / w) if w > 0 else None,
                min_all_losses=math.ceil(lower / l) if l < 0 else None)


def sprt_normal_mean(x, d_alt=0.5, alpha=0.05, beta=0.10, sigma=None):
    """Wald's SPRT for the mean of a normal with (estimated) known scale.

    Wald (1945) gives the SPRT for any family with a likelihood ratio; for
    x_i ~ N(mu, sigma^2) testing H0: mu = 0 against H1: mu = d_alt * sigma the
    accumulated log-likelihood ratio after n observations is

        LLR_n = (d_alt / sigma) * sum(x_i) - n * d_alt^2 / 2

    with the same boundaries as the binomial form,
    B = log((1-beta)/alpha) and A = log(beta/(1-alpha)).

    Why this exists alongside ``sprt_pairwise``: the fold-level outcome here is a
    *magnitude* (the first-order change in the pooled score contributed by that
    fold), and binarising it to a sign throws the magnitude away.  A candidate can
    win the pooled metric decisively -- the pooled metric is what the organizer
    scores -- while losing a majority of individual folds, by gaining a lot in a
    few truth-dense folds and a little in many truth-poor ones.  The sign test
    then rejects a real pooled improvement.  Both tests are reported; neither is
    allowed to be the only evidence.

    ``sigma`` may be supplied; otherwise it is estimated from the observed spread
    (Wald's sequential t-test form).  Estimating it is a declared deviation from
    the exactly-known-variance case and is recorded in the return value.
    """
    import numpy as _np
    v = _np.asarray(list(x), float)
    if v.size == 0:
        return dict(decision="continue", n=0, llr=0.0,
                    upper=math.log((1 - beta) / alpha),
                    lower=math.log(beta / (1 - alpha)), sigma=None,
                    sigma_estimated=True, mean=None, d_alt=d_alt)
    s = float(sigma) if sigma else float(v.std(ddof=1)) if v.size > 1 else float(abs(v.mean()) or 1.0)
    if not math.isfinite(s) or s <= 0:
        s = 1.0
    upper = math.log((1 - beta) / alpha)
    lower = math.log(beta / (1 - alpha))
    llr = 0.0
    decision = "continue"
    n = 0
    path = []
    run = 0.0
    for xi in v:
        n += 1
        run += xi
        llr = (d_alt / s) * run - n * d_alt * d_alt / 2.0
        path.append(round(llr, 6))
        if llr >= upper:
            decision = "accept_H1"
            break
        if llr <= lower:
            decision = "accept_H0"
            break
    return dict(decision=decision, n=n, llr=llr, upper=upper, lower=lower,
                sigma=s, sigma_estimated=sigma is None, d_alt=d_alt,
                mean=float(v.mean()), sum=float(v.sum()), path=path,
                total_n=int(v.size))
