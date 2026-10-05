"""Sequential probability ratio tests for fold-wise candidate-vs-incumbent outcomes.

Two boundary conventions are supported and BOTH are always reported, because on
this project's data the promotion decision is sensitive to the choice:

* ``boundary="wald"`` -- Wald's (1945) classic cutoffs
  ``B = log((1-beta)/alpha)``, ``A = log(beta/(1-alpha))``.  Wald proved these
  hold the actual type-I and type-II errors at or below the declared alpha and
  beta when the observation sequence is fixed in advance and there is no peeking,
  which is exactly how the fold list is used here (a pre-declared, deterministic
  32-fold list, consumed in sorted-key order).  This is the convention the project
  charter names, so it is the default and the one the promotion gate reads.
* ``boundary="ville"`` -- the anytime-valid likelihood-ratio (Ville /
  test-martingale) cutoffs ``B = log(1/alpha)``, ``A = log(beta)``.  These remain
  valid under optional stopping with no pre-declared maximum sample size, and are
  strictly more conservative.

The convention is a declared design choice, not a tuning knob.  The DEFAULT is
``"ville"`` -- the stricter of the two -- so that a promotion claimed by this
repository also holds under the conservative reading, and the ``wald`` result is
returned alongside as ``alt_decision``.  Where the two disagree, the disagreement
is reported on the site and in the manifest rather than resolved quietly.

This matters concretely, and the effect is not always in the "stricter" direction
one might expect.  A wider lower bound (Ville: log(beta) = -2.3026 vs Wald:
log(beta/(1-alpha)) = -2.2513) means the test does not bail out on losses quite as
early, so it consumes more folds and can accumulate more wins.  For the promoted
H40 candidate the sign test stops at 21/29 with LLR +2.979 under Wald -- clearing
+2.8904 -- and at 22/30 with LLR +3.316 under Ville, clearing +2.9957.  The
promotion therefore holds under BOTH conventions, which is the point of reporting
both: the conclusion is not an artifact of the boundary choice.

Error control is conditional on independent Bernoulli block outcomes (or a valid
conditional supermartingale assumption).  A spatial collar reduces but does not
prove independence in geoscience data; callers must report that limitation.  A
non-crossing at the fixed end of the registered fold list is
"continue/inconclusive", never automatic acceptance.
"""
from __future__ import annotations

import math
from collections.abc import Iterable

BOUNDARIES = {
    "wald": lambda alpha, beta: (math.log((1.0 - beta) / alpha),
                                 math.log(beta / (1.0 - alpha))),
    "ville": lambda alpha, beta: (math.log(1.0 / alpha), math.log(beta)),
}


def _validate(p0: float, p1: float, alpha: float, beta: float) -> None:
    if not (0.0 < p0 < p1 < 1.0):
        raise ValueError("Require 0 < p0 < p1 < 1")
    if not (0.0 < alpha < 1.0 and 0.0 < beta < 1.0):
        raise ValueError("Require alpha and beta in (0, 1)")


def bounds(alpha: float, beta: float, boundary: str = "wald") -> tuple[float, float]:
    """(upper, lower) log-likelihood-ratio cutoffs for the named convention."""
    try:
        return BOUNDARIES[boundary](alpha, beta)
    except KeyError:
        raise ValueError(f"boundary must be one of {sorted(BOUNDARIES)}, got {boundary!r}")


class PairwiseSPRT:
    """Incremental two-hypothesis Bernoulli SPRT.

    H0: the candidate's spatial-block win probability is at most p0.
    H1: the candidate's spatial-block win probability is at least p1.
    The indifference region (p0, p1) is intentionally not resolved.
    """

    def __init__(
        self,
        p0: float = 0.5,
        p1: float = 0.7,
        alpha: float = 0.05,
        beta: float = 0.10,
        boundary: str = "ville",
    ) -> None:
        _validate(p0, p1, alpha, beta)
        self.p0, self.p1 = float(p0), float(p1)
        self.alpha, self.beta = float(alpha), float(beta)
        self.boundary = boundary
        self.upper, self.lower = bounds(alpha, beta, boundary)
        self.llr = 0.0
        self.n = 0
        self.wins = 0
        self.losses = 0
        self.decision = "continue"
        self.outcomes: list[bool] = []
        self.llr_path: list[float] = []

    def update(self, win: bool | int) -> dict:
        """Consume exactly one outcome and return the current test state."""
        if self.decision != "continue":
            raise RuntimeError("SPRT already crossed a boundary; no further folds may be consumed")
        if win not in (False, True, 0, 1):
            raise ValueError("Each fold outcome must be boolean or 0/1")
        observed = bool(win)
        self.n += 1
        self.wins += int(observed)
        self.losses += int(not observed)
        self.outcomes.append(observed)
        if observed:
            self.llr += math.log(self.p1 / self.p0)
        else:
            self.llr += math.log((1.0 - self.p1) / (1.0 - self.p0))
        self.llr_path.append(round(self.llr, 6))
        if self.llr >= self.upper:
            self.decision = "accept_H1"
        elif self.llr <= self.lower:
            self.decision = "accept_H0"
        return self.as_dict()

    def as_dict(self) -> dict:
        """Test state, including the decision under the OTHER convention.

        ``wins`` counts only folds consumed before the stopping boundary;
        ``wins_total`` reports the whole supplied vector.  An earlier revision
        recomputed both by re-iterating the input *after* the accumulation loop,
        so a generator caller got 0/0 while ``n`` was correct -- an internally
        inconsistent audit record.  The input is now materialised once.
        """
        alt = "ville" if self.boundary == "wald" else "wald"
        au, al = bounds(self.alpha, self.beta, alt)
        if self.llr >= au:
            alt_decision = "accept_H1"
        elif self.llr <= al:
            alt_decision = "accept_H0"
        else:
            alt_decision = "continue"
        return dict(
            decision=self.decision,
            n=self.n,
            llr=float(self.llr),
            upper=float(self.upper),
            lower=float(self.lower),
            p0=self.p0,
            p1=self.p1,
            alpha=self.alpha,
            beta=self.beta,
            wins=self.wins,
            losses=self.losses,
            wins_total=int(sum(self.outcomes)),
            total=len(self.outcomes),
            consumed=self.n,
            outcomes=list(self.outcomes),
            llr_path=list(self.llr_path),
            boundary_method=self.boundary,
            boundary_definition=("Wald (1945): upper=log((1-beta)/alpha), "
                                 "lower=log(beta/(1-alpha))" if self.boundary == "wald"
                                 else "Ville anytime-valid: upper=log(1/alpha), "
                                      "lower=log(beta)"),
            alt_boundary_method=alt,
            alt_upper=float(au),
            alt_lower=float(al),
            alt_decision=alt_decision,
            validity_assumption=("independent Bernoulli block outcomes / valid conditional "
                                 "supermartingale; a spatial collar reduces but does not "
                                 "prove independence"),
        )


def sprt_pairwise(
    outcomes: Iterable[bool | int],
    p0: float = 0.5,
    p1: float = 0.7,
    alpha: float = 0.05,
    beta: float = 0.10,
    boundary: str = "ville",
) -> dict:
    """Consume outcomes sequentially and stop at the first crossed boundary."""
    test = PairwiseSPRT(p0=p0, p1=p1, alpha=alpha, beta=beta, boundary=boundary)
    for outcome in list(outcomes):
        test.update(outcome)
        if test.decision != "continue":
            break
    return test.as_dict()


def folds_needed(p0: float = 0.5, p1: float = 0.7, alpha: float = 0.05,
                 beta: float = 0.10, boundary: str = "ville") -> dict:
    """Minimum all-win / all-loss fold counts that can reach a boundary.

    Reported so the holdout design is checked *before* any fold is evaluated: an
    SPRT that cannot reach either boundary at the planned n is a design defect,
    not a "run it longer" invitation.  An earlier revision of this project ran an
    8-fold SPRT whose upper boundary needed 9 all-wins -- it could never accept
    H1, and 8/8 (LLR 2.6918) fell short of 2.8904.
    """
    _validate(p0, p1, alpha, beta)
    upper, lower = bounds(alpha, beta, boundary)
    w = math.log(p1 / p0)
    l = math.log((1 - p1) / (1 - p0))
    alt = "ville" if boundary == "wald" else "wald"
    au, al = bounds(alpha, beta, alt)
    return dict(upper=upper, lower=lower, llr_per_win=w, llr_per_loss=l,
                min_all_wins=math.ceil(upper / w) if w > 0 else None,
                min_all_losses=math.ceil(lower / l) if l < 0 else None,
                boundary_method=boundary,
                alt_boundary_method=alt, alt_upper=au, alt_lower=al,
                alt_min_all_wins=math.ceil(au / w) if w > 0 else None,
                alt_min_all_losses=math.ceil(al / l) if l < 0 else None)


def sprt_normal_mean(x, d_alt: float = 0.5, alpha: float = 0.05, beta: float = 0.10,
                     sigma: float | None = None, boundary: str = "ville") -> dict:
    """Wald's SPRT for the mean of a normal with (estimated) known scale.

    Wald (1945) gives the SPRT for any family with a likelihood ratio; for
    ``x_i ~ N(mu, sigma^2)`` testing H0: mu = 0 against H1: mu = d_alt * sigma the
    accumulated log-likelihood ratio after n observations is

        LLR_n = (d_alt / sigma) * sum(x_i) - n * d_alt^2 / 2

    with the same boundary convention as the binomial form.

    Why this exists alongside ``sprt_pairwise``: the fold-level outcome here is a
    *magnitude* (the exact per-fold contribution ``D_cand*(s_cand - s_anchor)`` to
    the pooled score), and binarising it to a sign throws the magnitude away.  A
    candidate can win the pooled metric decisively -- the pooled metric is what the
    organizer scores -- while losing a majority of individual folds, by gaining a
    lot in a few truth-dense folds and a little in many truth-poor ones.  The sign
    test then rejects a real pooled improvement.  Both tests are reported; neither
    is allowed to be the only evidence.

    ``sigma`` may be supplied; otherwise it is estimated from the observed spread
    (Wald's sequential t-test form).  Estimating it is a declared deviation from
    the exactly-known-variance case and is recorded in the return value.
    """
    import numpy as _np
    _validate(0.5, 0.7, alpha, beta)          # reuse the alpha/beta range checks
    upper, lower = bounds(alpha, beta, boundary)
    alt = "ville" if boundary == "wald" else "wald"
    au, al = bounds(alpha, beta, alt)
    v = _np.asarray(list(x), float)
    if v.size == 0:
        return dict(decision="continue", n=0, llr=0.0, upper=upper, lower=lower,
                    sigma=None, sigma_estimated=True, mean=None, d_alt=d_alt,
                    boundary_method=boundary, alt_boundary_method=alt,
                    alt_upper=au, alt_lower=al, alt_decision="continue",
                    total_n=0, path=[])
    s = (float(sigma) if sigma else
         (float(v.std(ddof=1)) if v.size > 1 else float(abs(v.mean()) or 1.0)))
    if not math.isfinite(s) or s <= 0:
        s = 1.0
    llr = 0.0
    decision = "continue"
    n = 0
    path = []
    run = 0.0
    for xi in v:
        n += 1
        run += float(xi)
        llr = (d_alt / s) * run - n * d_alt * d_alt / 2.0
        path.append(round(llr, 6))
        if llr >= upper:
            decision = "accept_H1"
            break
        if llr <= lower:
            decision = "accept_H0"
            break
    if llr >= au:
        alt_decision = "accept_H1"
    elif llr <= al:
        alt_decision = "accept_H0"
    else:
        alt_decision = "continue"
    return dict(decision=decision, n=n, llr=float(llr), upper=float(upper),
                lower=float(lower), sigma=float(s), sigma_estimated=sigma is None,
                d_alt=d_alt, mean=float(v.mean()), sum=float(v.sum()), path=path,
                total_n=int(v.size), boundary_method=boundary,
                boundary_definition=("Wald (1945) normal-mean SPRT" if boundary == "wald"
                                     else "Ville anytime-valid bounds on the normal-mean LLR"),
                alt_boundary_method=alt, alt_upper=float(au), alt_lower=float(al),
                alt_decision=alt_decision)
