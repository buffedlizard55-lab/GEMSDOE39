"""Sequential Probability Ratio Test for fold-wise candidate-vs-incumbent wins.

The outcome is Bernoulli: a pre-registered candidate wins a spatial block if
its DTI is strictly larger than the locked incumbent's DTI; ties are losses.
For the simple hypotheses p=p0 versus p=p1, the log-likelihood ratio is updated
one block at a time. The boundaries used here are the conservative, anytime-
valid likelihood-ratio boundaries upper=log(1/alpha), lower=log(beta), rather
than Wald's common approximate (1-beta)/alpha and beta/(1-alpha) cutoffs.

Error control is conditional on independent Bernoulli block outcomes (or a
valid conditional supermartingale assumption). A spatial collar reduces but
does not prove independence in geoscience data; callers must report that
limitation. A non-crossing at the fixed end of the registered fold list is
"continue/inconclusive", never automatic acceptance.
"""
from __future__ import annotations

import math
from collections.abc import Iterable


def _validate(p0: float, p1: float, alpha: float, beta: float) -> None:
    if not (0.0 < p0 < p1 < 1.0):
        raise ValueError("Require 0 < p0 < p1 < 1")
    if not (0.0 < alpha < 1.0 and 0.0 < beta < 1.0):
        raise ValueError("Require alpha and beta in (0, 1)")


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
    ) -> None:
        _validate(p0, p1, alpha, beta)
        self.p0, self.p1 = float(p0), float(p1)
        self.alpha, self.beta = float(alpha), float(beta)
        self.upper = math.log(1.0 / self.alpha)
        self.lower = math.log(self.beta)
        self.llr = 0.0
        self.n = 0
        self.wins = 0
        self.losses = 0
        self.decision = "continue"
        self.outcomes: list[bool] = []

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
        if self.llr >= self.upper:
            self.decision = "accept_H1"
        elif self.llr <= self.lower:
            self.decision = "accept_H0"
        return self.as_dict()

    def as_dict(self) -> dict:
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
            outcomes=list(self.outcomes),
            boundary_method="Ville-safe likelihood-ratio bounds",
            validity_assumption="independent Bernoulli block outcomes / valid conditional supermartingale",
        )


def sprt_pairwise(
    outcomes: Iterable[bool | int],
    p0: float = 0.5,
    p1: float = 0.7,
    alpha: float = 0.05,
    beta: float = 0.10,
) -> dict:
    """Consume outcomes sequentially and stop at the first crossed boundary."""
    test = PairwiseSPRT(p0=p0, p1=p1, alpha=alpha, beta=beta)
    for outcome in outcomes:
        test.update(outcome)
        if test.decision != "continue":
            break
    return test.as_dict()
