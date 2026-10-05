# GEMSDOE39 — research brief and decision register

**Date:** 2026-10-05. This is a hypothesis and measurement register. Scores attributed to
previous submissions in the project brief are **owner-reported**; no organiser receipt ties a
raster hash to a score, and nothing here is promoted from a proxy measurement to a leaderboard
claim.

---

## 1. What was measured before anything was built

Thirty-nine historical artifacts were restored from the sibling repositories and re-opened from
their bytes. Thirty-one of them carry an owner-reported score
([`registry/score_ledger.csv`](../registry/score_ledger.csv)). On those 30+ artifacts
([`registry/forensics_instruments.json`](../registry/forensics_instruments.json),
[`registry/instrument_calibration.json`](../registry/instrument_calibration.json)):

| statistic | Spearman vs owner-reported live score | p | n |
|---|---:|---:|---:|
| emitted pixel count | **−0.686** | < 1e-4 | 30 |
| fraction of dots on the known catalogue | **−0.598** | 5e-4 | 30 |
| DTI against the visible catalogue | −0.344 | 0.062 | 30 |
| DTI against off-catalogue USGS SGMC faults | −0.080 | 0.674 | 30 |
| catalogue-hidden spatially-blocked holdout DTI | −0.104 | 0.585 | 30 |

Two conclusions drive everything downstream:

1. **Emitted mass is the dominant controllable variable.** Score falls monotonically from 0.2708 at
   40,199 px to 0.03 at 517 k px.
2. **No local proxy ranks detector fields.** Every candidate instrument sits at |ρ| ≤ 0.37, and the
   catalogue-DTI is *negatively* correlated — it rewards re-predicting faults that are already mapped and
   therefore masked. GEMSDOE29 reports ρ = +0.709 for a catalogue-hidden proxy on 10 artifacts; on our 30
   the equivalent design gives −0.104. The calibration sets differ and neither is significant. This
   discordance is recorded, not resolved.

## 2. The metric, reduced to a decision rule

From the organiser's equations
([page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric)):

```
DTI = TP_w / (TP_w + 0.2*FP_w + 0.8*FN_w)
```

Adding one predicted pixel with expected kernel credit `w`: `ΔTP = w`, `ΔFP = 1−w`, `ΔFN = −w`, so the
denominator changes by `w + 0.2(1−w) − 0.8w = 0.2` — a constant. The pixel is worth adding iff

```
w > 0.2 * DTI          # 0.0556 at DTI = 0.2778
```

This is the **break-even rule** and it is exact. GEMSDOE32 arrived at the same bar empirically (0.0548)
and analytically (0.2 × 0.26 = 0.0520); our route is a third, derived directly from the published
formula.

Fitting `TP = A·n^γ` to the one clean within-lineage pair (60,069 px → 0.2477, 44,090 px → 0.2600) and
imposing `γ·TP/n = 0.2·DTI`:

| assumed truth size K | fitted γ | break-even n* | DTI at n* |
|---:|---:|---:|---:|
| 15,000 px | 0.305 | ≈ 26,800 | 0.282 |
| 30,000 px | 0.144 | ≈ 20,900 | 0.286 |
| 45,000 px | 0.066 | ≈ 13,400 | 0.292 |
| 60,000 px | 0.020 | ≈ 5,400 | 0.302 |

The DTI surface is nearly flat between roughly 13 k and 45 k px (0.281–0.289 at n = 30,000 under every
K). **Budget is a second-order decision; field quality is first-order**, because DTI is very nearly
linear in TP at fixed budget. The shipped budget is **24,000 px**, inside the flat region for every
assumed K and 36 % below the historical champion's pixel count.

## 3. Ranked hypotheses (registered before implementation)

Each names the layers, the physical signature, why it should catch a fault the catalogue lacks rather
than one it contains, how it differs from anything already in this repository family, expected DTI
direction and implementation cost.

| # | hypothesis | layers | signature | why off-catalogue | novelty | expected DTI / cost |
|---|---|---|---|---|---|---|
| 1 | **H39-A** blocked off-catalogue discriminant | all 19 bands + 21 derived transforms | supervised P(Quaternary fault), hard negatives, blocked OOF | trained on the catalogue but applied **only** where no mapped fault exists | prior supervised models in this family leaked `dist_to_catalogue` (train AUC 1.0; 65.95 % of dots within 300 m of a mapped fault) | **High** (measured 0.274 → 0.670 lift). Medium cost — **SHIPPED** |
| 2 | **H39-B** off-catalogue-fault discriminant | same stack | P(fault), trained on 62,703 px of USGS SGMC faults > 300 m from the catalogue | positives are literally faults the catalogue does not contain | no prior model in the family trains on an off-catalogue fault set | Moderate (lift 0.589). Low cost — built as channel B |
| 3 | **H39-1** range-front suppression | `det_elev`, `det_elev_slope` | large-scale range-front gradient field, then the residual lineament that survives it | mapped Quaternary catalogues are dominated by range-front scarps; the missing faults are disproportionately intra-basin / piedmont / antithetic | every prior detector in the family *amplifies* the strongest lineament; this inverts the weighting | Moderate, unproven (lift 0.627 as a channel). Low cost |
| 4 | **H39-C/D** cross-physics coherence and magnetic worms | `rtp`, `tmi`, `iso_grav_anom`, `cond_surf`, `tc` | directional 3-physics cross-gradient coherence; upward-continued analytic-signal ridge persistence | a buried fault can show a co-located edge in three independent physics with no scarp at all | prior work fuses scalar edge *amplitudes*; these use gradient *direction* and cross-scale amplitude ratio | Low–moderate (lift 0.645 / 0.634). Low cost |
| 5 | **H39-X** flight-line residual test | GeoDAWN Area-2 raw magnetic/radiometric CSV | does a candidate lineament reproduce across adjacent 400 m traverses? | removes survey-parallel artefacts that every gridded detector inherits | nobody in the family has tested it | Unknown. **BLOCKED** — needs official USGS ScienceBase bytes (3.74 GB / 427 MB); obtainability **not** verified |

**Selection rule, preregistered:** promote the candidate that the Wald SPRT accepts against the
topographic × magnetic control and that has the highest holdout lift. Result: H39-A
(12/12 cell wins, LLR +3.028 ≥ upper boundary 2.890 → accept H1, lift 0.6698 vs 0.2735).
All five results are reported in the build manifest, not just the winner.

## 4. Sequential test (Wald, 1945)

Declared before any fold was read:

```
alpha = 0.05    beta = 0.10    p0 = 0.50    p1 = 0.70
upper = log(1/alpha)  = +2.995732  ->  accept H1, stop   (anytime-valid / Ville bound)
lower = log(beta)     = -2.302585  ->  accept H0, stop
```

The repository's `sprt_select.py` uses the conservative anytime-valid likelihood-ratio bounds rather
than Wald's common approximations `log((1−β)/α) = +2.890` and `log(β/(1−α)) = −2.251`. The former are
strictly harder to cross, so the reported acceptance is conservative in both directions.

Twelve independent cells are available: 2 × 2 spatial blocks × 3 draws, hiding 25 % of catalogue
**connected components** with a 2 px collar. Each cell votes on whether the candidate's mean
percentile-rank lift on the hidden truth beats the control's; ties are losses. Each win adds
`log(0.7/0.5) = +0.3365` to the log-likelihood ratio.

H39-A won the first nine cells; the LLR reached **+3.028 ≥ +2.995732 on the 9th**, so the test
**stopped by its own rule** with three cells left unused. It did not keep consuming cells after the
boundary was crossed, which is exactly the "run it a bit longer" failure mode this design exists to
prevent. All five candidates were scored; only the highest-lift accepted one (H39-A, lift 0.6698 vs
0.2735 control) was promoted.

The SPRT controls the sequential-stopping error for one declared pairwise comparison. It does **not**
correct for searching across a portfolio of candidates; that is handled by preregistering the five
candidates with their layers and rationale before any of them was scored, and reporting all five.

## 5. Verified-source register

| claim | source |
|---|---|
| metric, α = 0.2, β = 0.8, R = 300 m, submission format | [DrivenData #306 problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) |
| known USGS/INGENIOUS fault pixels are masked from evaluation | [DrivenData community #11516, staff reply](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2) |
| rules: submission limits, external-data licence, AI disclosure | [NREL/DOE PDF 96647](https://docs.nlr.gov/docs/fy26osti/96647.pdf) · [rules page](https://www.drivendata.org/competitions/306/competition-doe-gems/rules/) |
| GeoDAWN provenance | [USGS GeoDAWN](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) · DOI [10.5066/P93LGLVQ](https://doi.org/10.5066/P93LGLVQ) · [ScienceBase item](https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7) |
| off-catalogue fault truth (channel B) | USGS SGMC faults, 62,703 px > 300 m from the given catalogue |
| geothermal wells and springs (available, unused) | [OpenEI GDR 1391](https://gdr.openei.org/submissions/1391) · DOI [10.15121/1881483](https://doi.org/10.15121/1881483) (CC BY 4.0) |
| 1 m DEM (used via the derived scarp layer) | [USGS 3DEP](https://www.usgs.gov/3d-elevation-program) |
| sequential test | Wald, A. (1945), *Sequential Tests of Statistical Hypotheses*, Ann. Math. Statist. 16(2) |
| reference solution | [drivendataorg/gems-prize-reference-solution](https://github.com/drivendataorg/gems-prize-reference-solution) |

The brief also prints `docs.nlr.gov`; the domain in the live link is `docs.nrel.gov`. `docs.nlr.gov`
resolves to the same document, but the spelling in the brief is treated as an irregularity and both are
recorded.

## 6. Irregularities

- **IR-39-01 — the sample submission is not what the page says.** The mirrored `sample_submission.tif`
  is described as "total fault absence" but contains 60,988 pixels equal to `1.0`, coinciding exactly
  with `labels.tif == 1`. Used only as a grid/footprint template; the footprint it implies (5,167,373 px)
  is independently confirmed by `labels.tif >= 0` (5,106,385 zeros + 60,988 ones), so the geometry is
  safe even though the file's semantics are not as documented.
- **IR-39-02 — `existing_faults.tif` and `labels.tif` are byte-identical** (sha256 `7ba308cc…`).
  `labels.tif` is treated as the single authoritative catalogue.
- **IR-39-03 — instrument discordance.** GEMSDOE32 reports the SGMC instrument at ρ = +0.54 (n = 12);
  we measure −0.080 (n = 30). Calibration sets differ; neither is significant.
- **IR-39-05 — the 0.2778 attribution is unverified.** The project brief credits 0.2778 to
  `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros`, but the GEMSDOE32 site marks that file "UNSCORED", and
  the public leaderboard shows a 0.2778 row with no filename attached. We treat 0.2778 as an official
  board row, not as a verified score for that specific file, and we say so on the site.
- **IR-39-06 — the brief's rules link is misspelled.** It prints `docs.nlr.gov`; the live host is
  `docs.nrel.gov`. Both currently resolve to the same document; both spellings are recorded.
- **IR-39-07 — two sessions edited this repository concurrently.** `main` advanced past this branch's
  base with three merges from `arena/01a109d9-gemsdoe39` (an H39X-01 strain-corridor candidate that was
  audited and *not* approved). The merge changed `grid.write_submission` to reject anything but
  `outside="nan"`, `validate_submission.py` to require NaN-outside, `emission.emit_fast` to stop
  under-filling its budget, and `sprt_select` to use Ville bounds. We merged that work in and adapted to
  it rather than reverting it, adding `src/gems39/io39.py` (a writer that emits **both** encodings,
  because the brief records a real portal rejection caused by the NaN-outside encoding) and a
  `--outside {auto,zeros,nan}` mode in the validator so each encoding is judged by the criteria that
  actually apply to it.

- **IR-39-04 (fixed this session) — manifest serialiser dropped nested per-cell tables.** The previous
  revision's filter kept only keys named `dti/tp/fp/fn/n_truth`, and fold names are not among them, so
  `catalogue_hidden_folds` and `catalogue_hidden_per_quadrant` were written as `{}`. The current
  manifest is `schema_version: 2` and stores the full per-cell table.

## 7. Current limitations

- The repository has no DrivenData credentials. Competition rasters are restored from owner-maintained
  GitHub mirrors and SHA-256 pinned; hash agreement proves mirror integrity, not organiser
  authentication.
- The private truth is inaccessible, so no field metric here is the scored objective.
- The public leaderboard scores only a public chunk of the new faults; the Initial Prize Round scores a
  private chunk; the Final Prize Round re-scores the same file against an expanded label set. Optimising
  hard for the public number is not the same as optimising for the prize.
- Hypothesis 5 is blocked on data, not on ideas, and its obtainability is deliberately unverified.

---

## 8. H39Y-01 preregistered portfolio and observed result (2026-10-05)

A separate, single-candidate cycle registered four geological hypotheses before implementation or H39Y scoring. The first-ranked and only candidate tested was **H39Y-01**, a co-located GeoDAWN `ThK`, `UK`, or `UTh` ratio-gradient corroborated by the `cond_surf` gradient. Its frozen detector was `sqrt(R*C) * abs(cos(theta_R-theta_C))`, with 2-pixel Gaussian gradients and footprint quantiles 0.02/0.995. Geological plausibility was informed by the official [USGS GeoDAWN release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) and [USGS Fact Sheet 2020-3055](https://pubs.usgs.gov/fs/2020/3055/fs20203055.pdf); these sources do not demonstrate candidate skill. The aligned/quantized GeoDAWN TIFF is an owner mirror, not raw USGS bytes retrieved for this run.

The frozen 8×8 spatial plan yielded 11 eligible even/even cells. The H39-A-model comparator used 40 checked-in fields, fold-specific training-only catalogue masks and negatives, and the existing HistGradientBoosting hyperparameters. Candidate and comparator used matched per-cell Poisson-disk emission (2.8-pixel spacing, 2-pixel prediction buffer, prorated 24,000-pixel budget). **H39Y-01 lost all first five eligible cells;** its SPRT crossed the lower boundary after tile 5 (`LLR = −2.554128`, lower boundary `−2.302585`, `accept_H0`). Pooled local proxy DTI over those five tiles was 0.035772 for H39Y-01 versus 0.168755 for H39-A-model. The candidate is rejected for this holdout; the incumbent is retained, no H39Y TIFF was created, and **no weekly submission slot is recommended**. No later eligible cell was fit or scored.

The scoring-mask details were amended before testing to match the [DrivenData staff clarification](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4): exact known-fault pixels are masked, but adjacent predictions are scored normally. The 2-pixel collar applies only to emission. Two pre-score addenda also specify fold seeds and remove hidden labels from negative-pool construction. The full plan, inputs, and results are archived in [`preregistration-h39y-20261005.md`](preregistration-h39y-20261005.md), [`preregistration-addendum-h39y-20261005.md`](preregistration-addendum-h39y-20261005.md), [`preregistration-addendum2-h39y-20261005.md`](preregistration-addendum2-h39y-20261005.md), [`h39y01-preflight-20261005.json`](reports/h39y01-preflight-20261005.json), and [`h39y01-validation-20261005.md`](reports/h39y01-validation-20261005.md).

**Reporting irregularity:** the evaluator correctly stopped at the H0 boundary, then its final report assembly raised a duplicate-key `TypeError`. The per-tile outcomes had already been atomically saved. The report was recovered from those outcomes and the committed pre-score plan without refitting, rerunning tile predictions, or repeating per-tile DTI scoring. Pooled DTI was re-aggregated from the saved per-tile TP/FP/FN counts. The start timestamp and candidate-field digest were not saved and are marked unavailable. This does not change the observed SPRT outcome, but the report-generation defect is recorded.

### Next steps after H39Y-01

- Do **not** tune or score H39Y-01 again on these exposed cells. A new recipe needs a fresh holdout or a predeclared family-wise alpha allocation.
- Keep the current H39-A-model/R39A incumbent as the downloadable artifact; this H39Y-01 result is not a leaderboard score and does not authorize spending a weekly slot.
- H39Y-02 (thermal-fluid well/spring geochemistry) and H39Y-03 (volcanic-vent alignment with magnetic lineaments) remain untested. H39Y-04 remains deferred until the official shallow-probe archive bytes and schema are available.
- Do not infer no geological relationship from this result. It rejects only the frozen detector relative to the incumbent under this public-catalogue proxy and the stated sequential model.
