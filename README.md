# GEMSDOE39 — Unique Submission for the DOE GEMS Prize (DrivenData #306)

> **Mission: Maximize P(Win). Own the Outcome.**
> Produce a *unique*, format-valid GeoTIFF submission for DrivenData competition
> **#306 — The Geologic Enhanced Mapping System (GEMS) Prize Challenge**, gated by a
> pre-declared **Wald sequential probability ratio test (1945)**, priced against a
> **calibration derived from organizer-returned scores**, and delivered with a
> 1-click download from GitHub Pages.

**Read this file at the start of every session.** The full standing request is
reproduced verbatim in [§ Permanent project charter](#permanent-project-charter-read-before-every-session).

---

## DOWNLOAD THE SUBMISSION

| | |
|---|---|
| **Primary artifact (submit this one)** | [`docs/downloads/gemsdoe39-h40-pfpt-20261005T030000Z-zeros.tif`](docs/downloads/gemsdoe39-h40-pfpt-20261005T030000Z-zeros.tif) |
| Site with the big yellow button | <https://buffedlizard55-lab.github.io/GEMSDOE39/> |
| How to upload, step by step | [docs/executive-summary.html](docs/executive-summary.html) |
| Submission **Name** | `gemsdoe39-h40-pfpt-playfairway-permeability` |
| Submission **Note** | see § Submission name and note below |
| NaN-outside twin (identical pixels) | `docs/downloads/gemsdoe39-h40-pfpt-20261005T030000Z-nan.tif` |
| Full machine-readable receipt | `docs/downloads/gemsdoe39-h40-pfpt-20261005T030000Z-manifest.json` |

Both twins carry **exactly the same predicted pixels**; they differ only in how
the area outside the study footprint is encoded (`0.0` vs `NaN`). Both satisfy the
organizer's format rule, *"data outside the bounds is null or nan"* and *"a single
layer of float32 values between 0 and 1"*
([page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)).

**Submit the `-zeros` twin.** Reason, from evidence rather than preference: the
highest-scoring artifact in this project's family (`h33-2-b2`, **0.2778**) is the
`-zeros` twin, and zeros-outside cannot trigger the
`Predicted values must be in range [0, 1]` validator error under any uploader
that mishandles NaN.

---

## Status (2026-10-05, session 2 — H40 round)

* **Unique artifact generated.** No pixel is copied from any previous submission.
  The live-scored artifacts are read only to (a) verify the nested pair that
  identifies `|G|` and (b) serve as a comparison arm on local instruments.
* **Format:** single-band float32 GeoTIFF, EPSG:32611, 3730 × 3292 @ 100 m,
  every in-footprint pixel finite and in `[0, 1]`, zero pixels on the masked
  catalogue.
* **`scripts/validate_submission.py` on both twins:** PASS (all 12 checks).
* **Live leaderboard read by hand on 2026-10-05**
  ([leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)):
  **#1 nchuzhoy 0.3262**, #2 kinghorton42 0.3222, #3 DARD 0.3195, …
  **#13 extradr19 0.2778** (this family's best), #14/#15 0.2708.
* **No score is claimed for the H40 artifact.** What is claimed is a *price* —
  a model-based prediction with a stated derivation and stated nuisance priors.

---

## Permanent project charter (read before every session)

Reproduced verbatim from the standing request:

```text
Review the repo.

MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION. DO NOT COPY A PREVIOUS
SUBMISSION UNLESS IT'S FOR LEARNING AND EDUCATION. BUT WE MUST GENERATE A UNIQUE
TIF SUBMISSION.

There should be an easy-to-download submission tif file as described by the
prompt. Read the entire prompt.

Use a formally valid sequential test so "run it a bit longer" doesn't quietly
inflate false positives. Wald's sequential probability ratio test (1945) is the
classical, formally valid answer: define error rates in advance, and let the
stopping decision be a principled test rather than a judgment call. Apply it
directly to holdout evaluation -- accumulate the log-likelihood ratio between
"this candidate beats the current best" and "it doesn't," fold by fold, and
stop only when it crosses one of the two pre-declared boundaries.

Before you spend a submission slot: generate 3-5 new candidate geological
hypotheses (not tried before). Each must name the layer(s) it uses, the physical
signature it targets, why that signature catches a fault the USGS/INGENIOUS
catalogue is missing, and how it differs from anything already implemented in
this repo. Rank them by expected DTI improvement against implementation cost.
Validate the top candidate on a spatially-blocked holdout first. If a hypothesis
needs new external data, name the specific free official source and confirm it
is obtainable.

Study the highest-scoring submission across the listed GEMSDOE...GEMSDOE39 sites
(h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros = 0.2778) and explain why/how it
scored highest, and whether we can beat it. The leaderboard target is above the
current best (0.3195 was quoted; the live page now shows 0.3262).

Put this prompt into the repo README and read it every session.

Work line by line from official verified trusted sources, provide links for
manual review, no manual input required from the user, flag irregularities,
NO hallucinations.

Run three passes: (1) implement and verify, (2) review for bugs, missing
requirements and edge cases and fix them, (3) re-check the whole implementation
against the original request and improve. Do not stop after pass 1.

Finally: create a pull request and merge it onto main, plus suggestions for
remaining work and the limitations blocking success.
```

---

## Why the previous best (0.2778) scored what it scored

Full derivation, links and sensitivity analysis:
[docs/research/h40-hypotheses.md § 1](docs/research/h40-hypotheses.md).

The artifact is `gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif`:
**37,654 dots**, Poisson-disk at 2.8 px, every dot more than 200 m from the
catalogue. It is the 40,199-dot `h27-4-r1-solo-d2-8` artifact (0.2708) with
**exactly 2,545 dots deleted** — every dot at 1.414–2.000 px from the catalogue.
Verified by re-reading both GeoTIFFs (`calibrate.verify_nested_pair`, pinned by
`tests/test_h40.py::test_nested_pair_is_exactly_as_documented`).

Because `TP_w + FN_w = |G|`, the metric collapses to
`DTI = TP_w / (0.2·TP_w + 0.2·FP_w + 0.8·|G|)`, and that nested pair identifies
`|G|`, the number of *active* scored-truth pixels:

```
|G| = 14,143            [13,987 .. 14,434] across a 5x4 (l0, kb) sensitivity grid
anchor hit rate h = 6.130%   (2,308 of 37,654 dots land within 300 m of scored truth)
forward model reproduces 0.2778 exactly, and 0.2707 vs the observed 0.2708
```

Denominator at the 0.2778 artifact: `0.2·TP_w` = 1,095 (5.6%) ·
`0.2·FP_w` = 7,300 (37.0%) · **`0.8·|G|` = 11,314 (57.4%)**.

So **57% of the denominator is a fixed cost nobody can reduce**, and of the
controllable 43%, false-positive mass costs 6.7× what true-positive mass
contributes. Break-even marginal hit probability `π* = 3.149%`.

**Answer to "can we beat it?"** Yes, and the calibration says how much is needed.
Reaching 0.3262 (live #1) at the anchor's FP mass requires `TP_w ≈ 6,489` instead
of 5,470 — a **+18.6% improvement in hit rate** (6.13% → ~7.3%). Reaching it at
the anchor's `TP_w` instead requires cutting `FP_w` from 36,379 to 21,770. Both
are the same lever: **more hits per dot.** Budget tuning cannot do it — the priced
optimum for the incumbent surface is ~25–28k dots for ~0.283.

Also settled analytically and pinned by
`tests/test_h40.py::test_binary_emission_is_never_worse_than_scaling_down`:
`0.8·|G|` is a constant additive term, so `DTI(c·p)` rises monotonically in `c`
and the inclusion rule is independent of `p`. **Binary {0,1} emission at full
scale is provably optimal; submitting soft probabilities is strictly worse.**

---

## The seven H40 hypotheses

Full register with layers, signatures, literature links, novelty statements and
rank justification: [docs/research/h40-hypotheses.md](docs/research/h40-hypotheses.md).

| Rank | ID | Hypothesis | Layers | Why it finds a *missing* fault | Novelty here |
|---|---|---|---|---|---|
| 1 | **H40-F** | Supervised propensity for faults *absent from the catalogue* | all 19 bands (σ=0, 3 px) + 7 external layers; positives = USGS SGMC faults >300 m from the catalogue | Training on the catalogue learns "what a *mapped* fault looks like" and over-weights geomorphic expression — the very property that got it mapped. Training on catalogue-*absent* faults inverts the bias. | Prior rounds used SGMC-off as an **evaluation instrument only**; this uses it as the **training target**, out-of-fold on the same block grid. |
| 2 | **H40-A** | Critically stressed lineament | USGS **Siler (2022)** slip/dilation tendency, DOI [10.5066/P9YL58W6](https://doi.org/10.5066/P9YL58W6) (37,811 in-footprint segments) + 7 competition bands | USGS built this dataset to "identify faults … likely to host **as-yet-undiscovered hydrothermal processes**" — the same population this competition scores. | Orientation vs a *measured local* stress field, using an **empirical** favourability curve (peak/trough = 6.9×). All prior work was orientation-blind amplitude. |
| 3 | **H40-C** | Deep reservoir temperature anomaly | GDR **1391** INGENIOUS, DOI [10.15121/1881483](https://doi.org/10.15121/1881483) — 27,092 well/spring rows, 1,244 with a fluid geothermometer | A hot reservoir cannot exist without a deep permeable pathway. Geothermometry measures the fluid that travelled the conduit; surface heat flow smears it. | First use of a *reservoir-depth* geothermometer as a primary targeting field. |
| 4 | **H40-G** | Tip-continuation / linkage corridor | catalogue geometry (skeleton tips + local strike) × H40-A fields | Staff: a new fault "can include newly mapped geometry of an existing fault system" ([11536](https://community.drivendata.org/t/11536/2)). A trace that stops is a tip or a mapping limit. | Derives continuation geometry **from the raster**, not from a precomputed NBMG link table, and gates it on the measured favourability curve while respecting the ≥2 px exclusion. |
| 5 | **H40-B** | Heat-flow anomaly | USGS **DeAngelo et al. (2022)**, DOI [10.5066/P9BZPVUC](https://doi.org/10.5066/P9BZPVUC) — 879 QC A/B/C wells | Convective heat loss along a permeable fault leaves a local positive residual that conduction-only mapping smooths away. | No prior round used a measured heat-flow release; heat was inferred from radiometrics. |
| 6 | **H40-D** | Strain-rate kink + seismic lineament | bands **4, 7, 8, 10, 16** — `geod_2ndinv`, `geod_shearrate`, `geod_dilaterate`, `deq_n100a15`, `ieq_n100a15` | QFFD requires geomorphic evidence; GPS/InSAR strain and microseismicity outline faults that deform without leaving a scarp. | These five bands were **entirely unreachable** before this session (see defect 1). |
| 7 | **H40-E** | Corrected tilt-angle zero-crossing | band **6 `tc`**, corroborated by bands 3 `tmi_hg` and 9 `tmi_vg` | The tilt angle crosses zero directly above a vertical contact — a sub-pixel locator, unlike the analytic-signal magnitude which peaks off-contact. | Repairs a band-index defect *and* replaces magnitude with a zero-crossing locator. |

Fusion is a **geometric mean** (`fuse_play_fairway`): a pixel scores high only if
several independent factors agree, so one noisy channel cannot carry it.

---

## Emission: the budget is a price, not a copied constant

* **Poisson-disk best-first thinning at 2.8 px (280 m).** Justified by the scored
  corpus: 2.8 px beat 1.5 px at the same surface (0.2600 vs 0.2477), consistent
  with redundancy — at 1.5 px several dots share one truth pixel, buying FP
  without buying TP.
* **Hard exclusion within 2 px (200 m) of the catalogue.** Not a masking rule —
  staff say masking is pixel-exact. It is a measured economic fact: deleting
  exactly those 2,545 dots moved 0.2708 → 0.2778.
* **Stop at the break-even price.** Emission proceeds down the ranked propensity
  until the *marginal*, windowed hit rate — transferred to the live scale by the
  anchor's measured instrument/live ratio — falls below `π* = 3.149%`. The
  marginal curve and the chosen cutoff are printed for every candidate, so the
  stopping point is auditable rather than asserted.

---

## Validation protocol (pre-declared before any fold was evaluated)

* **Blocks:** 4 × 6 = 24 rectangular spatial blocks × 2 seeds = **48 folds** per
  instrument; 6 px erosion inside, 10 px dilation collar for component selection.
* **Instrument I1** — hidden whole catalogue components (25% of block catalogue
  length, whole components only, none touching the collar). Emulates the
  organizer's protocol exactly: visible catalogue is `known`, hidden component is
  truth, scoring is on `active & ~known`.
* **Instrument I2** — USGS SGMC off-catalogue faults, **prevalence matched** to
  `|G| = 14,143` by dropping whole components. Without matching, I2 carries
  62,122 truth pixels and inflates every hit rate ~4×.
* **Wald SPRT:** `α = 0.05`, `β = 0.10`, `p0 = 0.5`, `p1 = 0.7`, declared up
  front. Boundaries `B = ln((1−β)/α) = +2.8904`, `A = ln(β/(1−α)) = −2.2513`.
  Deterministic fold order (sort on fold key). No peeking.
* **Power checked first:** `folds_needed()` reports the upper boundary needs **9**
  all-wins. The previous revision's 8-fold design could therefore *never* accept
  H1 (8/8 gives LLR 2.6918 < 2.8904). 48 folds makes the test decidable; pinned
  by `tests/test_h40.py::test_holdout_design_can_reach_a_boundary`.
* **Promotion gate (all four required):** `accept_H1` on the SPRT vs the
  live-scored anchor arm · live-priced hit rate ≥ 1.02 × 6.130% · leakage ≤ 2% ·
  priced score > 0.2778. If nothing clears all four, the **pre-declared fallback**
  emits the best-priced candidate at the anchor's validated budget rather than
  gambling a slot on an unvalidated idea.

---

## Submission name and note

**Name** (unique, not used on any listed GEMSDOE…GEMSDOE39 site):

```
gemsdoe39-h40-pfpt-playfairway-permeability
```

**Note** (paste into the DrivenData submission form; the authoritative string is
`receipt["note"]` in the manifest):

```
GEMSDOE39 H40 play-fairway | H40-PFPT | <N> dots | Poisson 2.8px, catalogue
exclusion 2px | budget set by break-even price pi*=3.15% from live-score
calibrated |G|=14143 | SPRT a=0.05 b=0.10 p0=0.5 p1=0.7
```

**How to upload** — six steps, with the exact validator wording and the reason
each format choice was made: [docs/executive-summary.html](docs/executive-summary.html).

---

## Quick start

```bash
# 0. Dependencies (system pip is PEP-668 blocked in this sandbox -> use a venv)
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# 1. Restore the SHA-256 pinned mirrors (needs `gh` authenticated; no DrivenData creds)
.venv/bin/python scripts/restore_data.py --group core
.venv/bin/python scripts/restore_data.py --group external
.venv/bin/python scripts/restore_data.py --group scored

# 2. Full H40 run: calibrate -> detectors -> holdout -> SPRT -> priced emission -> TIFs
.venv/bin/python scripts/build_h40.py --max-n 70000
#    faster smoke run without the supervised channel:
.venv/bin/python scripts/build_h40.py --fast --max-n 20000

# 3. Validate both twins (12 checks each)
.venv/bin/python scripts/validate_submission.py docs/downloads/<primary>-zeros.tif
.venv/bin/python scripts/validate_submission.py docs/downloads/<primary>-nan.tif

# 4. Tests
.venv/bin/python -m pytest tests/ -q
```

---

## Repository layout

```
data/                            competition rasters (restored, SHA-256 verified; gitignored)
data/external/                   GDR 1391 + GeoDAWN + SGMC derived layers
data/external2/                  USGS Siler 2022 slip tendency; DeAngelo 2022 heat flow
data/scored/                     12 previously scored artifacts (learning only)
data/calib/                      the nested calibration pair (0.2708 / 0.2778) + peers
docs/index.html                  GitHub Pages landing page (download button at top)
docs/executive-summary.html      how to submit + full source register
docs/research/h40-hypotheses.md  the seven-hypothesis register with links
docs/research-brief.md           round-1 brief (superseded by h40-hypotheses.md)
docs/downloads/                  submission GeoTIFF twins + manifest
registry/data_manifest.json      SHA-256 pins for every input raster
registry/live_scores.json        the 16-artifact scored corpus used for calibration
scripts/build_h40.py             end-to-end H40 pipeline (9 stages, fully logged)
scripts/build_pipeline.py        round-1 H39 pipeline (retained; see defect 6)
scripts/restore_data.py          authenticated, hash-verified mirror restore
scripts/validate_submission.py   independent 12-point format validator
scripts/sprt.py                  standalone Wald SPRT CLI
src/gems39/calibrate.py          inverse-DTI calibration against organizer scores
src/gems39/external.py           free official external layers + georeferencing
src/gems39/detectors40.py        H40-A..G channels and play-fairway fusion
src/gems39/holdout40.py          spatially blocked, prevalence-matched instruments
src/gems39/priced_emit.py        break-even-priced best-first Poisson-disk emission
src/gems39/{metric,grid,features,emission,holdout,sprt_select}.py
tests/test_h40.py                every documented claim pinned to a recomputed number
tests/test_pipeline_smoke.py     round-1 smoke test
```

---

## Core values, applied

* **Maximize P(Win).** The score is dominated by a fixed `0.8·|G|` term, so the
  only lever that matters is hits-per-dot. Every design choice in this round is
  aimed at that lever, and the size of the lever is quantified rather than hoped
  for: +18.6% hit rate ⇒ 0.3262.
* **Own the Outcome.** Nine defects in the previous revision were found, fixed and
  pinned by tests (§ Irregularities) rather than reported and left. The pipeline
  refuses to build on a guessed band mapping, refuses to place an external layer
  on an unverified transform, and refuses to promote a candidate that leaks.
* **No hallucinations.** Every external claim carries a link. Every number in the
  documents is either recomputed from bytes on disk or quoted verbatim from a
  fetched page. § "What is not claimed" lists the limits explicitly.
* **SPRT, not peeking.** α, β, p0, p1 declared before any fold is scored; fold
  order deterministic; power checked before the test is run; stopping only at a
  Wald boundary.

---

## Irregularities found and flagged this round

1. **Band-name defect (fixed).** `grid.read_all_bands` hardcoded
   `["mag_anom","rtp","b3",…,"b11","det_elev","iso_grav_anom","tmi",
   "depth_to_base_surf","b16","cond_surf","tc","det_elev_slope"]`. The raster's
   own `band_name` tags put **`tc` at band 6** and `iso_grav_anom_hg` at band 18,
   so every tilt-based detector silently ran on the gravity horizontal gradient —
   and **11 of 19 bands were unreachable** (all three geodetic strain bands, both
   seismic bands, `tmi_hg`, `tmi_vg`, and four gravity-gradient bands). Names are
   now read from the file and asserted; pinned by
   `tests/test_h40.py::test_band_order_is_read_from_the_file_not_hardcoded`.
2. **SPRT audit-record defect (fixed).** `sprt_pairwise` recomputed `wins`/`total`
   by re-iterating the input *after* the accumulation loop; a generator was
   already exhausted, so both were reported as 0 while `n` was correct. It also
   counted folds consumed *after* the stopping boundary.
3. **SPRT had no power (fixed).** 8 folds cannot cross a boundary that needs 9
   wins. Now 48.
4. **Stale timestamp in README (fixed).** It named
   `gemsdoe39-h39-b-20261005T012128Z-zeros.tif`; the file that exists is
   `…T020000Z-zeros.tif`.
5. **Malformed SHA-256 in `docs/index.html` (fixed).** The NaN variant was printed
   as `c985646d4a1bb6de975b65fdc07904faaa55d4739eb176942c537d63f9` — 56 hex
   characters. True digest:
   `c985646d4a1bb6de975290b75b65fdc07904faaa55d4739eb176942c537d63f9`.
6. **Manifest not reproducible from committed code.**
   `docs/downloads/gemsdoe39-h39-b-20261005T020000Z-manifest.json` reports
   `catalogue_hidden_folds: {}` for every candidate while claiming
   `sprt_decisions` with `n=8, wins=8`. `build_pipeline.py` derives `wins` by
   iterating that same dict, so it would have produced `n=0`. The committed run
   could not have produced that manifest.
7. **Holdout proxy off by ~60× (addressed).** The previous manifest's holdout DTI
   (0.0022–0.0051) sits two orders of magnitude below the organizer's 0.24–0.28
   for the same artifacts. Ranking on that proxy is close to ranking on noise;
   the inverse-DTI calibration replaces it.
8. **Emitter silently under-filled its budget (addressed).** `docs/index.html`
   reports 35,805 emitted pixels against `--budget 44090`. The H40 emitter reports
   pool size, chosen cutoff and the marginal hit rate at the cutoff.
9. **`labels.tif` and `existing_faults.tif` are byte-identical** (both 425,830 B,
   both SHA-256 `7ba308ccdc44…`). Flagged, not silently deduplicated.
10. **`sample_submission.tif` is not empty.** It contains 60,988 positive pixels —
    exactly the catalogue count — so the template is the supplied catalogue
    rasterised to float32, not a blank canvas.
11. **OOM in the first holdout draft (fixed).** The run was killed with exit 137
    at stage 5: 24 full-grid block masks plus a stored kernel field per cell
    exceeded the 3.9 GB sandbox, and `read_all_bands` materialised all 19 raw
    bands *and* all 19 cleaned copies at once. Both are now streamed/cropped.
12. **Network reachability.** This sandbox reaches only `github.com`,
    `api.github.com` and `pypi.org`. Verified 2026-10-05: HTTP `000` from
    `sciencebase.gov`, `gdr.openei.org`, `mrdata.usgs.gov`, `web2.nbmg.unr.edu`,
    `osti.gov`, `pangea.stanford.edu`, `community.drivendata.org`. External
    payloads therefore arrive through the owner's SHA-256-pinned GitHub mirrors:
    **integrity-pinned, NOT organizer-authenticated.** The *text* of the forum
    threads, the ScienceBase purpose statement and the live leaderboard were read
    through the page-fetch service and are quoted verbatim with links.

---

## Limitations blocking a guaranteed win

1. **No hidden labels.** The only organizer information available offline is the
   corpus of previously scored artifacts. Both local instruments are proxies.
2. **Single-anchor transfer.** The instrument→live transfer factor is estimated
   from **one** live-scored artifact. A second independent nested pair would
   tighten it; none is available offline.
3. **Nuisance priors in the calibration.** `l0 = 3.0` is exact geometry (the
   triangular-kernel sum along a straight trace) but `kb = 0.5` is an assumption.
   `|G|` moves only ±3% across the whole grid, which is why the result is usable,
   but it is not a measurement.
4. **Sparse external coverage.** 879 QC'd heat-flow wells over 51,674 km² and
   1,244 geothermometer values are thin; the interpolations are two-scale
   Gaussian-weighted local means, **not** kriging — no variogram is fitted.
5. **Submission rate limit.** 3 submissions per rolling 7-day window
   ([thread 11524](https://community.drivendata.org/t/11524)), so only a few
   priced candidates can be tested live per week.
6. **Phase 1 selection is blind and binding.** Competitors must pick one
   submission for *both* the Initial Prize Round (private set) and the Final Prize
   Round (expert-expanded labels from **all** Phase 1 submissions) without knowing
   private performance. Staff: "your fault predictions have an impact on final
   evaluation even if they are not the most performant in Phase 1"
   ([thread 11527](https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527/7)).
   A geologically defensible submission therefore has value beyond its Phase 1
   number — which is a deliberate design input to this round, not an afterthought.

---

## Source register (manual review links)

| # | Source | Link |
|---|---|---|
| 1 | DrivenData #306 problem description, metric, submission format | https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/ |
| 2 | Live leaderboard (read by hand 2026-10-05) | https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/ |
| 3 | Staff: catalogue pixels masked pixel-exactly, both rounds | https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516 |
| 4 | Staff: definition of "new fault" | https://community.drivendata.org/t/where-do-you-draw-the-line/11536 |
| 5 | Staff: no further test-set details; Phase 2 expert review of all submissions | https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527 |
| 6 | USGS Siler (2022) slip & dilation tendency, Quaternary faults, Great Basin | https://doi.org/10.5066/P9YL58W6 · https://www.sciencebase.gov/catalog/item/6296974dd34ec53d276bb33d |
| 7 | USGS DeAngelo et al. (2022) Great Basin heat flow | https://doi.org/10.5066/P9BZPVUC · https://www.sciencebase.gov/catalog/item/6297d2fad34ec53d276c5b28 |
| 8 | USGS Glen et al. (2022) regional geophysical maps of the Great Basin | https://www.sciencebase.gov/catalog/item/628d4fabd34ef70cdba3c4a4 |
| 9 | USGS Peacock & Bedrosian (2022) electrical conductance maps | https://www.sciencebase.gov/catalog/item/62979746d34ec53d276c113b |
| 10 | DOE GDR 1391 — INGENIOUS Great Basin compilation (wells, springs, geothermometers, Qfaults) | https://gdr.openei.org/submissions/1391 · https://doi.org/10.15121/1881483 |
| 11 | INGENIOUS project page | https://gbcge.org/current-projects/ingenious/ |
| 12 | USGS GeoDAWN airborne magnetic & radiometric surveys | https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and · https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7 |
| 13 | USGS Quaternary Fault and Fold Database | https://www.usgs.gov/programs/earthquake-hazards/science/quaternary-fault-and-fold-database-united-states |
| 14 | NBMG Qfaults_INGENIOUS ArcGIS REST layer 0 | https://web2.nbmg.unr.edu/arcgis/rest/services/Qfaults/Qfaults_INGENIOUS/MapServer/0 |
| 15 | USGS State Geologic Map Compilation | https://mrdata.usgs.gov/geology/state/ |
| 16 | USGS 3DEP 1 m DEM | https://registry.opendata.aws/usgs-lidar/ |
| 17 | Faulds & Hinz (2015) favourable settings of Great Basin geothermal systems | https://www.osti.gov/servlets/purl/1724082 · https://www.osti.gov/servlets/purl/1724109 |
| 18 | Barton, Zoback & Moos (1995) *Geology* 23:913 — critically stressed faults are conductive | https://doi.org/10.1130/0091-7613(1995)023<0913:UFFASA>2.3.CO;2 |
| 19 | Morris, Ferrizzoli & Zoback (1996) *Geology* 24:1107 — slip/dilation tendency | https://doi.org/10.1130/0091-7613(1996)024<1107:FPTSOT>2.3.CO;2 |
| 20 | Biasi & Wesnousky (2016) *BSSA* 106:1110 — steps and gaps in ground ruptures | https://doi.org/10.1785/0120150223 |
| 21 | Hermant, Kiersnowski & Bellanger (2025) deep-learning Quaternary fault mapping, Stanford SGP-TR-229 | https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2025/Hermant.pdf |
| 22 | Wald (1945) sequential probability ratio test | https://doi.org/10.1214/aoms/1177731008 |
| 23 | DrivenData reference solution repository (ships a Tversky *training loss*, not the metric) | https://github.com/drivendataorg/gems-prize-reference-solution |

---

## What is explicitly **not** claimed

* No claim that any H40 channel scores above 0.2778 on the organizer's private
  set. The priced number is model output from a two-equation calibration with
  stated nuisance priors and a single-anchor transfer factor.
* No claim that `|G| = 14,143` exactly; the defensible statement is
  `|G| ∈ [13,987, 14,434]` under the stated priors.
* No claim that the interpolations are kriging estimates.
* No claim that SGMC off-catalogue faults *are* the scored population.
* No claim that the external mirrors are the bytes USGS published — they are the
  bytes the owner pinned, and the pin is what is verified.
