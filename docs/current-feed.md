# Verified source / competition feed

**Last checked:** 2026-10-05 (UTC)

**Purpose:** a timestamped evidence snapshot with direct links for manual review. It is not an automated DrivenData scraper.

## Current public leaderboard snapshot

The public page was opened and its rendered table read on 2026-10-05:

| Rank | Participant shown by official page | Best public DW-Tversky |
|---:|---|---:|
| 1 | nchuzhoy | 0.3262 |
| 2 | kinghorton42 | 0.3222 |
| 3 | DARD | 0.3195 |
| 4 | xiaofanhu | 0.3060 |
| 5 | alexoktaba | 0.3042 |
| 13 | extradr19 | 0.2778 |

Source: [official DrivenData leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/). These are page values from the check date, not a guarantee that ranks remain the same. The page's rendered table does not identify the submitted TIFF filename for a participant.

### Historical 0.2778 attribution irregularity

- The supplied project history attributes 0.2778 to `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros`.
- The [GEMSDOE32 project page](https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html), which is an owner repository rather than an organizer source, labels the H33-2-B2 file “UNSCORED.”
- The official page shows a 0.2778 row at rank 13 but no filename-to-row crosswalk.
- Conclusion: 0.2778 is present on the official leaderboard snapshot; attribution to that exact H33 artifact is **not verified**. Do not report it as a verified file-level score.

The user-supplied claim that 0.3195 was the current maximum is stale for this checked snapshot: 0.3195 is rank 3 and 0.3262 is rank 1.

## Latest local validation — H39Y-01 (not a leaderboard result)

The single preregistered GeoDAWN radiometric-ratio / `cond_surf` candidate was evaluated once against the checked-in H39-A-model incumbent. The frozen 8×8 plan produced 11 eligible cells. The SPRT stopped after the first five eligible cells: **0 wins, 5 losses**, LLR **−2.554128** crossed the lower boundary **−2.302585** (`accept_H0`). Pooled DTI over those five public-catalogue proxy tiles was **0.035772** for H39Y-01 versus **0.168755** for the incumbent. No H39Y TIFF/ZIP was created and no weekly submission slot is recommended.

- [Detailed result and caveats](reports/h39y01-validation-20261005.md) · [machine-readable outcome report](reports/h39y01-validation-20261005.json)
- [Committed pre-score geometry/training-pool lock](reports/h39y01-preflight-20261005.json) · [candidate portfolio and frozen recipe](preregistration-h39y-20261005.md)
- This is not a DrivenData score. The sequential test's nominal error guarantees remain conditional on its Bernoulli-independence / conditional-supermartingale assumption; spatial buffering cannot prove that assumption.

## Source register

| Claim used in the repository | Source | Verification / evidence class |
|---|---|---|
| Task is to predict fault structures relevant to geothermal resources; training features include DEM, gravity, magnetic, conductivity, geodetic strain, and earthquake-derived layers. | [DrivenData problem page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Opened on 2026-10-05. Official competition description. Exact 19 band names are also read from the locally restored GeoTIFF metadata and are therefore mirror-level evidence, not organizer-authenticated metadata. |
| Metric is distance-weighted Tversky: α=0.2, β=0.8, 300 m triangular kernel; submission is one float32 GeoTIFF, EPSG:32611, 100 m, matching bounds, `[0,1]` values, null/NaN outside bounds. | [DrivenData problem page, metric & format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Official. **Verified from the bytes 2026-10-05:** `sample_submission.tif` is float32 with `nodata=nan` and all 7,111,787 outside cells NaN, so NaN-outside is the template's own convention. `training_features.tif` (19 bands) instead carries `nodata=-3.4028234663852886e+38`. **Both** encodings of our submission ship and `validate_submission.py` passes both (15/15 zeros, 14/14 nan); the all-finite file is recommended because it satisfies the range check under every reading. An earlier claim that the template stores `-3.4e38` was wrong. |
| Data download page is login-gated in this environment. | [DrivenData data page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) | Fetch resolves to the DrivenData login form. No competition username/password was requested or stored. |
| Core feature, label, and sample rasters used by this run. | [Local manifest](../registry/data_manifest.json) and [owner mirror](https://github.com/buffedlizard55-lab/GEMSDOE) | Restored through authenticated GitHub Contents API and matched manifest SHA-256 values. **Integrity checked, organizer provenance not authenticated.** Manifest notes labels and `existing_faults.tif` are byte-identical and the sample finite mask is the output footprint. |
| Current top public score is 0.3262; 0.3195 is third on check date. | [Official leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) | Rendered table checked on 2026-10-05. Time-sensitive snapshot, not a private score. |
| Submission limit is up to three scoring submissions per week; only one final submission is selected for prize rounds. | [Official GEMS rules PDF](https://www.nlr.gov/docs/fy26osti/96647.pdf), sections 3.2–3.5 | PDF parsed on 2026-10-05. This repository did not submit a file. |
| Generative AI is allowed but use must be indicated in the narrative; competitor remains responsible for accuracy/authorship. | [Official GEMS rules PDF](https://www.nlr.gov/docs/fy26osti/96647.pdf), section 3.2 | PDF parsed on 2026-10-05. Include disclosure in any prize package. |
| Known USGS/INGENIOUS faults are masked pixel-exactly; nearby predictions are still scored normally, and new-fault truth may occur within 300 m of a known trace. | [DrivenData staff clarification, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4) | DrivenData Staff explicitly confirmed the mask is identical to provided training-label pixels, with no 300 m scoring halo. H39Y-01 used exact known pixels in DTI scoring; its matched emitter alone used the preregistered 2-pixel exclusion. |
| Public and private leaderboard DTI pool pixels then compute one Tversky index. | [DrivenData staff reply](https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550/2) | Staff reply read on 2026-10-05. Our block-specific DTI is a diagnostic; the pooled four-block DTI is reported separately. |
| GeoDAWN includes airborne magnetic/radiometric data and supports geologic/geophysical mapping; USGS notes airborne K/Th/U radiometry maps surface geology and potassium enrichment may accompany hydrothermal alteration. | [USGS GeoDAWN release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) · [DOI](https://doi.org/10.5066/P93LGLVQ) · [USGS FS 2020-3055](https://pubs.usgs.gov/fs/2020/3055/fs20203055.pdf) | Official USGS sources opened 2026-10-05; the release identifies CC0 1.0 rights. H39Y-01 used owner-mirrored, quantized/aligned ratio layers, not newly retrieved raw USGS bytes. These sources support plausibility only, not measured fault-detection performance. |
| 3DEP provides public elevation products and identifies USGS download paths. | [USGS LiDAR FAQ](https://www.usgs.gov/faqs/what-lidar-data-and-where-can-i-download-it) | Official USGS page opened 2026-10-05. H39X-01 does not consume additional DEM tiles. |
| Great Basin geothermal systems are associated with structural discontinuities; fault stepovers and terminations can localize stress/strain. | [USGS publication: Structural discontinuities and hydrothermal systems](https://www.usgs.gov/publications/structural-discontinuities-and-their-control-hydrothermal-systems-great-basin-usa) | Official USGS publication page opened 2026-10-05. Supports plausibility of structural-strain hypotheses, not proof that the present transform predicts hidden faults. |
| **Label provenance.** Labels come from the USGS Quaternary Fault and Fold Database *plus* "a set of newly identified faults labeled by geology experts at the National Laboratory of the Rockies (NLR) and USGS". Total pool $300,000. Three submissions per week. Single GeoTIFF, single raster layer, 100 m. Generative AI allowed but must be indicated in the narrative. Second-round rankings run against "the complete updated test set created by expert review". | [GEMS Prize Official Rules, September 2026](https://www.nlr.gov/docs/fy26osti/96647.pdf), sections 2 and 3.2 | PDF read 2026-10-05. Official. Directly relevant: expert-interpreted labels mean a *geologically legible* emission is worth more in round 2 than an isotropic scatter, and the round-2 test set is partly built from reviewing our own predictions. |
| DrivenData prohibits automated robot/spider access for monitoring or copying site content. | [DrivenData Terms of Use](https://www.drivendata.org/termsofuse/) | Terms page opened 2026-10-05; prohibited-uses section explicitly names robots, spiders, and automatic devices. Do not build an automated leaderboard scraper. |

## New measurements added 2026-10-05 (session 2), all recomputed from local bytes

| Measurement | Value | Where |
|---|---|---|
| Catalogue pixels within 300 m of any SGMC fault | 34.4 % (footprint base 10.9 %, enrichment 3.16×) | `scripts/run_gems40.py` stage 4 |
| Off-catalogue fault relative density, 200–250 m / 300–400 m / >10 km | 2.43× / 1.58× / 0.46× | `src/gems39/h40.rel_density_profile` |
| ρ(SGMC adjacency, live score) over 12 historical artifacts | +0.698, p = 0.012 | `registry/h40_report_h40e-30k.json` |
| Calibrated instrument fit (n, w) → live score over 12 artifacts | RMSE 0.0364, LOO-RMSE 0.0397, ρ +0.923 / +0.902 LOO | same |
| Blocked out-of-fold AUC of the off-catalogue discriminant | 0.7489 over 11 spatial blocks | same |
| Shipped file: dots, w, predicted score | 30,000, 0.1461, 0.3006 (extrapolation) | same |

These are **local** measurements against a surrogate. None of them is an organiser score.

## Competition-data geometry observed in the restored mirror

The sample and feature rasters read locally on 2026-10-05 have a common grid: shape 3,730 rows × 3,292 columns, EPSG:32611, affine transform `(100, 0, 243350, 0, -100, 4508550)`. The sample has 5,167,373 finite cells and 7,111,787 NaN cells. The restored feature raster has 19 bands with per-band descriptions/tags; H39Y-01 reads `cond_surf` by metadata. The separate owner-mirrored `geodawn_extensions_u8.tif` has four uint8 bands (`ThK`, `UK`, `UTh`, `TMI_up150`) on the same grid; the frozen H39Y-01 formula uses the first three. These observations are reproducible from the local mirror and are not independently organizer-authenticated.

## Known irregularities to retain for review

1. The original project brief says 0.3195 was the current maximum. The 2026-10-05 official page instead shows 0.3262 at rank 1; refresh the dated snapshot when future work is done.
2. The user-reported 0.2778 H33-2-B2 file-to-score mapping conflicts with the owner's GEMSDOE32 page marking that file “UNSCORED.” The 0.2778 board row is real in the rendered snapshot; exact file identity remains unknown.
3. **Resolved 2026-10-05, reversing an earlier statement.** The official format page describes null/NaN outside, and `sample_submission.tif` is indeed NaN outside with `nodata=nan` (verified from the bytes). But `validate_submission.py` does **not** enforce NaN-only: it accepts both encodings. This repository therefore ships **both** files and recommends the all-finite `…-zeros.tif` for upload, because a prior attempt was rejected with `Predicted values must be in range [0, 1]` and an all-finite file satisfies that check under every reading. The server-side validator cannot be inspected locally, so the causal mechanism is **not** established; the two candidate mechanisms are a leaked `-3.4e38` sentinel from `training_features.tif`, or NaN failing an IEEE-754 range comparison.
4. Earlier README content pointed to a nonexistent timestamp (`T012128Z`) while the HTML used another TIFF timestamp. Current README and site should link only to the current candidate and its manifest.
5. Earlier documentation said `docs.nlr.gov` was a typo and `docs.nrel.gov` was correct. **Re-verified 2026-10-05:** `https://www.nlr.gov/docs/fy26osti/96647.pdf` resolves (redirecting to `docs.nlr.gov`) and returns *Geologic Enhanced Mapping System (GEMS) Prize Official Rules, September 2026*. NLR = National Laboratory of the Rockies, the prize administrator. The `nrel.gov` variant is **not** needed; use `nlr.gov`.
6. The exact H39X-01 detector formula was entered into the hypothesis register only after its holdout scores were observed. Fixed SPRT settings do not retroactively preregister the candidate. Its holdout result is exploratory, not confirmatory; the observed SPRT was `continue`, so it is not approved for a weekly slot.
7. Duplicate comparison for the shipped H40-E file covered 15 same-grid rasters (12 organiser-scored historical artifacts + the 3 files this repository previously shipped); the zeros/nan twins of the same run are excluded as vacuous. It found no exact pixelwise match, but does not establish uniqueness across all sibling-site/public artifacts.
8. H39Y-01 was preregistered and evaluated on 2026-10-05; it lost all first five eligible tile comparisons and crossed the SPRT lower boundary (`accept_H0`). No H39Y submission artifact was generated; see the [full report](reports/h39y01-validation-20261005.md).
9. The evaluator's post-stop report assembly initially failed on a duplicate JSON key after it correctly stopped at the H0 boundary. The final report was reconstructed from the committed pre-score lock and saved per-tile outcomes; no fitting, tile predictions, or per-tile DTI scoring was rerun. Pooled DTI was re-aggregated from the saved per-tile TP/FP/FN counts, with no candidate emission. The recovery is detailed in the report.
