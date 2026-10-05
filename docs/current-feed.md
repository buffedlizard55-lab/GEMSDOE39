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
| Metric is distance-weighted Tversky: α=0.2, β=0.8, 300 m triangular kernel; submission is one float32 GeoTIFF, EPSG:32611, 100 m, matching bounds, `[0,1]` values, null/NaN outside bounds. | [DrivenData problem page, metric & format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Opened both content chunks on 2026-10-05. Official. The currently featured incumbent uses all-finite zero outside to avoid the reported range rejection and includes a NaN twin; portal acceptance of either local file has not been independently verified. |
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
| DrivenData prohibits automated robot/spider access for monitoring or copying site content. | [DrivenData Terms of Use](https://www.drivendata.org/termsofuse/) | Terms page opened 2026-10-05; prohibited-uses section explicitly names robots, spiders, and automatic devices. Do not build an automated leaderboard scraper. |

## Competition-data geometry observed in the restored mirror

The sample and feature rasters read locally on 2026-10-05 have a common grid: shape 3,730 rows × 3,292 columns, EPSG:32611, affine transform `(100, 0, 243350, 0, -100, 4508550)`. The sample has 5,167,373 finite cells and 7,111,787 NaN cells. The restored feature raster has 19 bands with per-band descriptions/tags; H39Y-01 reads `cond_surf` by metadata. The separate owner-mirrored `geodawn_extensions_u8.tif` has four uint8 bands (`ThK`, `UK`, `UTh`, `TMI_up150`) on the same grid; the frozen H39Y-01 formula uses the first three. These observations are reproducible from the local mirror and are not independently organizer-authenticated.

## Known irregularities to retain for review

1. The original project brief says 0.3195 was the current maximum. The 2026-10-05 official page instead shows 0.3262 at rank 1; refresh the dated snapshot when future work is done.
2. The user-reported 0.2778 H33-2-B2 file-to-score mapping conflicts with the owner's GEMSDOE32 page marking that file “UNSCORED.” The 0.2778 board row is real in the rendered snapshot; exact file identity remains unknown.
3. Earlier GEMSDOE39 README/HTML described the all-finite zero-outside incumbent as the only safe encoding. The official format page says null/NaN outside, while the portal has previously returned `Predicted values must be in range [0, 1]`. The repository retains both validated local encodings; actual organizer acceptance of these specific files is not confirmed.
4. Earlier README content pointed to a nonexistent timestamp (`T012128Z`) while the HTML used another TIFF timestamp. Current README and site should link only to the current candidate and its manifest.
5. Earlier documentation said `docs.nlr.gov` was a typo and `docs.nrel.gov` was correct. The supplied PDF resolves from the official `nlr.gov` domain (National Laboratory of the Rockies); the old typo warning was incorrect and has been removed.
6. The exact H39X-01 detector formula was entered into the hypothesis register only after its holdout scores were observed. Fixed SPRT settings do not retroactively preregister the candidate. Its holdout result is exploratory, not confirmatory; the observed SPRT was `continue`, so it is not approved for a weekly slot.
7. Duplicate comparison covered 17 local/retrieved TIFFs only. It found no exact pixelwise match, but does not establish uniqueness across all sibling-site/public artifacts.
8. H39Y-01 was preregistered and evaluated on 2026-10-05; it lost all first five eligible tile comparisons and crossed the SPRT lower boundary (`accept_H0`). No H39Y submission artifact was generated; see the [full report](reports/h39y01-validation-20261005.md).
9. The evaluator's post-stop report assembly initially failed on a duplicate JSON key after it correctly stopped at the H0 boundary. The final report was reconstructed from the committed pre-score lock and saved per-tile outcomes; no fitting, DTI scoring, or candidate emission was rerun. The recovery is detailed in the report.
