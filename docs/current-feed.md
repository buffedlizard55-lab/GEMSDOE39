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

## Source register

| Claim used in the repository | Source | Verification / evidence class |
|---|---|---|
| Task is to predict fault structures relevant to geothermal resources; training features include DEM, gravity, magnetic, conductivity, geodetic strain, and earthquake-derived layers. | [DrivenData problem page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Opened on 2026-10-05. Official competition description. Exact 19 band names are also read from the locally restored GeoTIFF metadata and are therefore mirror-level evidence, not organizer-authenticated metadata. |
| Metric is distance-weighted Tversky: α=0.2, β=0.8, 300 m triangular kernel; submission is one float32 GeoTIFF, EPSG:32611, 100 m, matching bounds, `[0,1]` values, null/NaN outside bounds. | [DrivenData problem page, metric & format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Opened both content chunks on 2026-10-05. Official. This is why NaN—not zero—is written outside the footprint. |
| Data download page is login-gated in this environment. | [DrivenData data page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) | Fetch resolves to the DrivenData login form. No competition username/password was requested or stored. |
| Core feature, label, and sample rasters used by this run. | [Local manifest](../registry/data_manifest.json) and [owner mirror](https://github.com/buffedlizard55-lab/GEMSDOE) | Restored through authenticated GitHub Contents API and matched manifest SHA-256 values. **Integrity checked, organizer provenance not authenticated.** Manifest notes labels and `existing_faults.tif` are byte-identical and the sample finite mask is the output footprint. |
| Current top public score is 0.3262; 0.3195 is third on check date. | [Official leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) | Rendered table checked on 2026-10-05. Time-sensitive snapshot, not a private score. |
| Submission limit is up to three scoring submissions per week; only one final submission is selected for prize rounds. | [Official GEMS rules PDF](https://www.nlr.gov/docs/fy26osti/96647.pdf), sections 3.2–3.5 | PDF parsed on 2026-10-05. This repository did not submit a file. |
| Generative AI is allowed but use must be indicated in the narrative; competitor remains responsible for accuracy/authorship. | [Official GEMS rules PDF](https://www.nlr.gov/docs/fy26osti/96647.pdf), section 3.2 | PDF parsed on 2026-10-05. Include disclosure in any prize package. |
| Known USGS/INGENIOUS fault pixels are masked from scoring. | [DrivenData staff reply](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2) | Staff reply read on 2026-10-05. This supports masking exact known cells, not an unverified 200 m exclusion halo. |
| Public and private leaderboard DTI pool pixels then compute one Tversky index. | [DrivenData staff reply](https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550/2) | Staff reply read on 2026-10-05. Our block-specific DTI is a diagnostic; the pooled four-block DTI is reported separately. |
| GeoDAWN includes airborne magnetic/radiometric data, coordinated with LiDAR, and supports geologic/geophysical mapping. | [USGS GeoDAWN release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) · [DOI](https://doi.org/10.5066/P93LGLVQ) | Official USGS page and DOI landing page reachable 2026-10-05; USGS page identifies CC0 1.0 rights. This source was not needed to compute H39X-01. |
| 3DEP provides public elevation products and identifies USGS download paths. | [USGS LiDAR FAQ](https://www.usgs.gov/faqs/what-lidar-data-and-where-can-i-download-it) | Official USGS page opened 2026-10-05. H39X-01 does not consume additional DEM tiles. |
| Great Basin geothermal systems are associated with structural discontinuities; fault stepovers and terminations can localize stress/strain. | [USGS publication: Structural discontinuities and hydrothermal systems](https://www.usgs.gov/publications/structural-discontinuities-and-their-control-hydrothermal-systems-great-basin-usa) | Official USGS publication page opened 2026-10-05. Supports plausibility of structural-strain hypotheses, not proof that the present transform predicts hidden faults. |
| DrivenData prohibits automated robot/spider access for monitoring or copying site content. | [DrivenData Terms of Use](https://www.drivendata.org/termsofuse/) | Terms page opened 2026-10-05; prohibited-uses section explicitly names robots, spiders, and automatic devices. Do not build an automated leaderboard scraper. |

## Competition-data geometry observed in the restored mirror

The sample and feature rasters read locally on 2026-10-05 have a common grid: shape 3,730 rows × 3,292 columns, EPSG:32611, affine transform `(100, 0, 243350, 0, -100, 4508550)`. The sample has 5,167,373 finite cells and 7,111,787 NaN cells. The restored feature raster has 19 bands with per-band `band_name` tags/descriptions. The feature metadata lists `geod_2ndinv`, `geod_shearrate`, and `geod_dilaterate`, which are the three bands used by H39X-01. These observations are reproducible from the local mirror and should not be described as independently organizer-authenticated.

## Known irregularities to retain for review

1. The original project brief says 0.3195 was the current maximum. The 2026-10-05 official page instead shows 0.3262 at rank 1; refresh the dated snapshot when future work is done.
2. The user-reported 0.2778 H33-2-B2 file-to-score mapping conflicts with the owner's GEMSDOE32 page marking that file “UNSCORED.” The 0.2778 board row is real in the rendered snapshot; exact file identity remains unknown.
3. Earlier GEMSDOE39 README/HTML recommended a zeros-outside TIFF and asserted it was safer. The official format page specifies null/NaN outside; this repository now writes NaN outside and the validator enforces it.
4. Earlier README content pointed to a nonexistent timestamp (`T012128Z`) while the HTML used another TIFF timestamp. Current README and site should link only to the current candidate and its manifest.
5. Earlier documentation said `docs.nlr.gov` was a typo and `docs.nrel.gov` was correct. The supplied PDF resolves from the official `nlr.gov` domain (National Laboratory of the Rockies); the old typo warning was incorrect and has been removed.
6. The exact H39X-01 detector formula was entered into the hypothesis register only after its holdout scores were observed. Fixed SPRT settings do not retroactively preregister the candidate. Its holdout result is exploratory, not confirmatory; the observed SPRT was `continue`, so it is not approved for a weekly slot.
7. Duplicate comparison covered 17 local/retrieved TIFFs only. It found no exact pixelwise match, but does not establish uniqueness across all sibling-site/public artifacts.
