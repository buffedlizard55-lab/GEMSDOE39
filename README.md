# GEMSDOE39 — auditable GEMS submission lab

> **Mission:** maximize P(Win) and own the outcome while producing a *unique* geothermal-emission submission. Every candidate must be reproducible, bounded to `[0,1]`, spatially blocked, and evaluated with a pre-declared sequential test before a competition submission.

## Status (2026-10-05)

This checkout is a clean scaffold; it contains no competition rasters and no DrivenData credentials. Consequently it cannot honestly claim a leaderboard score or create a competition-compatible GeoTIFF yet. The public Dropbox URLs in the project brief were not used as authoritative competition data. Put the official files in `data/` (see below), then run the pipeline. The included `artifacts/demo_unique_submission.tif` is deliberately **not** a competition submission: it is a one-pixel smoke-test artifact and is labeled as such.

## Quick start

```bash
python scripts/validate_submission.py artifacts/demo_unique_submission.tif
# With official sample_submission.tif present:
python scripts/generate_submission.py --reference data/sample_submission.tif --output artifacts/unique_submission.tif
python scripts/validate_submission.py artifacts/unique_submission.tif --reference data/sample_submission.tif
```

`generate_submission.py` never copies prediction values from a previous submission. It computes a deterministic, novel bounded score from the supplied feature stack (or uses a clearly marked demo fallback), preserves the reference CRS/shape/geotransform, writes a single-band float32 GeoTIFF, and emits a manifest with SHA-256 and method parameters. Do not upload anything unless validation says `PASS`.

## Executive summary / submission procedure

1. Obtain `training_features.tif`, `labels.tif`, and `sample_submission.tif` from the official DrivenData competition data page after authenticating in your own browser. No credentials are stored here.
2. Place them in `data/`; run `python scripts/prepare_data.py` to inspect dimensions, CRS, transform, bands, finite ranges, and label geometry. Missing or mismatched files fail closed.
3. Implement/evaluate a candidate using the spatial blocks and `scripts/sprt.py`. The SPRT uses pre-declared alpha/beta and stops only at Wald's upper/lower boundary; no visual peeking or “run longer” rule is permitted.
4. Generate the candidate GeoTIFF using the reference grid, validate range `[0,1]`, single band, finite values, CRS, shape and transform, then download `artifacts/unique_submission.tif`.
5. On DrivenData choose **New submission → File to submit**, upload that file, and add a unique note such as `gemsdoe39-curvature-residual-sprt-a05-b10-20261005`. Scores are unknown until the competition evaluates the file; record the returned score in `reports/submission_log.csv`.

## Candidate hypotheses (ranked before validation)

| Rank | Hypothesis / layers and signature | Why it may find catalogue-missing faults | Difference / cost |
|---|---|---|---|
| 1 | DEM-derived **multi-scale Laplacian-of-Gaussian curvature residual**, masked by existing-fault distance | A short-wavelength break in slope can indicate a concealed structural boundary; distance masking avoids rewarding catalogue lines | No external data; distinct from raw DEM/ridge scores; medium |
| 2 | DEM **positive/negative openness asymmetry** plus local relief, with catalogue-distance exclusion | Opposing terrain openness may expose narrow scarps not represented as a mapped line | No external data; distinct transform; low |
| 3 | DEM **directional relief anisotropy**: maximum directional gradient minus median directional gradient | A lineament has orientation coherence even where absolute relief is weak | No external data; distinct directional statistic; medium/high |
| 4 | Geothermal proxy layer × DEM curvature, only outside catalogue buffers | Thermal/structural coincidence is more specific than either alone, but depends on an official proxy band being present | Uses existing competition features only; high |
| 5 | External USGS 3DEP DEM cross-resolution disagreement | Independent acquisition artefacts can identify real narrow scarps missed by one DEM | Requires free official 3DEP access and reprojection audit; high |

These are hypotheses, not results. The top candidate must beat the current holdout best under the locked protocol before any submission slot is spent. The cited competition score history in the brief is user-provided and is not independently verified here.

## Verification and limitations

- Official competition overview: https://www.drivendata.org/competitions/306/competition-doe-gems/
- Official problem description: https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Official data page: https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Official reference solution: https://github.com/drivendataorg/gems-prize-reference-solution
- Official competition PDF link supplied in the brief: https://docs.nrel.gov/docs/fy26osti/96647.pdf (availability should be checked manually; the brief says `nlr.gov`, which is flagged as an irregularity)
- USGS 3DEP landing page: https://www.usgs.gov/3d-elevation-program

No score, data download, external research result, or “highest submission” claim is fabricated. The SPRT is for sequential holdout comparison; it does not correct selection among many candidate hypotheses. Use a pre-registered candidate family and account for portfolio selection separately.

## Research register and full project brief

The ranked hypotheses, source register, validation gate, and limitations are maintained in [docs/research-brief.md](docs/research-brief.md). The original task brief is preserved in the session request; claims are intentionally not promoted to facts unless locally reproducible or linked to an official source.

## Core values

**Maximize P(Win. Own the Outcome.** Every failure is explicit, logged, and actionable; no silent fallback is allowed for a real submission.
