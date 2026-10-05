# Research brief and decision register

**Date:** 2026-10-05. This is a hypothesis register, not evidence of leaderboard performance. The competition score list in the user brief is not independently verified in this checkout.

## Ranked candidates (pre-registration)

1. **Multi-scale DEM curvature residual (medium cost).** Layers: DEM bands; signature: Laplacian-of-Gaussian at several pixel scales, then remove a buffer around `existing_faults`. Target: short-wavelength breaks in slope that are not already catalogued. Distinct from raw ridge/topographic scoring in the prior project records. No external data required.
2. **Positive/negative openness asymmetry (low cost).** DEM-derived sky-view/openess difference across opposing azimuths, excluding catalogue buffers. Target: asymmetric scarps whose line is absent from a catalogue. Distinct transform from elevation/ridge features. No external data required.
3. **Directional relief anisotropy (medium cost).** DEM layer; maximum directional gradient minus median directional gradient over azimuths. Target: coherent linear terrain texture below absolute-relief thresholds. Distinct from isotropic curvature. No external data required.
4. **Thermal/geophysical proxy × curvature (high cost).** Available competition proxy bands crossed with DEM curvature outside catalogue buffers. Target: structural and geothermal coincidence. Viable only if the official feature raster documents such a band; otherwise do not invent one.
5. **Independent USGS 3DEP cross-resolution disagreement (high cost).** Official 3DEP DEM versus competition DEM after CRS/grid audit. Target: acquisition-independent narrow scarps. Source: https://www.usgs.gov/3d-elevation-program . Requires network access and a reproducible tile selection.

## Validation gate

No candidate should consume a weekly submission slot until it beats the locked spatial-block holdout best. Compare one aggregate score per independent spatial block; use `scripts/sprt.py` with alpha and beta declared before reading outcomes. Wald SPRT controls sequential stopping for the specified pair of hypotheses; it does **not** correct selection across a portfolio of candidates. Record all tried candidates and use a separate multiplicity policy.

## Verified-source register

- Competition overview/data: https://www.drivendata.org/competitions/306/competition-doe-gems/ and https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Reference solution: https://github.com/drivendataorg/gems-prize-reference-solution
- USGS 3DEP: https://www.usgs.gov/3d-elevation-program
- NREL report URL supplied by the brief: https://docs.nrel.gov/docs/fy26osti/96647.pdf . The brief also printed `docs.nlr.gov`; that domain spelling is an irregularity and is not treated as verified.

## Current limitations

The repository has no competition rasters and no DrivenData credentials. Consequently no real submission can be generated or claimed scored here. The downloadable demo is intentionally a format smoke test only. A real run needs the official `training_features.tif`, `labels.tif`, and `sample_submission.tif`, plus a documented holdout evaluator. External data must be independently downloaded, licensed, checksum-recorded, reprojected, and tested for leakage.
