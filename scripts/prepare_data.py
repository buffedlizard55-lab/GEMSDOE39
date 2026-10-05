#!/usr/bin/env python3
from pathlib import Path
for n in ('training_features.tif','labels.tif','sample_submission.tif'):
 p=Path('data')/n
 if not p.exists(): raise SystemExit(f'MISSING data/{n}; download from the official DrivenData data page')
 print('FOUND',p,p.stat().st_size,'bytes')
print('Placement check passed. Install requirements.txt and add rasterio-based inspection next.')
