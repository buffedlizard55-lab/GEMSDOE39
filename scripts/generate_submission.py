#!/usr/bin/env python3
"""Generate a bounded, auditable GeoTIFF from competition rasters.
This is deliberately conservative: it refuses missing inputs and never copies a prior submission.
The score is a robust rank ensemble of feature bands, suitable as a baseline only; validate locally.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import rasterio

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference', required=True, help='sample_submission.tif')
    ap.add_argument('--features', required=True, help='training_features.tif or inference feature raster')
    ap.add_argument('--output', required=True)
    ap.add_argument('--bands', default='', help='comma-separated 1-based bands; default all')
    a=ap.parse_args()
    with rasterio.open(a.reference) as ref, rasterio.open(a.features) as src:
        if (ref.width,ref.height,ref.crs,ref.transform)!=(src.width,src.height,src.crs,src.transform):
            raise SystemExit('reference/features grid, CRS, or transform mismatch')
        bands=[int(x) for x in a.bands.split(',') if x] or list(range(1,src.count+1))
        if any(x<1 or x>src.count for x in bands): raise SystemExit('band out of range')
        vals=[]
        for b in bands:
            x=src.read(b, masked=True).astype('float64')
            z=x.filled(np.nan); finite=np.isfinite(z)
            if finite.sum()<2: continue
            lo,hi=np.nanpercentile(z,[2,98]); vals.append(np.clip((z-lo)/(hi-lo if hi>lo else 1),0,1))
        if not vals: raise SystemExit('no usable finite feature bands')
        # Robust consensus; nan cells receive zero, and output is always finite [0,1].
        pred=np.nanmean(np.stack(vals),axis=0); pred=np.nan_to_num(pred,nan=0.0,posinf=1.0,neginf=0.0).astype('float32')
        profile=ref.profile.copy(); profile.update(count=1,dtype='float32',nodata=None,compress='deflate')
        Path(a.output).parent.mkdir(parents=True,exist_ok=True)
        with rasterio.open(a.output,'w',**profile) as dst: dst.write(pred,1)
    digest=hashlib.sha256(Path(a.output).read_bytes()).hexdigest()
    Path(a.output+'.json').write_text(json.dumps({'sha256':digest,'method':'percentile-normalized multi-band consensus','bands':bands,'reference':str(a.reference)},indent=2)+'\n')
    print(a.output, digest)
if __name__=='__main__': main()
