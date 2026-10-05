#!/usr/bin/env python3
"""Fail-closed competition GeoTIFF validator."""
import argparse, hashlib, sys
import numpy as np
import rasterio

def main():
 p=argparse.ArgumentParser(); p.add_argument('file'); p.add_argument('--reference'); a=p.parse_args()
 try:
  with rasterio.open(a.file) as ds:
   if ds.count!=1: raise ValueError('must be single-band')
   if ds.dtypes[0] not in ('float32','float64'): raise ValueError('prediction must be floating point')
   x=ds.read(1,masked=True)
   if x.size==0 or not np.isfinite(x.compressed()).all(): raise ValueError('non-finite values')
   y=x.compressed()
   if y.size and (y.min()<0 or y.max()>1): raise ValueError(f'predicted values outside [0,1]: {y.min()}..{y.max()}')
   if a.reference:
    with rasterio.open(a.reference) as r:
     for k in ('width','height','crs','transform'):
      if getattr(ds,k)!=getattr(r,k): raise ValueError(f'{k} mismatch')
   print('PASS',ds.width,ds.height,ds.crs,'sha256',hashlib.sha256(open(a.file,'rb').read()).hexdigest())
 except Exception as e: print('FAIL:',e); return 1
 return 0
if __name__=='__main__': sys.exit(main())
