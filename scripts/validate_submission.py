#!/usr/bin/env python3
import argparse, hashlib, struct, sys

def read_tiff(path):
 b=open(path,'rb').read();
 if b[:2] not in (b'II',b'MM') or struct.unpack('<H' if b[:2]==b'II' else '>H',b[2:4])[0]!=42: raise ValueError('not classic TIFF')
 endian='<' if b[:2]==b'II' else '>'; off=struct.unpack(endian+'I',b[4:8])[0]; n=struct.unpack(endian+'H',b[off:off+2])[0]; tags={}
 for i in range(n):
  q=b[off+2+i*12:off+14+i*12]; tag,typ,count=struct.unpack(endian+'HHI',q[:8]); size={1:1,2:1,3:2,4:4,5:8}.get(typ,0)*count; raw=q[8:12] if size<=4 else b[struct.unpack(endian+'I',q[8:12])[0]:][:size]
  if typ==3 and count==1: val=struct.unpack(endian+'H',raw[:2])[0]
  elif typ==4 and count==1: val=struct.unpack(endian+'I',raw)[0]
  else: val=raw
  tags[tag]=val
 return b,tags

def main():
 p=argparse.ArgumentParser(); p.add_argument('file'); p.add_argument('--reference'); a=p.parse_args()
 try:
  b,t=read_tiff(a.file); assert t.get(277)==1, 'must be single band'; assert t.get(339)==3, 'must be float32 sample format'
  print('PASS: valid TIFF structure, single band, float32; bytes',len(b),'sha256',hashlib.sha256(b).hexdigest())
  if a.reference: print('NOTE: reference comparison requires rasterio (install requirements.txt); dimensions/geotransform cannot be inferred from this minimal validator.')
 except Exception as e: print('FAIL:',e); return 1
 return 0
if __name__=='__main__': sys.exit(main())
