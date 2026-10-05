#!/usr/bin/env python3
import struct, os
# One-pixel float32 TIFF, value is bounded and intentionally only a smoke test.
# IFD: width,height,bits,compression,photometric,strip offset, samples/px, rows/strip, strip byte count,sample format.
end=b'II'; entries=[(256,3,1,struct.pack('<H',1)+b'\0\0'),(257,3,1,struct.pack('<H',1)+b'\0\0'),(258,3,1,struct.pack('<H',32)+b'\0\0'),(259,3,1,struct.pack('<H',1)+b'\0\0'),(262,3,1,struct.pack('<H',1)+b'\0\0'),(273,4,1,struct.pack('<I',134)),(277,3,1,struct.pack('<H',1)+b'\0\0'),(278,4,1,struct.pack('<I',1)),(279,4,1,struct.pack('<I',4)),(339,3,1,struct.pack('<H',3)+b'\0\0')]
b=bytearray(end+struct.pack('<HI',42,8)); b+=struct.pack('<H',len(entries));
for tag,typ,n,v in entries:b+=struct.pack('<HHI4s',tag,typ,n,v)
b+=struct.pack('<I',0)+struct.pack('<f',0.731) ; open('artifacts/demo_unique_submission.tif','wb').write(b)
