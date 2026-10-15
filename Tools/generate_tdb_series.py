"""Extract the full Fairhead series from a caller-supplied ERFA 2.0.1 dtdb.c.
No source download or C build is part of the Swift package.
"""
import argparse,hashlib,re
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('source',type=Path)
args=parser.parse_args()
data=args.source.read_bytes()
text=re.sub(r'/\*.*?\*/','',data.decode(),flags=re.S)
body=re.search(r'fairhd\[787\]\[3\]\s*=\s*\{(.*?)\};',text,re.S).group(1)
rows=re.findall(r'\{([^{}]+)\}',body)
assert len(rows)==787
out=['// Copyright (C) 2013-2023, NumFOCUS Foundation.',
 '// Full Fairhead–Bretagnon series adapted from ERFA 2.0.1 (NumFOCUS).',
 '// BSD terms and SOFA heritage: Sources/ERFASwift/LICENSE.txt.',
 '// Source SHA256: '+hashlib.sha256(data).hexdigest(),
 'internal enum FairheadTDBSeries {']
chunks=[rows[i:i+32] for i in range(0,len(rows),32)]
out.append('  static let terms: [[Double]] = ['+','.join('chunk'+str(i) for i in range(len(chunks)))+'].flatMap { $0 }')
for i,chunk in enumerate(chunks):
 out.append('  private static let chunk'+str(i)+': [[Double]] = [')
 out.extend('    ['+row.strip()+'],' for row in chunk)
 out.append('  ]')
out.append('}')
Path('Sources/ERFASwift/FairheadTDBSeries.swift').write_text('\n'.join(out)+'\n')
