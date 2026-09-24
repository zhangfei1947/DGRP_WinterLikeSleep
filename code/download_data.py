#!/usr/bin/env python3
"""Download the fixed v1.0.0 assets, verify SHA-256, and safely extract them."""
from pathlib import Path
import argparse
import hashlib
import json
import tarfile
import urllib.request

root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--asset',action='append',help='Download only a named archive (repeatable); default: all')
args=p.parse_args()
assets=json.loads((root/'ASSETS.json').read_text())
if args.asset:
    known={a['name'] for a in assets}
    if not set(args.asset)<=known:p.error('Unknown asset name')
    assets=[a for a in assets if a['name'] in args.asset]
cache=root/'.downloads';cache.mkdir(exist_ok=True)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
for a in assets:
    dest=cache/a['name']
    if not dest.exists() or sha(dest)!=a['sha256']:
        temp=dest.with_suffix(dest.suffix+'.partial')
        print('Downloading',a['name'],flush=True)
        urllib.request.urlretrieve(a['url'],temp)
        if sha(temp)!=a['sha256']:raise ValueError('Checksum mismatch: '+a['name'])
        temp.replace(dest)
    with tarfile.open(dest,'r:gz') as tar:
        members=tar.getmembers()
        for m in members:
            q=(root/m.name).resolve()
            if not q.is_relative_to(root) or not (m.isfile() or m.isdir()):
                raise ValueError('Unsafe archive member: '+m.name)
        tar.extractall(root,members=members,filter='data')
    print('Verified and extracted',a['name'],flush=True)
