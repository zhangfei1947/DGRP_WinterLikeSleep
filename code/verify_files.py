#!/usr/bin/env python3
"""Verify release files after downloading and extracting all assets."""
from pathlib import Path
import csv
import hashlib
import sys

root=Path(__file__).resolve().parents[1]
failed=[]
with (root/'MANIFEST.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
    p=root/r['path']
    if not p.is_file():
        failed.append((r['path'],'missing'));continue
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    if p.stat().st_size!=int(r['bytes']) or h.hexdigest()!=r['sha256']:
        failed.append((r['path'],'checksum or size mismatch'))
for x in failed:print(*x,sep=': ')
print(f'{len(rows)-len(failed)}/{len(rows)} files verified')
sys.exit(bool(failed))
