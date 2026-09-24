#!/usr/bin/env python3
"""Copy verbatim DAM rows within metadata start minus 10 min through stop.

Usage: python3 code/crop_dam.py ORIGINAL_BATCH_ROOT [RELEASE_ROOT]
Original files are opened read-only. Dates are recording wall-clock labels;
no timezone or daylight-saving conversion is applied.
"""
import argparse
import csv
import datetime as dt
import hashlib
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('release', nargs='?', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    out = args.release / 'data/raw_dam'
    records = []
    dates = {}
    for batch in sorted(args.source.glob('Batch*_DGRP')):
        loading = batch / f'loadinginfo_{batch.name}.csv'
        if not loading.exists():
            continue
        rows = list(csv.DictReader(loading.open()))
        groups = {}
        for r in rows:
            groups.setdefault(r['file'], set()).add((r['start_datetime'], r['stop_datetime']))
        target = out / batch.name
        target.mkdir(parents=True, exist_ok=True)
        for name in [loading.name, f'summary_{batch.name}.csv', f'removed_list1_{batch.name}.csv', f'removed_list2_{batch.name}.csv']:
            p = batch / name
            if p.exists():
                shutil.copy2(p, target / name)
        for name, windows in sorted(groups.items()):
            assert len(windows) == 1, (batch.name, name, windows)
            start, stop = next(iter(windows))
            lower = (dt.datetime.fromisoformat(start) - dt.timedelta(minutes=10)).isoformat(sep=' ').encode()
            upper = stop.encode()
            src, dst = batch / name, target / name
            source_sha, crop_sha = hashlib.sha256(), hashlib.sha256()
            source_md5 = hashlib.md5()
            n_source = n_crop = n_blank = 0
            first = last = None
            with src.open('rb') as inp, dst.open('wb') as output:
                for line in inp:
                    source_sha.update(line); source_md5.update(line)
                    if not line.strip():
                        n_blank += 1
                        continue
                    fields = line.split(b'\t', 3)
                    if len(fields) != 4:
                        raise ValueError(f'Malformed DAM row in {src}: {n_source + 1}')
                    n_source += 1
                    date = fields[1].strip()
                    if date not in dates:
                        dates[date] = dt.datetime.strptime(date.decode(), '%d %b %y').strftime('%Y-%m-%d').encode()
                    stamp = dates[date] + b' ' + fields[2].strip()
                    if lower <= stamp <= upper:
                        output.write(line); crop_sha.update(line); n_crop += 1
                        first = stamp if first is None else min(first, stamp)
                        last = stamp if last is None else max(last, stamp)
            if n_crop == 0:
                raise ValueError(f'No rows within window: {src}')
            records.append(dict(batch=batch.name, file=name,
                source_relative_path=f'{batch.name}/{name}',
                source_bytes=src.stat().st_size, source_rows=n_source,
                source_blank_lines=n_blank, source_sha256=source_sha.hexdigest(), source_md5=source_md5.hexdigest(),
                analysis_start=start, metadata_stop=stop, crop_start_inclusive=lower.decode(), crop_stop_inclusive=upper.decode(),
                crop_first_observed=first.decode(), crop_last_observed=last.decode(),
                crop_rows=n_crop, crop_bytes=dst.stat().st_size, crop_sha256=crop_sha.hexdigest()))
        print(f'{batch.name}: {len(groups)} monitors cropped', flush=True)
    manifest = args.release / 'data/metadata/raw_crop_manifest.tsv'
    with manifest.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(records)
    print(f'Finished {len(records)} monitor files; {sum(r["crop_bytes"] for r in records):,} bytes', flush=True)


if __name__ == '__main__':
    main()
