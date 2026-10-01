#!/usr/bin/env python3
"""Verify frozen sources, or download missing sources without overwriting existing files."""
import argparse,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

p=argparse.ArgumentParser();p.add_argument('--download-missing',action='store_true');args=p.parse_args()
manifest=json.loads((ROOT/'sources.lock.json').read_text())
for row in manifest['sources']:
    path=ROOT/row['path']
    if not path.exists():
        if not args.download_missing:raise SystemExit(f'Missing {path}; use --download-missing')
        path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.download')
        subprocess.run(['curl','-fLsS','--connect-timeout','15','--max-time','600','--retry','2',row['url'],'-o',str(temp)],check=True)
        actual=sha(temp)
        if actual!=row['sha256']:
            raise SystemExit(f'Source changed: {row["path"]}; kept candidate at {temp}. Expected {row["sha256"]}; got {actual}. Live signal inventory cannot reproduce an old snapshot; keep the archived original.')
        temp.rename(path)
    if sha(path)!=row['sha256']:raise SystemExit(f'Checksum mismatch: {path}')
    print(f'OK {row["path"]}',flush=True)
