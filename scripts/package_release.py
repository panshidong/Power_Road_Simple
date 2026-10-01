#!/usr/bin/env python3
"""Package prepared tables and small source snapshots; keep the 320 MB official ZIP separate."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import json
ROOT=Path(__file__).resolve().parents[1]
assert json.loads((ROOT/'reports/release_verification.json').read_text())['passed']
paths=[]
for name in ['data/processed','data/raw/road_coordinates','data/scenarios','scripts','vendor/tap-b']:
    paths += [p for p in (ROOT/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.o','.d']]
paths += [p for p in (ROOT/'reports').glob('*') if p.suffix in ['.json','.html','.csv']]
paths += [ROOT/n for n in ['README.md','DATA_SCHEMA.md','Makefile','requirements.txt','sources.lock.json','release_manifest.json','data/raw/road/Austin_sdb_net.txt','data/raw/road/Austin_sdb_trips.txt','data/raw/signals/traffic_signals.json','data/raw/power/Travis150_Electric_Gas.zip']]
paths += [ROOT/name for name in ['runs/tapb_baseline/tapb.log','runs/tapb_baseline/link_flows.csv.gz','runs/tapb_signal_outage/tapb.log','runs/tapb_signal_outage/link_flows.csv.gz']]
dest=ROOT/'Austin_v1_georeferenced.zip'
with ZipFile(dest,'w',compression=ZIP_DEFLATED,compresslevel=5) as z:
    for p in sorted(paths):z.write(p,'Austin/'+str(p.relative_to(ROOT)))
    z.writestr('Austin/PACKAGE_CONTENTS.txt','Prepared tables, scripts, reports, TAP-B and small source snapshots. The original 320,294,926-byte syn-Austin-TDgrid-v03.zip is kept separately at data/raw/power/ in the full workspace. To reconstruct it in a fresh unpacked folder, run: python3 scripts/fetch_sources.py --download-missing. Then create .venv with make setup. Do not interpret this package as a completed road-power coupling or operating-limit validation.\n')
with ZipFile(dest) as z:assert z.testzip() is None
print(dest, dest.stat().st_size, 'bytes')
