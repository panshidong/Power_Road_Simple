#!/usr/bin/env python3
"""Independent consistency checks across delivered tables, frozen input and solver evidence."""
import csv,gzip,hashlib,json,math,subprocess,sys
from collections import Counter,defaultdict
from pathlib import Path
from build_dataset import ROOT,OUT,dump

def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return list(csv.DictReader(f))

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

sources=json.loads((ROOT/'sources.lock.json').read_text())
for r in sources['sources']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
oldhash=json.loads((ROOT/'reports/processed_sha256.json').read_text())
for path,expected in oldhash.items():assert sha(ROOT/path)==expected,path
report=json.loads((ROOT/'reports/dataset_validation.json').read_text())
feeders={r['feeder_id']:r for r in read(OUT/'power/distribution_v03/feeders.csv')}
assert len(feeders)==report['distribution']['feeder_count']==448
quality=read(OUT/'power/distribution_v03/solver_quality.csv')
assert {r['feeder_id'] for r in quality}==set(feeders)
assert all(r['eligibility']=='import_solve_ready' for r in quality)
rawcheck=json.loads((ROOT/'reports/opendss_validation.json').read_text())
assert rawcheck['tested_feeders']==448
failures={r['feeder_id'] for r in rawcheck['results'] if not r['passed']}
assert failures=={'p3uhs2_1247--p3udt31865','p5uhs1_1247--p5udt4629'}
audit=json.loads((ROOT/'reports/power_exception_audit.json').read_text())
assert len(audit['voltagebase_corrections'])==11 and all(r['passed'] and r['max_physical_voltage_change_v']==0 for r in audit['voltagebase_corrections'])
assert {r['region'] for r in audit['original_regional_solves']}=={'P3U','P5U'}
assert all(r['passed'] for r in audit['original_regional_solves'])
assert all(not r['unenergized_loads'] for r in rawcheck['results'] if r['passed'])

signals=read(OUT/'signals/inventory.csv');candidates=read(OUT/'coupling/signal_power_candidates.csv')
assert len(signals)==1345
assert all(r['cycle_seconds']=='' and r['green_ratio']=='' for r in signals)
active={r['signal_id'] for r in signals if r['use_in_v0']=='1'}
assert len(active)==965
assert Counter(r['signal_id'] for r in candidates)==Counter({i:3 for i in active})
selected=[r for r in candidates if r['selected_provisional']=='1']
assert len(selected)==889
assert all(r['candidate_rank']=='1' and float(r['distance_m'])<=200 for r in selected)
signal_by_id={r['signal_id']:r for r in signals}
assert all(r['real_utility_connection_verified']=='0' and r['road_node_id']==signal_by_id[r['signal_id']]['road_node_id'] for r in candidates)
needed={(r['feeder_id'],r['load_id']) for r in candidates}
found={};counts=Counter();kw=defaultdict(float);all_load_ids=set()
with gzip.open(OUT/'power/distribution_v03/loads.csv.gz','rt') as f:
    for r in csv.DictReader(f):
        key=(r['feeder_id'],r['load_id']);assert key not in all_load_ids;all_load_ids.add(key)
        counts[r['feeder_id']]+=1;kw[r['feeder_id']]+=float(r['kw'])
        if key in needed:found[key]=r
assert len(all_load_ids)==307236
for feeder,r in feeders.items():
    assert counts[feeder]==int(r['load_count'])
    assert math.isclose(kw[feeder],float(r['load_kw']),rel_tol=1e-10)
for r in candidates:
    load=found[(r['feeder_id'],r['load_id'])]
    assert load['bus_id']==r['distribution_bus_id'] and load['substation_id']==r['substation_id']
road=read(OUT/'road/links.csv');nodes=read(OUT/'road/nodes.csv')
assert len(road)==18710 and len(nodes)==7466
assert len({r['link_id'] for r in road})==len(road)
with (ROOT/'data/raw/road_coordinates/tap_demand_Austin_sdb_node.txt').open() as f:
    source_nodes={r['Node']:r for r in csv.DictReader(f,delimiter='\t')}
assert set(source_nodes)=={r['node_id'] for r in nodes}
for n in nodes:
    src=source_nodes[n['node_id']]
    assert float(n['longitude'])==float(src['X']) and float(n['latitude'])==float(src['Y'])
    assert n['coordinate_status']=='source_crosschecked_datum_undeclared'
coordinate_audit=json.loads((ROOT/'reports/road_coordinate_audit.json').read_text())
assert coordinate_audit['passed'] and coordinate_audit['secondary_road_fft_differences']==286
road_candidates=read(OUT/'coupling/signal_road_candidates.csv')
assert Counter(r['signal_id'] for r in road_candidates)==Counter({i:3 for i in active})
chosen={r['signal_id']:r for r in road_candidates if r['selected_provisional']=='1'}
assert len(chosen)==428 and len({r['road_node_id'] for r in chosen.values()})==428
nodes_by_id={r['node_id']:r for r in nodes}
neighbours=defaultdict(set); incoming=defaultdict(set)
for link in road:
    if link['centroid_connector']=='0':
        neighbours[link['from_node']].add(link['to_node']);neighbours[link['to_node']].add(link['from_node'])
        incoming[link['to_node']].add(link['link_id'])
for r in road_candidates:
    n=nodes_by_id[r['road_node_id']];sig=signal_by_id[r['signal_id']]
    assert n['node_type']=='physical' and len(neighbours[r['road_node_id']])>=3
    lon1,lat1,lon2,lat2=map(math.radians,[float(n['longitude']),float(n['latitude']),float(sig['longitude']),float(sig['latitude'])])
    h=math.sin((lat1-lat2)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin((lon1-lon2)/2)**2
    d=2*6371008.8*math.asin(math.sqrt(h));assert abs(d-float(r['distance_m']))<.001
    assert r['manual_verified']=='0'
for sid,r in chosen.items():
    assert r['candidate_rank']=='1' and float(r['distance_m'])<=50 and float(r['second_minus_first_distance_m'])>=20
    assert signal_by_id[sid]['road_node_id']==r['road_node_id']
assert all((r['road_node_id']!='')==(r['signal_id'] in chosen) for r in signals)
approaches=read(OUT/'coupling/signal_approach_links.csv')
assert len(approaches)==1427 and len({r['link_id'] for r in approaches})==1427
by_signal=defaultdict(set)
for r in approaches:
    assert r['signal_id'] in chosen and r['to_node']==chosen[r['signal_id']]['road_node_id']
    assert r['manual_verified']=='0';by_signal[r['signal_id']].add(r['link_id'])
assert all(by_signal[sid]==incoming[r['road_node_id']] for sid,r in chosen.items())
joint=read(OUT/'coupling/signal_road_power_provisional.csv')
assert len(joint)==410
assert {r['signal_id'] for r in joint}==set(chosen)&{r['signal_id'] for r in selected}
assert all(r['road_manual_verified']=='0' and r['real_utility_connection_verified']=='0' for r in joint)
access=read(OUT/'coupling/substation_road_access_candidates.csv')
assert len(access)==384 and all(r['actual_access_verified']=='0' and nodes_by_id[r['road_node_id']]['node_type']=='physical' for r in access)

pairs=Counter((r['from_node'],r['to_node']) for r in road);assert sum(v>1 for v in pairs.values())==7
check=json.loads((ROOT/'reports/tapb_validation.json').read_text())
assert check['passed'] and check['relative_gap']<=1e-4 and check['link_flow_rows']==18710
assert check['network_sha256']==sha(OUT/'road/Austin_net.tntp') and check['trips_sha256']==sha(OUT/'road/Austin_trips.tntp')
assert check['max_node_flow_balance_error']<1e-4 and check['parallel_links_preserved_in_row_order']
scenario=json.loads((ROOT/'reports/tapb_signal_outage_validation.json').read_text())
assert scenario['passed'] and scenario['relative_gap']<=1e-4 and scenario['only_selected_capacity_fields_changed']
assert scenario['network_sha256']==sha(ROOT/'data/scenarios/signal_outage_demo/Austin_signal_outage_net.tntp')
assert scenario['trips_sha256']==check['trips_sha256'] and scenario['max_node_flow_balance_error']<1e-4
assert not scenario['real_event_or_calibrated_outage']

vendor=json.loads((ROOT/'vendor/tap-b/SOURCE.json').read_text())
assert sha(ROOT/'vendor/tap-b/bin/tap')==vendor['binary_sha256']
assert check['binary_sha256']==vendor['binary_sha256']
for name,h in vendor['source_sha256'].items():assert sha(ROOT/'vendor/tap-b'/name)==h
# Validate inline JavaScript without requiring network or a browser.
html=(ROOT/'reports/geographic_preview.html').read_text();script=html.split('<script>',1)[1].split('</script>',1)[0]
js=ROOT/'runs/preview-syntax-check.js';js.write_text(script)
subprocess.run(['node','--check',str(js)],check=True)
for path in (ROOT/'scripts').glob('*.py'):compile(path.read_text(),str(path),'exec')
summary={'passed':True,'source_files_verified':len(sources['sources']),'processed_files_verified_against_build_hashes':len(oldhash),'road_links':len(road),'loads':len(all_load_ids),'feeders_import_solve_ready':len(quality),'raw_standalone_feeders_passed':448-len(failures),'regional_dependency_exceptions_verified':len(failures),'voltagebase_metadata_corrections_verified':11,'signal_power_provisional_matches':len(selected),'source_road_nodes_with_coordinates':len(nodes),'provisional_road_signal_matches':len(chosen),'provisional_joint_signal_road_power_matches':len(joint),'actual_road_signal_matches_verified':0,'tapb_signal_outage_scenario_passed':scenario['passed'],'html_javascript_syntax_checked':True,'full_geographic_coupling_ready':False,'operating_limits_and_restoration_validated':False}
dump(ROOT/'reports/release_verification.json',summary)
files=[p for folder in [OUT,ROOT/'data/scenarios',ROOT/'scripts'] for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
files += [ROOT/n for n in ['README.md','DATA_SCHEMA.md','Makefile','requirements.txt','sources.lock.json']]
dump(ROOT/'release_manifest.json',{'release':'austin-v1-georeferenced','files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(files)},'verification':summary})
print(json.dumps(summary,indent=2,ensure_ascii=False))
