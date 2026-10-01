#!/usr/bin/env python3
"""Build auditable Austin tables from frozen source files; never invent road coordinates."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import shlex
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
OUT = ROOT / 'data/processed'


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Fixed gzip timestamp and no embedded filename make table hashes reproducible.
    if path.suffix == '.gz':
        raw = path.open('wb')
        stream = io.TextIOWrapper(gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0), encoding='utf-8', newline='')
    else:
        raw = None
        stream = path.open('w', encoding='utf-8', newline='')
    try:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)
    finally:
        stream.close()
        if raw:
            raw.close()


def metadata(text):
    return {key: value.strip() for key, value in re.findall(r'^\s*<([^>]+)>\s*([^\n]*)', text, re.M) if key != 'END OF METADATA'}


def build_road():
    folder = OUT / 'road'
    source = (RAW / 'road/Austin_sdb_net.txt').read_text()
    meta = metadata(source)
    n, zones, thru = (int(meta[k]) for k in ('NUMBER OF NODES', 'NUMBER OF ZONES', 'FIRST THRU NODE'))
    columns = ['from_node', 'to_node', 'capacity_source', 'length_source', 'free_flow_time_source', 'bpr_alpha', 'bpr_beta', 'speed_source', 'toll_source', 'link_type']
    rows, raw_rows = [], []
    for line in source.split('<END OF METADATA>', 1)[1].splitlines():
        if not line.strip() or line.lstrip().startswith('~'):
            continue
        values = line.split(';', 1)[0].split()
        if len(values) != 10:
            raise ValueError(f'Unexpected TNTP row: {line}')
        raw_rows.append(values)
        row = dict(zip(columns, values))
        u, v = int(values[0]), int(values[1])
        assert 1 <= u <= n and 1 <= v <= n
        assert float(values[2]) > 0 and float(values[4]) >= 0
        row.update(link_id=len(rows) + 1, centroid_connector=int(u < thru or v < thru))
        rows.append(row)
    assert len(rows) == int(meta['NUMBER OF LINKS'])
    pairs = Counter((int(r['from_node']), int(r['to_node'])) for r in rows)
    write_csv(folder / 'links.csv', ['link_id'] + columns + ['centroid_connector'], rows)
    from road_geography import load_coordinates
    coordinates = load_coordinates(ROOT, n)
    write_csv(folder / 'nodes.csv', ['node_id', 'node_type', 'longitude', 'latitude', 'coordinate_status'],
              ({'node_id': i, 'node_type': 'centroid' if i < thru else 'physical', 'longitude': coordinates[i][0], 'latitude': coordinates[i][1], 'coordinate_status': 'source_crosschecked_datum_undeclared'} for i in range(1, n + 1)))
    # Existing road_util assumes eight header lines. Emit exactly eight without changing any link values.
    header = [f'<NUMBER OF ZONES> {zones}', f'<NUMBER OF NODES> {n}', f'<FIRST THRU NODE> {thru}', f'<NUMBER OF LINKS> {len(rows)}', '<END OF METADATA>', '', '~ Austin_sdb; source IDs and numerical values preserved', '~ Tail Head Capacity Length FFT B Power Speed Toll Type ;']
    (folder / 'Austin_net.tntp').write_text('\n'.join(header) + '\n' + '\n'.join('\t'.join(v) + '\t;' for v in raw_rows) + '\n')
    trip_text = (RAW / 'road/Austin_sdb_trips.txt').read_text()
    trip_meta = metadata(trip_text)
    demand = []
    for origin, body in re.findall(r'Origin\s+(\d+)(.*?)(?=Origin\s+\d+|\Z)', trip_text, re.S):
        for dest, volume in re.findall(r'(\d+)\s*:\s*([+\-\d.eE]+)\s*;', body):
            o, d, q = int(origin), int(dest), float(volume)
            assert 1 <= o <= zones and 1 <= d <= zones and q >= 0 and math.isfinite(q)
            demand.append({'origin': o, 'destination': d, 'demand_source': q, 'intrazonal': int(o == d)})
    assert len({(d['origin'], d['destination']) for d in demand}) == len(demand)
    total = math.fsum(d['demand_source'] for d in demand)
    assert abs(total - float(trip_meta['TOTAL OD FLOW'])) < 1e-6
    write_csv(folder / 'od.csv.gz', ['origin', 'destination', 'demand_source', 'intrazonal'], demand)
    (folder / 'Austin_trips.tntp').write_text(trip_text)
    graph = nx.DiGraph()
    graph.add_nodes_from(range(thru, n + 1))
    graph.add_edges_from((u, v) for u, v in pairs if u >= thru and v >= thru)
    components = list(nx.strongly_connected_components(graph))
    # SCC condensation includes one-way appendages; centroids are endpoints only.
    connected = max(components, key=len)
    dag=nx.condensation(graph,components); component=dag.graph['mapping']
    reachable={i:{i}|nx.descendants(dag,i) for i in dag}
    outbound=defaultdict(set);inbound=defaultdict(set)
    for u,v in pairs:
        if u<thru and v>=thru:outbound[u].update(reachable[component[v]])
        if v<thru and u>=thru:inbound[v].add(component[u])
    unreachable=[r for r in demand if r['demand_source']>0 and not r['intrazonal'] and (r['origin'],r['destination']) not in pairs and not outbound[r['origin']].intersection(inbound[r['destination']])]
    assert not unreachable,f'Unreachable OD pairs: {unreachable[:3]}'
    stats = {'source_variant': 'spartalab/tap-b Austin_sdb', 'nodes': n, 'zones': zones, 'physical_nodes': n-thru+1, 'first_thru_node': thru, 'directed_links': len(rows), 'centroid_connectors': sum(r['centroid_connector'] for r in rows), 'duplicate_directed_pairs': {f'{u},{v}': c for (u,v),c in pairs.items() if c>1}, 'physical_scc_count': len(components), 'largest_physical_scc_nodes': len(connected), 'all_positive_interzonal_od_reachable_without_centroid_through_traffic': True, 'od_records': len(demand), 'positive_od_pairs': sum(r['demand_source']>0 for r in demand), 'total_demand_source': total, 'intrazonal_demand_source': sum(r['demand_source'] for r in demand if r['intrazonal']), 'coordinate_status': 'all_source_ids_crosschecked', 'coordinate_reference': 'longitude/latitude decimal degrees; source geodetic datum undeclared', 'coordinate_audit': 'reports/road_coordinate_audit.json', 'units': 'Numerical source units preserved. Time/length/period labels need independent source confirmation; no OD rescaling.'}
    dump(folder / 'metadata.json', stats)
    return stats


def parse_aux(text):
    result = defaultdict(list)
    for match in re.finditer(r'DATA\s*\(\s*(\w+)\s*,\s*\[(.*?)\]\s*\)\s*\{(.*?)\}', text, re.S):
        typ, names, body = match.groups()
        fields = [s.strip() for s in names.split(',')]
        for line in body.splitlines():
            if not line.strip():
                continue
            values = shlex.split(line)
            if len(values) != len(fields):
                raise ValueError(f'AUX {typ}: expected {len(fields)} values, got {len(values)}')
            result[typ].append(dict(zip(fields, values)))
    return result


def build_transmission():
    with ZipFile(RAW / 'power/Travis150_Electric_Gas.zip') as z:
        text = z.read('Travis150_Electric_Data.aux').decode('utf-8-sig')
    folder = OUT / 'power/transmission_travis150_aux'
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'Travis150_Electric_Data.aux').write_text(text)
    tables = parse_aux(text)
    for name, rows in tables.items():
        fields = list(dict.fromkeys(k for row in rows for k in row))
        write_csv(folder / f'{name}.csv', fields, rows)
    subs = {r['SubNum']: r for r in tables['Substation']}
    buses = {r['BusNum']: r for r in tables['Bus']}
    out = []
    for r in tables['Bus']:
        sub = subs.get(r['SubNum'], {})
        out.append({'bus_id': int(r['BusNum']), 'bus_name': r['BusName'], 'base_kv': float(r['BusNomVolt']), 'substation_id': r['SubNum'], 'longitude': r['Longitude'] or sub.get('Longitude', ''), 'latitude': r['Latitude'] or sub.get('Latitude', ''), 'coordinate_source': 'bus' if r['Longitude'] else ('substation' if sub else 'missing'), 'slack': int(r['BusSlack'].strip() == 'YES')})
    write_csv(folder / 'bus_geography.csv', list(out[0]), out)
    for r in tables['Branch']:
        assert r['BusNum'] in buses and r['BusNum:1'] in buses
    for name in ['Gen', 'Load']:
        assert all(r['BusNum'] in buses for r in tables[name])
    stats = {'case_label': 'Travis150_Electric_Gas release; not silently substituted for v03 PWB', 'bus_records': len(buses), 'substation_records': len(subs), 'two_terminal_branch_records': len(tables['Branch']), 'three_winding_transformer_records': len(tables['3WXFormer']), 'generator_records': len(tables['Gen']), 'load_records': len(tables['Load']), 'active_constant_power_load_mw': sum(float(r['LoadSMW']) for r in tables['Load'] if r['LoadStatus']=='Closed'), 'conversion': 'Lossless AUX object tables; branch circuits and 3-winding objects retained separately. No unverified MATPOWER conversion or transmission power-flow claim.'}
    dump(folder / 'metadata.json', stats)
    return stats, out


def commands(text):
    """DSS commands in this frozen archive; retains repeated winding properties verbatim."""
    current = ''
    for line in text.splitlines():
        line = line.split('!', 1)[0].strip()
        if not line:
            continue
        if line.startswith('~'):
            current += ' ' + line[1:].strip()
        else:
            if current:
                yield current
            current = line
    if current:
        yield current


def properties(cmd):
    return {k.lower(): v.strip('"\'') for k,v in re.findall(r'([\w%]+)\s*=\s*("[^"]*"|\[[^]]*\]|\([^)]*\)|\S+)', cmd)}


def unit_xyz(lon_lat):
    a = np.radians(np.asarray(lon_lat, dtype=float))
    return np.column_stack([np.cos(a[:,1])*np.cos(a[:,0]), np.cos(a[:,1])*np.sin(a[:,0]), np.sin(a[:,1])])


def build_distribution():
    folder = OUT / 'power/distribution_v03'
    folder.mkdir(parents=True, exist_ok=True)
    archive = ZipFile(RAW / 'power/syn-Austin-TDgrid-v03.zip')
    names = set(archive.namelist())
    masters = sorted(n for n in names if '--' in n and n.endswith('/Master.dss'))
    subrows = list(csv.DictReader(io.StringIO(archive.read('syn-austin-D_only-v03/all_substations.csv').decode())))
    write_csv(folder / 'substations.csv', list(subrows[0]), subrows)
    subs = {r['Name']: r for r in subrows}
    summaries, load_points, load_keys = [], [], []
    streams = {}
    fields = {
        'loads': ['feeder_id','substation_id','load_id','bus_id','bus_terminal','phases','connection','kv','kw','kvar','longitude','latitude'],
        'buses': ['feeder_id','bus_id','longitude','latitude'],
        'lines': ['feeder_id','line_id','bus1','bus2','bus1_terminal','bus2_terminal','phases','length','length_units','linecode','switch','enabled'],
        'transformers': ['feeder_id','transformer_id','winding_buses_json','raw_dss_command'],
        'linecodes': ['feeder_id','linecode_id','raw_dss_command'],
        'controls': ['feeder_id','source_file','raw_dss_command'],
    }
    for name, cols in fields.items():
        raw = (folder / (name+'.csv.gz')).open('wb')
        stream = io.TextIOWrapper(gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0),encoding='utf-8',newline='')
        w=csv.DictWriter(stream,fieldnames=cols);w.writeheader();streams[name]=(raw,stream,w)
    def write(name, row):streams[name][2].writerow(row)
    try:
        for i, master in enumerate(masters):
            prefix = master.rsplit('/',1)[0]
            feeder = prefix.rsplit('/',1)[-1]
            sub = feeder.split('--')[0]
            assert sub in subs
            coords = {}
            for line in archive.read(prefix+'/Long_lat_buscoords.txt').decode().splitlines():
                bus, lon, lat = line.split()
                coord=(float(lon),float(lat));assert -99<coord[0]<-96 and 29<coord[1]<32
                assert bus.lower() not in coords
                coords[bus.lower()] = coord
                write('buses',dict(feeder_id=feeder,bus_id=bus.lower(),longitude=lon,latitude=lat))
            kw=kvar=0.0; count=0; graph=nx.Graph(); graph.add_nodes_from(coords)
            for cmd in commands(archive.read(prefix+'/Loads.dss').decode()):
                assert cmd.lower().startswith('new load.'), cmd
                p=properties(cmd);bus=p['bus1'].split('.')[0].lower();lon,lat=coords[bus]
                load_id=cmd.split()[1].split('.',1)[1].lower()
                write('loads',dict(feeder_id=feeder,substation_id=sub,load_id=load_id,bus_id=bus,bus_terminal=p['bus1'],phases=p['phases'],connection=p['conn'],kv=p['kv'],kw=p['kw'],kvar=p['kvar'],longitude=lon,latitude=lat))
                kw+=float(p['kw']);kvar+=float(p['kvar']);count+=1
                load_keys.append((feeder,sub,load_id,bus,lon,lat));load_points.append((lon,lat))
            line_count=0
            for cmd in commands(archive.read(prefix+'/Lines.dss').decode()):
                if cmd.lower().startswith('new fuse.'):
                    write('controls',dict(feeder_id=feeder,source_file='Lines.dss',raw_dss_command=cmd))
                    continue
                assert cmd.lower().startswith('new line.'),cmd
                p=properties(cmd);u=p['bus1'].split('.')[0].lower();v=p['bus2'].split('.')[0].lower()
                write('lines',dict(feeder_id=feeder,line_id=cmd.split()[1].split('.',1)[1].lower(),bus1=u,bus2=v,bus1_terminal=p['bus1'],bus2_terminal=p['bus2'],phases=p['phases'],length=p['length'],length_units=p['units'],linecode=p['linecode'],switch=p.get('switch','n'),enabled=p.get('enabled','y')))
                if p.get('enabled','y').lower() in ('y','yes','true'):graph.add_edge(u,v)
                line_count+=1
            transformer_count=0
            for cmd in commands(archive.read(prefix+'/Transformers.dss').decode()):
                assert cmd.lower().startswith('new transformer.'),cmd
                terminals=re.findall(r'\bbus\s*=\s*(\S+)',cmd,re.I)
                if not terminals:
                    array=properties(cmd).get('buses','').strip('()[]')
                    terminals=[s for s in re.split(r'[,\s]+',array) if s]
                assert len(terminals)>=2,cmd
                buses=[s.split('.')[0].lower() for s in terminals]
                graph.add_edges_from((buses[0], b) for b in buses[1:])
                write('transformers',dict(feeder_id=feeder,transformer_id=cmd.split()[1].split('.',1)[1].lower(),winding_buses_json=json.dumps(terminals),raw_dss_command=cmd))
                transformer_count+=1
            for cmd in commands(archive.read(prefix+'/LineCodes.dss').decode()):
                assert cmd.lower().startswith('new linecode.'),cmd
                write('linecodes',dict(feeder_id=feeder,linecode_id=cmd.split()[1].split('.',1)[1].lower(),raw_dss_command=cmd))
            for filename in ['Capacitors.dss','Regulators.dss']:
                if prefix+'/'+filename in names:
                    for cmd in commands(archive.read(prefix+'/'+filename).decode()):
                        write('controls',dict(feeder_id=feeder,source_file=filename,raw_dss_command=cmd))
            mc=list(commands(archive.read(master).decode()))
            circuit=next(c for c in mc if c.lower().startswith('new circuit.'))
            source_bus=properties(circuit)['bus1'].split('.')[0].lower()
            reachable=nx.node_connected_component(graph,source_bus)
            current_loads=load_keys[-count:] if count else []
            assert all(r[3] in reachable for r in current_loads),f'Disconnected loads in {feeder}'
            summaries.append(dict(feeder_id=feeder,substation_id=sub,region=master.split('/')[1],source_bus=source_bus,source_kv=properties(circuit)['basekv'],bus_coordinate_count=len(coords),load_count=count,load_kw=kw,load_kvar=kvar,line_count=line_count,transformer_count=transformer_count,longitude=subs[sub]['Longitude'],latitude=subs[sub]['Latitude'],master_in_zip=master,all_loads_topologically_connected=1))
            if (i+1)%64==0: print(f'Parsed {i+1}/{len(masters)} feeders',flush=True)
    finally:
        for raw,stream,_ in streams.values():stream.close();raw.close()
        archive.close()
    write_csv(folder/'feeders.csv',list(summaries[0]),summaries)
    stats={'feeder_count':len(summaries),'substation_count':len(subs),'load_count':len(load_keys),'load_kw':sum(r['load_kw'] for r in summaries),'load_kvar':sum(r['load_kvar'] for r in summaries),'feeder_bus_coordinate_records':sum(r['bus_coordinate_count'] for r in summaries),'line_records':sum(r['line_count'] for r in summaries),'transformer_records':sum(r['transformer_count'] for r in summaries),'all_feeder_loads_topologically_connected':True,'scope':'Leaf feeders only; upstream substation circuits and interconnections remain in original ZIP. Bus coordinates count buses, not phase nodes. Transformer winding records and line impedance matrices retained verbatim.'}
    dump(folder/'metadata.json',stats)
    return stats,load_keys,load_points,summaries,subrows


def build_signals(load_keys,load_points,feeders,subrows,bus_geo):
    rows=json.loads((RAW/'signals/traffic_signals.json').read_text())
    assert len({r['signal_id'] for r in rows})==len(rows)
    inventory=[]
    for r in rows:
        xy=r.get('location',{}).get('coordinates',[])
        active=(r.get('signal_type')=='TRAFFIC' and r.get('signal_status')=='TURNED_ON' and r.get('control')=='PRIMARY')
        inventory.append(dict(signal_id=r['signal_id'],name=r.get('location_name',''),signal_type=r.get('signal_type',''),status=r.get('signal_status',''),control=r.get('control',''),longitude=xy[0] if len(xy)==2 else '',latitude=xy[1] if len(xy)==2 else '',turn_on_date=r.get('turn_on_date',''),modified_date=r.get('modified_date',''),use_in_v0=int(active and len(xy)==2),road_node_id='',road_match_status='pending_austin_sdb_coordinates',cycle_seconds='',green_ratio='',timing_status='not_in_inventory'))
    write_csv(OUT/'signals/inventory.csv',list(inventory[0]),inventory)
    active=[r for r in inventory if r['use_in_v0']]
    tree=cKDTree(unit_xyz(load_points))
    dist,idx=tree.query(unit_xyz([(r['longitude'],r['latitude']) for r in active]),k=3)
    distances=2*6371008.8*np.arcsin(np.minimum(dist/2,1))
    matches=[]
    for row,ds,indices in zip(active,distances,idx):
        for rank,(d,i) in enumerate(zip(ds,indices),1):
            feeder,sub,load,bus,lon,lat=load_keys[int(i)]
            matches.append(dict(signal_id=row['signal_id'],candidate_rank=rank,feeder_id=feeder,substation_id=sub,load_id=load,distribution_bus_id=bus,distance_m=round(float(d),3),candidate_within_200m=int(d<=200),selected_provisional=int(rank==1 and d<=200),mapping_kind='synthetic_nearest_load_assumption',real_utility_connection_verified=0,road_node_id='',road_match_status='pending_coordinates'))
    write_csv(OUT/'coupling/signal_power_candidates.csv',list(matches[0]),matches)
    # Cross-release spatial crosswalk is audited, never assumed from matching ordinal IDs.
    geo=[r for r in bus_geo if r['longitude']!='' and r['latitude']!='']
    t=cKDTree(unit_xyz([(r['longitude'],r['latitude']) for r in geo]))
    cross=[]
    for sub in subrows:
        point=unit_xyz([(sub['Longitude'],sub['Latitude'])])[0]
        candidates=t.query_ball_point(point,2*math.sin(10/(2*6371008.8)))
        ids=sorted(r['bus_id'] for i in candidates if abs((r:=geo[i])['base_kv']-float(sub['kV']))<1e-5)
        cross.append(dict(distribution_substation_id=sub['Name'],connection_node=sub['Connection Node'],candidate_transmission_bus_ids=json.dumps(ids),match_status='unique_location_voltage_candidate' if len(ids)==1 else ('ambiguous' if ids else 'unmatched'),cross_release_verified=0))
    write_csv(OUT/'coupling/distribution_transmission_candidates.csv',list(cross[0]),cross)
    nearest=distances[:,0]
    stats={'inventory_records':len(rows),'active_primary_traffic_signals_with_coordinates':len(active),'provisional_power_matches_within_200m':int(sum(nearest<=200)),'unmatched_power_over_200m':int(sum(nearest>200)),'nearest_power_distance_m_quantiles':dict(zip(['min','p50','p90','p95','max'],map(float,np.quantile(nearest,[0,.5,.9,.95,1])))),'accepted_road_matches':0,'cross_release_substation_match_counts':dict(Counter(r['match_status'] for r in cross)),'model_assumptions':{'signal_supply':'Nearest synthetic customer load within 200 m; three candidates retained; not measured utility service. Threshold requires sensitivity analysis.','signal_road_mapping':'Blocked by missing Austin_sdb road-node geography.','signal_timing':'Not supplied; left empty. Existing aggregate capacity-loss model does not require phase timing.','signal_outage_capacity_factor':0.5,'capacity_factor_status':'Scenario assumption inherited from small-case design; not calibrated Austin evidence.','date_alignment':'Current signal snapshot and historical network are not contemporaneous.'}}
    dump(OUT/'signals/metadata.json',stats)
    return stats


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    road=build_road();print('Road tables done',flush=True)
    power,buses=build_transmission();print('Transmission AUX tables done',flush=True)
    distribution,keys,points,feeders,subs=build_distribution()
    signals=build_signals(keys,points,feeders,subs,buses)
    from road_geography import build_road_coupling
    geography=build_road_coupling(ROOT)
    signals=json.loads((OUT/'signals/metadata.json').read_text())
    report={'release':'austin-v1-georeferenced','geography':geography,'road':road,'transmission_aux':power,'distribution':distribution,'signals':signals,'readiness':{'tapb_road_input':True,'distribution_opendss_source':True,'full_geographic_road_power_coupling':False,'full_abc_restoration_experiments':False},'remaining_inputs':['Coordinate geodetic datum confirmation and actual road polylines if needed','Reviewed road-node-to-signal assignment and repair access mapping','A chosen transmission/distribution model and validated restoration adapter','Critical destinations, depots, equity regions and repair scenario assumptions']}
    dump(ROOT/'reports/dataset_validation.json',report)
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob('*')) if p.is_file()}
    dump(ROOT/'reports/processed_sha256.json',hashes)
    print(json.dumps(report,indent=2,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
