#!/usr/bin/env python3
"""Audit the Austin_sdb coordinate crosswalk and build provisional spatial joins."""
from __future__ import annotations
import csv, hashlib, json, math, re
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

RADIUS = 6371008.8
ROAD_LIMIT_M = 50.0
SEPARATION_M = 20.0


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def xyz(points):
    p = np.radians(np.asarray(points, dtype=float))
    return np.column_stack((np.cos(p[:, 1])*np.cos(p[:, 0]), np.cos(p[:, 1])*np.sin(p[:, 0]), np.sin(p[:, 1])))


def metres(chord):
    return 2*RADIUS*np.arcsin(np.minimum(np.asarray(chord)/2, 1))


def quantiles(values):
    return dict(zip(['min','p50','p90','p95','max'], map(float, np.quantile(values,[0,.5,.9,.95,1]))))


def net_rows(path):
    return [tuple(float(x) for x in line.split(';')[0].split())
            for line in path.read_text().split('<END OF METADATA>',1)[1].splitlines()
            if line.strip() and not line.lstrip().startswith('~')]


def od_rows(path):
    rows = [(int(o), int(d), float(q))
            for o, body in re.findall(r'Origin\s+(\d+)(.*?)(?=Origin\s+\d+|\Z)',path.read_text(), re.S)
            for d,q in re.findall(r'(\d+)\s*:\s*([+\-\d.eE]+)\s*;', body)]
    assert len(rows)==len({r[:2] for r in rows})
    return {(o,d):q for o,d,q in rows}


def load_coordinates(root, n):
    raw=root/'data/raw'; folder=raw/'road_coordinates'
    def coords(name):
        with (folder/name).open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        result={int(r['Node']):(float(r['X']),float(r['Y'])) for r in rows}
        assert len(rows)==len(result)==n
        assert set(result)==set(range(1,n+1))
        assert all(-99<lon<-96 and 29<lat<32 for lon,lat in result.values())
        return result
    primary=coords('tap_demand_Austin_sdb_node.txt')
    secondary=coords('dstap_Austin_sdb_node.txt')
    assert primary==secondary, 'Independent repository coordinate tables disagree'
    local=net_rows(raw/'road/Austin_sdb_net.txt')
    other=net_rows(folder/'tap_demand_Austin_sdb_net.txt')
    assert local==other, 'Coordinate source road rows differ from local TAP-B road'
    od=od_rows(raw/'road/Austin_sdb_trips.txt')
    assert od==od_rows(folder/'tap_demand_Austin_sdb_trips.txt')
    dstap=net_rows(folder/'dstap_Austin_sdb_net.txt')
    assert [r[:2] for r in local]==[r[:2] for r in dstap]
    changed=[dict(link_id=i,from_node=int(a[0]),to_node=int(a[1]),local_fft=a[4],dstap_fft=b[4])
             for i,(a,b) in enumerate(zip(local,dstap),1) if a!=b]
    assert all(a[:4]+a[5:]==b[:4]+b[5:] for a,b in zip(local,dstap))
    assert od==od_rows(folder/'dstap_Austin_sdb_trips.txt')
    write(root/'reports/coordinate_source_road_differences.csv',changed)
    report={'passed':True,'node_count':n,'all_ids_covered':True,
            'primary':'spartalab/TAP_Demand@9e7463fbe2911924767b08a7e7429178efdb1e71',
            'crosscheck':'venktesh22/DSTAP@2e2d95234f3c580e26683231832866678f31385b',
            'coordinate_tables_numerically_identical':True,
            'primary_road_all_rows_all_ten_fields_same_order':True,
            'primary_and_secondary_od_equal_local':True,
            'secondary_road_endpoint_order_identical':True,
            'secondary_road_fft_differences':len(changed),
            'secondary_road_parameters_adopted':False,
            'bounds_lon_lat':[list(np.min(list(primary.values()),axis=0)),list(np.max(list(primary.values()),axis=0))],
            'coordinate_interpretation':'X=longitude, Y=latitude in decimal degrees; source does not declare a geodetic datum. WGS84-compatible plotting and sphere distances are an explicit assumption.',
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file()}}
    dump(root/'reports/road_coordinate_audit.json',report)
    return primary


def build_road_coupling(root):
    out=root/'data/processed'
    nodes=read(out/'road/nodes.csv'); links=read(out/'road/links.csv')
    byid={int(n['node_id']):n for n in nodes}
    physical={int(n['node_id']) for n in nodes if n['node_type']=='physical'}
    coords={i:(float(n['longitude']),float(n['latitude'])) for i,n in byid.items()}
    neighbours=defaultdict(set); incoming=defaultdict(list); outgoing=defaultdict(list)
    for link in links:
        u,v=int(link['from_node']),int(link['to_node'])
        if u in physical and v in physical:
            neighbours[u].add(v); neighbours[v].add(u)
            incoming[v].append(link); outgoing[u].append(link)
    # No centroids, dead ends, or two-neighbour geometry nodes are signal candidates.
    eligible=sorted(i for i in physical if len(neighbours[i])>=3 and incoming[i] and outgoing[i])
    tree=cKDTree(xyz([coords[i] for i in eligible]))
    inventory=read(out/'signals/inventory.csv')
    active=[s for s in inventory if s['use_in_v0']=='1']
    dist,ix=tree.query(xyz([(s['longitude'],s['latitude']) for s in active]),k=3)
    ds=metres(dist)
    prelim={s['signal_id']:eligible[ix[k,0]] for k,s in enumerate(active)
            if ds[k,0]<=ROAD_LIMIT_M and ds[k,1]-ds[k,0]>=SEPARATION_M}
    collisions=Counter(prelim.values())
    chosen={s:n for s,n in prelim.items() if collisions[n]==1}
    matches=[]; best={}; approaches=[]
    for k,s in enumerate(active):
        sid=s['signal_id']; n=eligible[ix[k,0]]
        status=('provisional_spatial_candidate' if sid in chosen else
                'over_50m' if ds[k,0]>ROAD_LIMIT_M else
                'ambiguous_nearby_nodes' if ds[k,1]-ds[k,0]<SEPARATION_M else 'multiple_signals_same_node')
        best[sid]=(chosen.get(sid,''),status)
        for rank in range(3):
            node=eligible[ix[k,rank]]
            matches.append(dict(signal_id=sid,candidate_rank=rank+1,road_node_id=node,
                road_longitude=coords[node][0],road_latitude=coords[node][1],distance_m=round(float(ds[k,rank]),3),
                second_minus_first_distance_m=round(float(ds[k,1]-ds[k,0]),3),physical_neighbour_count=len(neighbours[node]),
                incoming_physical_links=len(incoming[node]),selected_provisional=int(rank==0 and sid in chosen),
                road_match_status=status,manual_verified=0,mapping_kind='nearest_intersection_candidate'))
        if sid in chosen:
            for link in incoming[n]:
                approaches.append(dict(signal_id=sid,road_node_id=n,link_id=link['link_id'],from_node=link['from_node'],to_node=n,
                    original_capacity_source=link['capacity_source'],outage_capacity_factor=0.5,
                    mapping_kind='all_incoming_physical_links_scenario_assumption',manual_verified=0))
    write(out/'coupling/signal_road_candidates.csv',matches)
    write(out/'coupling/signal_approach_links.csv',approaches)
    for s in inventory:
        s['road_node_id'],s['road_match_status']=best.get(s['signal_id'],('','not_in_active_primary_subset'))
    write(out/'signals/inventory.csv',inventory)
    power=read(out/'coupling/signal_power_candidates.csv')
    for p in power:p['road_node_id'],p['road_match_status']=best[p['signal_id']]
    write(out/'coupling/signal_power_candidates.csv',power)
    combined=[]
    for p in power:
        if p['selected_provisional']=='1' and p['signal_id'] in chosen:
            combined.append(dict(signal_id=p['signal_id'],road_node_id=chosen[p['signal_id']],
                feeder_id=p['feeder_id'],substation_id=p['substation_id'],load_id=p['load_id'],distribution_bus_id=p['distribution_bus_id'],
                power_distance_m=p['distance_m'],road_distance_m=next(r['distance_m'] for r in matches if r['signal_id']==p['signal_id'] and r['candidate_rank']==1),
                relation_status='provisional_for_scenario_only',real_utility_connection_verified=0,road_manual_verified=0))
    write(out/'coupling/signal_road_power_provisional.csv',combined)
    physical_ids=sorted(i for i in physical if incoming[i] and outgoing[i])
    road_tree=cKDTree(xyz([coords[i] for i in physical_ids]))
    access=[]
    for s in read(out/'power/distribution_v03/substations.csv'):
        d,j=road_tree.query(xyz([(s['Longitude'],s['Latitude'])])[0],k=3)
        for rank,(distance,index) in enumerate(zip(metres(d),j),1):
            node=physical_ids[index]
            access.append(dict(substation_id=s['Name'],candidate_rank=rank,road_node_id=node,
                longitude=coords[node][0],latitude=coords[node][1],distance_m=round(float(distance),3),
                candidate_within_500m=int(distance<=500),actual_access_verified=0,
                mapping_kind='nearest_bidirectional_physical_node_not_verified_driveway'))
    write(out/'coupling/substation_road_access_candidates.csv',access)
    # Straight endpoint geometry is for QA, never substituted for model link lengths.
    features=[]; geometry_audit=[]
    for link in links:
        u,v=int(link['from_node']),int(link['to_node'])
        a,b=coords[u],coords[v]
        chord=float(metres(np.linalg.norm(xyz([a])[0]-xyz([b])[0])))
        inferred=float(link['length_source'])*.3048
        ratio=chord/inferred if inferred else None
        features.append(dict(type='Feature',id=int(link['link_id']),properties={
            'link_id':int(link['link_id']),'from_node':u,'to_node':v,'centroid_connector':int(link['centroid_connector']),
            'geometry_kind':'straight_between_source_nodes'},geometry=dict(type='LineString',coordinates=[a,b])))
        geometry_audit.append(dict(link_id=link['link_id'],centroid_connector=link['centroid_connector'],
            endpoint_distance_m=round(chord,3),source_length_if_feet_m=round(inferred,3),
            chord_to_length_if_feet_ratio=round(ratio,6) if ratio is not None else '',
            flag_ratio_over_1_2=int(ratio is not None and ratio>1.2)))
    dump(out/'road/links.geojson',dict(type='FeatureCollection',features=features))
    dump(out/'road/nodes.geojson',dict(type='FeatureCollection',features=[dict(type='Feature',id=i,
        properties=dict(node_id=i,node_type=byid[i]['node_type']),geometry=dict(type='Point',coordinates=coords[i])) for i in sorted(coords)]))
    write(out/'road/geometry_audit.csv',geometry_audit)
    with (out/'road/Austin_node.tntp').open('w') as f:
        f.write('Node\tX\tY\n')
        for i in sorted(coords):f.write(f'{i}\t{coords[i][0]:.6f}\t{coords[i][1]:.6f}\n')
    stats={'coordinate_nodes':len(nodes),'eligible_intersection_nodes':len(eligible),'active_signals':len(active),
        'road_provisional_matches':len(chosen),'joint_road_power_provisional_matches':len(combined),
        'approach_link_rows':len(approaches),'road_match_status_counts':dict(Counter(v[1] for v in best.values())),
        'nearest_eligible_road_distance_m':quantiles(ds[:,0]),
        'sensitivity_before_collision_exclusion':{str(t):int(sum((ds[:,0]<=t)&((ds[:,1]-ds[:,0])>=SEPARATION_M))) for t in [25,50,75,100,150,200]},
        'rule':{'max_distance_m':ROAD_LIMIT_M,'min_second_candidate_margin_m':SEPARATION_M,'min_distinct_physical_neighbours':3,'multiple_signals_at_same_node':'retain candidates, leave unselected'},
        'substation_access_candidates':len(access),'nearest_substation_access_within_500m':sum(r['candidate_rank']==1 and r['candidate_within_500m']==1 for r in access),
        'physical_links_chord_longer_than_1_2_times_length_if_feet':sum(r['centroid_connector']=='0' and r['flag_ratio_over_1_2'] for r in geometry_audit),
        'manual_road_matches_verified':0,'source_datum_declared':False,'road_polylines_available':False,
        'scope':'Geographic candidates and an aggregate signal-outage scenario; not observed wiring, confirmed turning movements, or a coupled power-restoration model.'}
    dump(root/'reports/road_geography_validation.json',stats)
    meta=json.loads((out/'signals/metadata.json').read_text())
    meta.update(accepted_road_matches=0,provisional_road_matches=len(chosen),joint_road_power_provisional_matches=len(combined))
    meta['model_assumptions']['signal_road_mapping']='Nearest physical intersection within 50 m, >=20 m separation from second candidate, no multi-signal collisions; all are unreviewed candidates.'
    dump(out/'signals/metadata.json',meta)
    return stats
