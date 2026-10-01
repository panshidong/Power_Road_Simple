#!/usr/bin/env python3
"""Run a TAP-B baseline and independent OpenDSS feeder checks in isolated folders."""
from __future__ import annotations
import argparse,csv,gzip,json,math,re,shutil,subprocess,time,hashlib
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile
from build_dataset import ROOT,RAW,OUT,dump,write_csv


def tapb(executable, network=None, run_name='tapb_baseline', report_name='tapb_validation.json'):
    path=ROOT/'runs'/run_name;path.mkdir(parents=True,exist_ok=True)
    network=Path(network or OUT/'road/Austin_net.tntp').resolve()
    trip=OUT/'road/Austin_trips.tntp'
    report_path=ROOT/'reports'/report_name
    executable=Path(executable).resolve()
    args=[str(executable),'0.0001','1',str(network),str(trip),'4']
    start=time.monotonic()
    try:
        with (path/'tapb.log').open('w') as out:
            cp=subprocess.run(args,cwd=path,stdout=out,stderr=subprocess.STDOUT,timeout=180,check=False)
    except subprocess.TimeoutExpired:
        dump(report_path,{'passed':False,'reason':'timeout','command':args});raise
    log=(path/'tapb.log').read_text()
    gaps=re.findall(r'Iteration\s+(\d+):\s+gap\s+([\d.eE+\-]+)',log)
    assert cp.returncode==0 and gaps,log[-2000:]
    iteration,gap=gaps[-1];assert 0<=float(gap)<=0.0001
    flows=[]
    for line in (path/'s.txt').read_text().splitlines():
        m=re.match(r'\((\d+),(\d+)\)\s+(\S+)\s+(\S+)',line.strip())
        if m:
            u,v,f,c=m.groups();flows.append((int(u),int(v),float(f),float(c)))
    assert all(math.isfinite(f) and math.isfinite(c) and f>=-1e-8 and c>=0 for _,_,f,c in flows)
    links=list(csv.DictReader((OUT/'road/links.csv').open()))
    assert len(flows)==len(links)
    assert [(u,v) for u,v,_,_ in flows]==[(int(r['from_node']),int(r['to_node'])) for r in links], 'TAP-B output row order changed'
    balance=defaultdict(float); demand_balance=defaultdict(float)
    for u,v,f,c in flows:balance[u]+=f;balance[v]-=f
    with gzip.open(OUT/'road/od.csv.gz','rt') as inp:
        for row in csv.DictReader(inp):
            o,d,q=int(row['origin']),int(row['destination']),float(row['demand_source'])
            demand_balance[o]+=q;demand_balance[d]-=q
    max_balance_error=max(abs(balance[n]-demand_balance[n]) for n in range(1,7467))
    assert max_balance_error<1e-4, max_balance_error
    write_csv(path/'link_flows.csv.gz',['link_id','from_node','to_node','flow_source','cost_source'],
        [dict(link_id=row['link_id'],from_node=u,to_node=v,flow_source=f,cost_source=c) for row,(u,v,f,c) in zip(links,flows)])
    observed=sum(f*c for _,_,f,c in flows)
    tstt=float(re.findall(r'Aggregate TSTT:\s+([\d.eE+\-]+)',log)[-1])
    # s.txt is rounded to six decimals; compare aggregate with a rounding tolerance.
    assert math.isclose(observed,tstt,rel_tol=1e-6)
    stats={'passed':True,'network_sha256':hashlib.sha256(network.read_bytes()).hexdigest(),'trips_sha256':hashlib.sha256(trip.read_bytes()).hexdigest(),'max_node_flow_balance_error':max_balance_error,'parallel_links_preserved_in_row_order':True,'binary_sha256':hashlib.sha256(executable.read_bytes()).hexdigest(),'command':args,'iterations':int(iteration),'relative_gap':float(gap),'tstt_source_units':tstt,'tstt_from_rounded_link_file':observed,'link_flow_rows':len(flows),'wall_seconds':time.monotonic()-start,'scope':'Intact road assignment only; no coupled restoration claim.'}
    dump(report_path,stats);print(json.dumps(stats),flush=True)
    return stats


def opendss(all_feeders=False,allow_known=False):
    from opendssdirect import dss
    rows=list(csv.DictReader((OUT/'power/distribution_v03/feeders.csv').open()))
    # Median-size feeder from each of six region partitions; reproducible selection.
    selected=rows if all_feeders else [sorted([r for r in rows if r['region']==region],key=lambda r:(int(r['load_count']),r['feeder_id']))[sum(r['region']==region for r in rows)//2] for region in sorted({r['region'] for r in rows})]
    results=[]
    with ZipFile(RAW/'power/syn-Austin-TDgrid-v03.zip') as z:
        for row in selected:
            start=time.monotonic();folder=ROOT/'runs/opendss'/row['feeder_id'];folder.mkdir(parents=True,exist_ok=True)
            prefix=row['master_in_zip'].rsplit('/',1)[0]+'/'
            for name in z.namelist():
                if name.startswith(prefix) and '/' not in name[len(prefix):] and name.endswith('.dss'):
                    (folder/Path(name).name).write_bytes(z.read(name))
            result={'feeder_id':row['feeder_id'],'region':row['region']}
            try:
                dss('Clear')
                dss(f'Compile "{folder / "Master.dss"}"')
                # These are independent feeders with the original ideal source at 1.03 pu.
                # Keep original control settings; do not silently modify failed cases.
                voltages=dss.Circuit.AllBusMagPu()
                kw=kvar=actual_kw=0.0
                dark_loaded_terminals=[]
                for load in dss.Loads:
                    kw+=load.kW();kvar+=load.kvar()
                    actual_kw+=sum(dss.CktElement.Powers()[::2])
                    vmags=dss.CktElement.VoltagesMagAng()[::2]
                    if any(v<1.0 for node,v in zip(dss.CktElement.NodeOrder(),vmags) if node!=0):
                        dark_loaded_terminals.append(load.Name())
                counts_match=dss.Loads.Count()==int(row['load_count'])
                nominal_match=math.isclose(kw,float(row['load_kw']),rel_tol=1e-9) and math.isclose(kvar,float(row['load_kvar']),rel_tol=1e-9)
                converged=dss.Solution.Converged()
                finite=bool(voltages) and all(math.isfinite(v) and v>=0 for v in voltages)
                energized=[v for v in voltages if v>0.1]
                dark_nodes=[n for n,v in zip(dss.Circuit.AllNodeNames(),voltages) if v<=0.1]
                result.update(passed=bool(converged and finite and counts_match and nominal_match and not dark_loaded_terminals),converged=bool(converged),loads=dss.Loads.Count(),buses=dss.Circuit.NumBuses(),phase_nodes=dss.Circuit.NumNodes(),nominal_load_kw=kw,nominal_load_kvar=kvar,table_nominal_loads_match=nominal_match,table_load_counts_match=counts_match,min_voltage_pu_including_unloaded_stubs=min(voltages),min_energized_node_voltage_pu=min(energized) if energized else None,unenergized_phase_nodes=dark_nodes,unenergized_loads=dark_loaded_terminals,measured_load_kw=actual_kw,measured_to_nominal_load_ratio=actual_kw/kw if kw else None,max_voltage_pu=max(voltages),losses_kw=dss.Circuit.Losses()[0]/1000,source_power_kw=-dss.Circuit.TotalPower()[0])
            except Exception as e:
                result.update(passed=False,error=str(e))
            result['wall_seconds']=time.monotonic()-start;results.append(result);print(json.dumps(result),flush=True)
    report={'engine':dss.Basic.Version(),'tested_feeders':len(results),'total_available_feeders':len(rows),'passed':all(r['passed'] for r in results),'results':results,'scope':'Original independent feeder AC solves; convergence is not proof of acceptable operating limits or full T+D co-simulation.'}
    dump(ROOT/'reports/opendss_validation.json',report)
    if not report['passed']:
        failures=[r for r in results if not r['passed']]
        known={'p3uhs2_1247--p3udt31865','p5uhs1_1247--p5udt4629'}
        if not (allow_known and {r['feeder_id'] for r in failures}==known and all('Element is not set' in r.get('error','') for r in failures)):
            raise SystemExit('Some OpenDSS checks failed; see report')
        print('Two known cross-feeder references require check_power_exceptions.py; original report retains passed=false.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tapb');p.add_argument('--opendss',action='store_true');p.add_argument('--all-feeders',action='store_true');p.add_argument('--allow-known-source-exceptions',action='store_true');a=p.parse_args()
    if a.tapb:tapb(a.tapb)
    if a.opendss:opendss(a.all_feeders,a.allow_known_source_exceptions)
    if not a.tapb and not a.opendss:p.error('Choose --tapb executable and/or --opendss')
