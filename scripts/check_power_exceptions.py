#!/usr/bin/env python3
"""Audit original source exceptions; derived voltage-base fixes are recorded separately."""
import json,re,time,csv
from zipfile import ZipFile
from pathlib import Path
from opendssdirect import dss
from build_dataset import ROOT,RAW,OUT,dump,properties,commands,write_csv

baseline=json.loads((ROOT/'reports/opendss_validation.json').read_text())['results']
rows=list(csv.DictReader((OUT/'power/distribution_v03/feeders.csv').open()))
byid={r['feeder_id']:r for r in rows}
fixes=[];regional=[]
for r in baseline:
    if not r.get('passed') or r.get('max_voltage_pu',0)<=1.1:continue
    folder=ROOT/'runs/opendss'/r['feeder_id'];master=folder/'Master.dss';text=master.read_text()
    circuit=next(c for c in commands(text) if c.lower().startswith('new circuit.'))
    base=float(properties(circuit)['basekv'])
    m=re.search(r'(?i)(Set\s+Voltagebases\s*=\s*\[)([^]]*)(\])',text)
    bases=[float(x) for x in m.group(2).split(',')]
    if any(abs(x-base)<1e-8 for x in bases):continue
    dss(f'Compile "{master}"');before=dict(zip(dss.Circuit.AllNodeNames(),dss.Circuit.AllBusVMag()))
    changed=text[:m.start(2)]+m.group(2)+', '+str(base)+text[m.end(2):]
    derived=folder/'Master_voltagebases_fixed.dss';derived.write_text(changed)
    dss(f'Compile "{derived}"');after=dict(zip(dss.Circuit.AllNodeNames(),dss.Circuit.AllBusVMag()))
    delta=max(abs(after[k]-v) for k,v in before.items());volts=[v for v in dss.Circuit.AllBusMagPu() if v>0.1]
    result={'feeder_id':r['feeder_id'],'change':'Append source basekV to Set Voltagebases for correct per-unit reporting','original_voltagebases':bases,'added_voltagebase_kv':base,'derived_master':str(derived.relative_to(ROOT)),'converged':bool(dss.Solution.Converged()),'max_physical_voltage_change_v':delta,'min_energized_voltage_pu':min(volts),'max_voltage_pu':max(volts),'passed':bool(dss.Solution.Converged() and delta<0.1 and max(volts)<1.1)}
    fixes.append(result);print(json.dumps(result),flush=True)

with ZipFile(RAW/'power/syn-Austin-TDgrid-v03.zip') as z:
    for region in sorted({r['region'] for r in baseline if not r['passed']}):
        prefix=f'syn-austin-D_only-v03/{region}/base/opendss/'
        folder=ROOT/'runs/opendss_regions'/region;folder.mkdir(parents=True,exist_ok=True)
        for name in z.namelist():
            if name.startswith(prefix) and name.endswith('.dss'):
                relative=Path(name[len(prefix):]);assert not relative.is_absolute() and '..' not in relative.parts
                dest=folder/relative;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(name))
        result={'region':region,'scope':'Original regional OpenDSS network including upstream 230/69 kV elements and cross-feeder controls; independent regional ideal source, not full TAMU T+D co-simulation'}
        start=time.monotonic()
        try:
            dss(f'Compile "{folder / "Master.dss"}"')
            expected=sum(float(r['load_kw']) for r in rows if r['region']==region)
            kw=sum(load.kW() for load in dss.Loads)
            dark=[]
            for load in dss.Loads:
                if any(v<1 for node,v in zip(dss.CktElement.NodeOrder(),dss.CktElement.VoltagesMagAng()[::2]) if node):dark.append(load.Name())
            volts=dss.Circuit.AllBusMagPu()
            result.update(passed=bool(dss.Solution.Converged() and not dark and abs(kw-expected)<1e-5),converged=bool(dss.Solution.Converged()),loads=dss.Loads.Count(),nominal_kw=kw,expected_kw=expected,unenergized_loads=dark,min_energized_voltage_pu=min(v for v in volts if v>0.1),max_voltage_pu=max(volts))
        except Exception as e:result.update(passed=False,error=str(e))
        result['wall_seconds']=time.monotonic()-start;regional.append(result);print(json.dumps(result),flush=True)

dump(ROOT/'reports/power_exception_audit.json',{'voltagebase_corrections':fixes,'original_regional_solves':regional,'raw_sources_modified':False})
fixed={r['feeder_id']:r for r in fixes};regions={r['region']:r for r in regional}
quality=[]
for r in baseline:
    f=fixed.get(r['feeder_id']);region=regions.get(r['region'])
    quality.append(dict(feeder_id=r['feeder_id'],raw_standalone_converged=int(r.get('passed',False)),requires_regional_model=int(not r.get('passed',False)),regional_model_passed=int(bool(region and region['passed'])),voltagebase_fix_required=int(bool(f)),voltagebase_fix_passed=int(bool(f and f['passed'])),recommended_entry=('regional:'+r['region']) if not r.get('passed',False) else ('Master_voltagebases_fixed.dss' if f else 'Master.dss'),eligibility='import_solve_ready' if ((r.get('passed') and (not f or f['passed'])) or (not r.get('passed') and region and region['passed'])) else 'needs_review'))
write_csv(OUT/'power/distribution_v03/solver_quality.csv',list(quality[0]),quality)
print('ready',sum(r['eligibility']=='import_solve_ready' for r in quality),'of',len(quality),flush=True)

if any(r['eligibility']!='import_solve_ready' for r in quality):raise SystemExit('Unresolved model exceptions remain; see reports')
