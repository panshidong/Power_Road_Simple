from pathlib import Path
import os, sys, json, time, argparse

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'case_worker'))
os.chdir(ROOT/'case_worker')
import road_util
import resilience_measurement as rm
from task_a_criticality import _sorted_by_score

ap=argparse.ArgumentParser()
ap.add_argument('--penalties',nargs='+',type=float,default=[9999,1000,100000])
ap.add_argument('--first-only',action='store_true')
args=ap.parse_args()
selection=json.loads((ROOT/'case_selection.json').read_text())
tables=json.loads(Path('/home/workenv/TaskA_v4/results/tables.json').read_text())['checkpoints']['400']['tables']
original_adjustment=road_util.capacity_adjustment
penalty=9999

def adjust(*a,**kw):
    original_adjustment(*a,**kw)
    if penalty != 9999:
        path=Path(kw.get('output_file',a[1] if len(a)>1 else ''))
        lines=path.read_text().splitlines(keepends=True)
        for i,line in enumerate(lines):
            parts=line.split()
            if len(parts)>=10 and parts[0].isdigit() and parts[4]=='9999':
                parts[4]=str(penalty)
                lines[i]='\t'.join(parts)+'\n'
        path.write_text(''.join(lines))

road_util.capacity_adjustment=adjust
rm.capacity_adjustment=adjust

def key(a):
    return f'road:{min(a)}-{max(a)}' if isinstance(a,tuple) else f'power:{a}'

out=ROOT/'case_results';out.mkdir(exist_ok=True)
for penalty in args.penalties:
    for label,case in selection['cases'].items():
        scenario=selection['scenarios'][case]
        assets=scenario['broken_buses']+[(l['u'],l['v']) for l in scenario['broken_links']]
        factors={(l['u'],l['v']):l['remaining_capacity_factor'] for l in scenario['broken_links']}
        for strategy in ['CEN','JSH','IJSH']:
            dest=out/f'{case}_{strategy}_penalty_{penalty:g}.json'
            if dest.exists():continue
            sequence=_sorted_by_score(assets,{a:tables[strategy][key(a)] for a in assets})
            start=time.monotonic()
            result=rm.run_model_multi(sequence,result_root=str(out),run_dir=str(out/f'{case}_{strategy}_{penalty:g}'),message='v5 representative recovery analysis',Scenario=f'{case}_{strategy}',strict=True,save_artifacts=False,crew_mode='specialized',power_crews=1,road_crews=1,preserve_sequence_order=True,broken_link_factors=factors,objective='triangle')
            record={k:result[k] for k in ['triangle_area','sequence','power_sequence','road_sequence','timeline','event_log']}
            record.update(case=case,label=label,strategy=strategy,penalty=penalty,elapsed_s=time.monotonic()-start)
            events=record['event_log']
            p_loss=r_loss=0.0
            for current,nxt in zip(events,events[1:]):
                dt=nxt['time']-current['time']
                pf=rm.eval_power_resilience(current['broken_buses'])
                current['power_func']=pf
                p_loss+=(1-pf)*dt;r_loss+=(1-current['state']['road_func'])*dt
            events[-1]['power_func']=rm.eval_power_resilience(events[-1]['broken_buses'])
            record.update(power_loss=p_loss,road_loss=r_loss)
            if penalty==9999:
                original=float(selection['rows'][case][strategy]['triangle_area'])
                record['difference_from_v4']=result['triangle_area']-original
                assert abs(record['difference_from_v4'])<1e-4,(case,strategy,record['difference_from_v4'])
            dest.write_text(json.dumps(record,indent=2))
            print(case,strategy,penalty,round(result['triangle_area'],3),round(record['elapsed_s'],2),flush=True)
            if args.first_only:sys.exit(0)
