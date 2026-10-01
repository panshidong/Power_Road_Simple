from pathlib import Path
import os, sys, json, csv, argparse, shutil
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--worker',type=int,required=True);args=ap.parse_args()
worker=ROOT/f'closure_worker_{args.worker}'
if not worker.exists():
    shutil.copytree('/home/workenv/TaskA_v4/template',worker,ignore=shutil.ignore_patterns('__pycache__','work','nctcog','s.txt','tapb_fallback.log'))
os.chdir(worker);sys.path.insert(0,str(worker))
import road_util
import resilience_measurement as rm
from task_a_criticality import _sorted_by_score
original=road_util.capacity_adjustment
def adjust(*a,**kw):
    original(*a,**kw)
    path=Path(kw.get('output_file',a[1] if len(a)>1 else ''))
    lines=path.read_text().splitlines(keepends=True)
    for i,line in enumerate(lines):
        parts=line.split()
        if len(parts)>=10 and parts[0].isdigit() and parts[4]=='9999':
            parts[4]='1000';lines[i]='\t'.join(parts)+'\n'
    path.write_text(''.join(lines))
road_util.capacity_adjustment=adjust;rm.capacity_adjustment=adjust
base=Path('/home/workenv/TaskA_v4')
tables=json.loads((base/'results/tables.json').read_text())['checkpoints']['400']['tables']
manifest=json.loads(next((base/'results/evaluation').glob('worker_*/evaluation_disasters.json')).read_text())
out=ROOT/'closure_check';out.mkdir(exist_ok=True)
for i,s in enumerate(manifest['scenarios']):
    if i%4!=args.worker:continue
    dest=out/f"{s['scenario_id']}.json"
    if dest.exists():continue
    assets=s['broken_buses']+[(l['u'],l['v']) for l in s['broken_links']]
    def key(a):return f'road:{min(a)}-{max(a)}' if isinstance(a,tuple) else f'power:{a}'
    factors={(l['u'],l['v']):l['remaining_capacity_factor'] for l in s['broken_links']}
    rows=[]
    for strategy in ['CEN','JSH','IJSH']:
        sequence=_sorted_by_score(assets,{a:tables[strategy][key(a)] for a in assets})
        result=rm.run_model_multi(sequence,result_root=str(out),run_dir=str(worker/'scratch_run'),message='v5 closure diagnostic',Scenario=s['scenario_id'],strict=True,save_artifacts=False,crew_mode='specialized',power_crews=1,road_crews=1,preserve_sequence_order=True,broken_link_factors=factors,objective='triangle')
        events=result['event_log'];pl=rl=0
        for e,n in zip(events,events[1:]):
            dt=n['time']-e['time'];pl+=(1-rm.eval_power_resilience(e['broken_buses']))*dt;rl+=(1-e['state']['road_func'])*dt
        rows.append(dict(scenario_id=s['scenario_id'],strategy_id=strategy,triangle_area=result['triangle_area'],power_loss=pl,road_loss=rl,penalty=1000))
    dest.write_text(json.dumps(rows,indent=2))
    if (i//4)%15==0:print('worker',args.worker,'completed',i//4+1,'/75',flush=True)
print('worker',args.worker,'DONE',flush=True)
