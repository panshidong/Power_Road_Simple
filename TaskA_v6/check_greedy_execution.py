"""Replay sampled greedy orders with the unmodified dispatcher behavior."""
from pathlib import Path
import json, csv, os, sys
ROOT=Path(__file__).resolve().parent
worker=ROOT/'greedy_worker_0'
os.chdir(worker);sys.path.insert(0,str(worker))
import resilience_measurement as rm
import scheduler
from task_a_criticality import _sorted_by_score
assert scheduler.PRIORITY_SELECTOR is None
base=Path('/home/workenv/TaskA_v4')
manifest=json.loads(next((base/'results/evaluation').glob('worker_*/evaluation_disasters.json')).read_text())
scenarios={s['scenario_id']:s for s in manifest['scenarios']}
checks=[]
for sid in ['scenario_001','scenario_046','scenario_220']:
    s=scenarios[sid];saved=json.loads((ROOT/'greedy_check'/f'{sid}.json').read_text())
    sequence=[tuple(d['selected']) if isinstance(d['selected'],list) else d['selected'] for d in saved['decisions']]
    factors={(l['u'],l['v']):l['remaining_capacity_factor'] for l in s['broken_links']}
    result=rm.run_model_multi(sequence,result_root=str(ROOT/'analysis'),run_dir=str(worker/'replay_run'),
        message='greedy order replay',Scenario=sid,strict=True,save_artifacts=False,
        crew_mode='specialized',power_crews=1,road_crews=1,preserve_sequence_order=True,
        broken_link_factors=factors,objective='triangle')
    error=result['triangle_area']-saved['triangle_area']
    assert abs(error)<1e-6,(sid,error)
    checks.append(dict(scenario_id=sid,test='greedy dispatch order replay',loss=result['triangle_area'],difference=error))
s=scenarios['scenario_001']
tables=json.loads((base/'results/tables.json').read_text())['checkpoints']['400']['tables']
assets=s['broken_buses']+[(l['u'],l['v']) for l in s['broken_links']]
def key(a):return f'road:{min(a)}-{max(a)}' if isinstance(a,tuple) else f'power:{a}'
sequence=_sorted_by_score(assets,{a:tables['JSH'][key(a)] for a in assets})
factors={(l['u'],l['v']):l['remaining_capacity_factor'] for l in s['broken_links']}
result=rm.run_model_multi(sequence,result_root=str(ROOT/'analysis'),run_dir=str(worker/'control_run'),
    message='original service Shapley control',Scenario=s['scenario_id'],strict=True,save_artifacts=False,
    crew_mode='specialized',power_crews=1,road_crews=1,preserve_sequence_order=True,
    broken_link_factors=factors,objective='triangle')
original=next(r for r in csv.DictReader(next((base/'results/evaluation').glob('worker_*/rows_scenario_001.csv')).open()) if r['strategy_id']=='JSH')
error=result['triangle_area']-float(original['triangle_area'])
assert abs(error)<1e-6,error
checks.append(dict(scenario_id='scenario_001',test='original service Shapley control',loss=result['triangle_area'],difference=error))
(ROOT/'analysis/greedy_execution_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
