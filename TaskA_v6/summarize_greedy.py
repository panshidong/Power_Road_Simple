from pathlib import Path
import json, csv
import numpy as np

ROOT=Path(__file__).resolve().parent
records=[json.loads(p.read_text()) for p in sorted((ROOT/'greedy_check').glob('scenario_*.json'))]
assert len(records)==300 and len({r['scenario_id'] for r in records})==300
original={}
for p in Path('/home/workenv/TaskA_v4/results/evaluation').glob('worker_*/rows_scenario_*.csv'):
    for r in csv.DictReader(p.open()):
        key=(r['scenario_id'],r['strategy_id'])
        assert key not in original
        original[key]=float(r['triangle_area'])
assert len(original)==900
greedy=np.array([r['triangle_area'] for r in records])
rng=np.random.default_rng(2026090701)
indices=rng.integers(0,300,size=(20000,300))
summary={'n':300,'greedy_mean':float(greedy.mean()),'greedy_median':float(np.median(greedy)),
    'comparisons':[], 'definition':'Positive differences favor the named Shapley or structural rule over greedy.'}
for strategy in ['CEN','JSH','IJSH']:
    candidate=np.array([original[(r['scenario_id'],strategy)] for r in records])
    delta=greedy-candidate
    sampled_greedy=greedy[indices].mean(axis=1)
    boot_delta=sampled_greedy-candidate[indices].mean(axis=1)
    summary['comparisons'].append(dict(strategy=strategy,mean=float(candidate.mean()),
        mean_reduction=float(delta.mean()),percent_reduction=float(100*delta.mean()/greedy.mean()),
        ci95=np.quantile(boot_delta,[.025,.975]).tolist(),
        percent_ci95=np.quantile(100*boot_delta/sampled_greedy,[.025,.975]).tolist(),
        wins=int((delta>1e-9).sum()),losses=int((delta < -1e-9).sum()),ties=int((abs(delta)<=1e-9).sum())))
for r in records:
    for d in r['decisions']:
        assert abs(d['gain']-max(c['gain'] for c in d['candidates']))<1e-12
out=ROOT/'analysis';out.mkdir(exist_ok=True)
(out/'greedy_summary.json').write_text(json.dumps(summary,indent=2))
with (out/'greedy_comparison.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['scenario_id','greedy_loss','structural_loss','service_shapley_loss','access_shapley_loss'])
    for r in records:w.writerow([r['scenario_id'],r['triangle_area']]+[original[(r['scenario_id'],s)] for s in ['CEN','JSH','IJSH']])
print(json.dumps(summary,indent=2))
