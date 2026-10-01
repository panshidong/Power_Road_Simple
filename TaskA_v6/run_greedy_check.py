"""Exploratory event-specific greedy service-gain comparator, private simulator copies.

At each dispatch, maximize immediate joint-service gain among pending jobs an
available specialized crew can perform. In-progress repairs remain damaged.
No access score, travel-time denominator, or future completed repairs enter gain.
All crew execution and loss evaluation follow the original main experiment.
"""
from pathlib import Path
import os, sys, json, argparse, shutil, math

ROOT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--worker', type=int, required=True)
ap.add_argument('--limit', type=int, default=0)
args = ap.parse_args()
worker = ROOT / f'greedy_worker_{args.worker}'
if not worker.exists():
    shutil.copytree('/home/workenv/TaskA_v4/template', worker,
        ignore=shutil.ignore_patterns('__pycache__', 'work', 'nctcog', 's.txt', 'tapb_fallback.log'))
    path = worker / 'scheduler.py'
    code = path.read_text()
    needle = '            idx = _next_dispatchable_asset_index(pending, free_crews)'
    assert code.count(needle) == 1
    code = code.replace(needle,
        '            idx = (PRIORITY_SELECTOR(pending, free_crews, broken_buses, broken_links, current_time)\n'
        '                   if PRIORITY_SELECTOR is not None\n'
        '                   else _next_dispatchable_asset_index(pending, free_crews))')
    code += '\n# Optional selector for this private exploratory experiment only.\nPRIORITY_SELECTOR = None\n'
    path.write_text(code)
os.chdir(worker)
sys.path.insert(0, str(worker))
import scheduler
import resilience_measurement as rm
from task_a_criticality import FullFunctionalityValue, TaskACriticalityConfig, _asset_sort_key

base = Path('/home/workenv/TaskA_v4')
manifest = json.loads(next((base/'results/evaluation').glob('worker_*/evaluation_disasters.json')).read_text())
out = ROOT/'greedy_check'
out.mkdir(exist_ok=True)
completed_count = 0
for i, s in enumerate(manifest['scenarios']):
    if i % 4 != args.worker:
        continue
    dest = out / f"{s['scenario_id']}.json"
    if dest.exists():
        continue
    roads = [(l['u'], l['v']) for l in s['broken_links']]
    assets = s['broken_buses'] + roads
    factors = {(l['u'], l['v']): l['remaining_capacity_factor'] for l in s['broken_links']}
    value = FullFunctionalityValue(power_assets=s['broken_buses'], road_assets=roads,
        link_factors=factors, cfg=TaskACriticalityConfig())
    decisions = []
    def select(pending, free_crews, broken_buses, broken_links, current_time):
        broken = set(broken_buses) | set(broken_links)
        repaired = set(assets) - broken
        base_value = value.integrated_value(repaired)
        eligible = [j for j, a in enumerate(pending) if any(c.can_repair(a) for c in free_crews)]
        if not eligible:
            return None
        gains = {j: value.integrated_value(repaired | {pending[j]}) - base_value for j in eligible}
        selected = min(eligible, key=lambda j: (-gains[j], _asset_sort_key(pending[j])))
        assert math.isfinite(gains[selected])
        assert gains[selected] >= max(gains.values())
        decisions.append(dict(time=current_time, selected=pending[selected],
            gain=gains[selected], completed_count=len(repaired),
            candidates=[dict(asset=pending[j], gain=gains[j]) for j in eligible]))
        return selected
    scheduler.PRIORITY_SELECTOR = select
    result = rm.run_model_multi(sorted(assets, key=_asset_sort_key), result_root=str(out),
        run_dir=str(worker/'scratch_run'), message='v6 exploratory greedy service-gain check',
        Scenario=s['scenario_id'], strict=True, save_artifacts=False,
        crew_mode='specialized', power_crews=1, road_crews=1,
        preserve_sequence_order=True, broken_link_factors=factors, objective='triangle')
    scheduler.PRIORITY_SELECTOR = None
    assert len(decisions) == len(assets)
    events = result['event_log']
    assert not events[-1]['broken_buses'] and not events[-1]['broken_links']
    pl = rl = 0
    for e, n in zip(events, events[1:]):
        dt = n['time'] - e['time']
        pl += (1-rm.eval_power_resilience(e['broken_buses'])) * dt
        rl += (1-e['state']['road_func']) * dt
    assert abs(pl + rl - result['triangle_area']) < 1e-6
    dest.write_text(json.dumps(dict(scenario_id=s['scenario_id'], strategy_id='GREEDY',
        triangle_area=result['triangle_area'], power_loss=pl, road_loss=rl,
        decisions=decisions, event_log=events), indent=2))
    completed_count += 1
    if completed_count % 10 == 0 or completed_count == 1:
        print('worker', args.worker, 'completed', completed_count, flush=True)
    if args.limit and completed_count >= args.limit:
        break
print('worker', args.worker, 'DONE', flush=True)
