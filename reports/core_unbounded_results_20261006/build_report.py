"""Offline result export. Reads the frozen run; never imports or starts a solver."""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime
from zoneinfo import ZoneInfo
import base64, csv, hashlib, html, io, json, math, shutil, statistics, textwrap, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'runtime/output/core-unbounded-12-20261006'
OUT = Path(__file__).resolve().parent
DATA = OUT / 'data'
FIG = OUT / 'figures'
for p in (DATA, FIG): p.mkdir(exist_ok=True)
read = lambda p: json.loads(p.read_text())
dump = lambda p, x: p.write_text(json.dumps(x, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
cfg = read(RUN/'control/launch.json')['config']
validation = read(RUN/'validation.json')
cat = read(RUN/'snapshot/runtime/prepared/catalog.json')
tables = read(RUN/'tables.json')['checkpoints']
paired = read(RUN/'analysis/paired_statistics.json')
fingerprint = read(RUN/'run_manifest.json')['fingerprint']
H = 1440.0  # Former horizon, used only for offline diagnostic comparisons
assert cfg['recovery']['stop_at_horizon'] is False
font_path = '/mnt/c/Windows/Fonts/msyh.ttc'
fontManager.addfont(font_path)
fontManager.addfont('/mnt/c/Windows/Fonts/msyhbd.ttc')
plt.rcParams.update({'font.family': FontProperties(fname=font_path).get_name(), 'axes.unicode_minus':False,
                     'font.size':10, 'figure.dpi':150, 'savefig.dpi':160, 'axes.spines.top':False, 'axes.spines.right':False})

def value(v):
    return json.dumps(v, ensure_ascii=False, separators=(',',':'), allow_nan=False) if isinstance(v,(dict,list,tuple)) else v

def csvout(name, rows, fields=None):
    if fields is None: fields = list(dict.fromkeys(k for r in rows for k in r))
    with (DATA/(name+'.csv')).open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in rows: w.writerow({k:value(r.get(k,'')) for k in fields})
    return fields

def csvread(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

docs=[]; original_hashes={}
for stage in ['construct','a-main','b-main','c-main']:
    for p in sorted((RUN/'results'/stage).glob('*.json')):
        d=read(p); assert d['fingerprint']==fingerprint
        docs.append((p,d)); original_hashes[str(p.relative_to(RUN))]=sha(p)
assert len(docs)==288
evaluations=[(p,d) for p,d in docs if d['job']['stage']!='construct']
assert len(evaluations)==264
rows=[]; unfinished=[]; remaining_assets=[]; dispatch_rows=[]; event_rows=[]; zone_rows=[]
sa=[]; optimization=[]; construction=[]; construction_scores=[]; damage=[]; checks=[]
metric_names=sorted(evaluations[0][1]['payload']['metrics'])
unique_scenarios={}; by_stage={}

def prefix(d):
    j=d['job']; v=d['payload']
    return dict(stage=j['stage'],job_id=j['id'],strategy=j['strategy'],scenario_id=v['scenario_id'],scenario_seed=v['scenario_seed'])

for p,d in docs:
    j=d['job'];v=d['payload'];stage=j['stage']
    if stage=='construct':
        s=v['scenario']; initial=s['damage']; group='construction'
        construction.append(dict(job_id=j['id'],scenario_id=s['id'],scenario_seed=s['seed'],
            damaged_power=sum(a.startswith('power:') for a in initial),damaged_road=sum(a.startswith('road:') for a in initial),
            permutations=v['permutations'],wall_seconds=d['resources'].get('wall_seconds'),method=v['method']))
        for strategy, scores in v['scores'].items():
            for a,score in scores.items():
                se=v['mc_standard_errors_joint_and_access'].get(a,[None,None])
                construction_scores.append(dict(scenario_id=s['id'],scenario_seed=s['seed'],strategy=strategy,asset=a,
                    score=score,mc_se_joint=se[0],mc_se_access=se[1],permutations=v['permutations']))
        unique_scenarios[(group,s['seed'])]=(s['id'],initial)
        continue
    key=prefix(d);initial=v['events'][0]['remaining'];last=v['events'][-1];end=last['time']
    group='task_b' if stage=='b-main' else 'evaluation_ac'
    unique_scenarios[(group,v['scenario_seed'])]=(v['scenario_id'],initial)
    by_stage[(stage,v['scenario_seed'],j['strategy'])]=v
    completed={l['asset'] for l in v['dispatch'] if l['finish']<=end+1e-8}
    assert set(v['remaining'])==set(initial)-completed
    assert set(v['remaining'])==set(last['remaining'])
    assert v['complete']==(not v['remaining'])
    assert v['stop_reason']==('completed' if v['complete'] else 'horizon')
    assert v['complete'] and v['recovery_mode']=='until_complete' and v['observation_horizon_minutes'] is None
    times=np.array([e['time'] for e in v['events']]);dt=np.diff(times)
    assert np.all(dt>=0)
    e=np.array([z['electric'] for z in v['events']]);a=np.array([z['access'] for z in v['events']])
    cri=cfg['equity']['electric_weight']*e+cfg['equity']['access_weight']*a
    avg=(cri[:-1]*dt[:,None]).sum(axis=0)/end
    integrations={
        'power_loss':sum((1-z['power_func'])*t for z,t in zip(v['events'],dt)),
        'road_loss':sum((1-z['road_func'])*t for z,t in zip(v['events'],dt)),
        'weighted_triangle_area':sum((3-z['road_func']-z['weighted_power_func']-z['shelter_access'])*t for z,t in zip(v['events'],dt)),
        'min_time_avg_cri':float(min(avg)), 'mean_time_avg_cri':float(np.mean(avg)),
    }
    integrations['triangle_area']=integrations['power_loss']+integrations['road_loss']
    for k,n in integrations.items():assert math.isclose(n,v['metrics'][k],rel_tol=1e-10,abs_tol=1e-8),(j['id'],k,n,v['metrics'][k])
    base=dict(**key,complete=v['complete'],stop_reason=v['stop_reason'],
        initial_power=sum(x.startswith('power:') for x in initial),initial_road=sum(x.startswith('road:') for x in initial),
        remaining_power=sum(x.startswith('power:') for x in v['remaining']),remaining_road=sum(x.startswith('road:') for x in v['remaining']),
        remaining_assets=v['remaining'],recovery_mode=v['recovery_mode'],observation_horizon_minutes=v['observation_horizon_minutes'],sequence=v['sequence'],fallback_power=sum(x.startswith('power:') for x in v['fallback_assets']),
        fallback_road=sum(x.startswith('road:') for x in v['fallback_assets']),fallback_assets=v['fallback_assets'],
        final_power_relative=last['power_func'],final_power_absolute=last['power_absolute'],final_road_relative=last['road_func'],
        healthy_served_fraction=v['healthy_served_fraction'],**v['metrics'],**d.get('resources',{}))
    rows.append(base)
    for i,z in enumerate(v['events']):
        event_rows.append(dict(**key,event_index=i,**{k:val for k,val in z.items() if k not in ['electric','access','power_operation']},
            mean_cri=float(cri[i].mean()),min_cri=float(cri[i].min())))
    for i,l in enumerate(v['dispatch']):
        elapsed=max(0.,min(end,l['finish'])-l['time'])
        travel_elapsed=min(elapsed,l['travel_minutes']);repair_elapsed=max(0.,elapsed-l['travel_minutes'])
        dispatch_rows.append(dict(**key,dispatch_index=i,**l,repair_minutes=l['finish']-l['time']-l['travel_minutes'],
            completed_by_end=l['asset'] in completed,travel_elapsed_by_end=travel_elapsed,repair_elapsed_by_end=repair_elapsed))
    for z,zone in enumerate(v['zones']):
        ch=np.flatnonzero(cri[:,z]>=cfg['equity']['cri_threshold']);ah=np.flatnonzero(a[:,z]>=cfg['equity']['access_threshold'])
        zone_rows.append(dict(**key,zone=zone,end_minutes=end,initial_electric=float(e[0,z]),final_electric=float(e[-1,z]),
            initial_access=float(a[0,z]),final_access=float(a[-1,z]),initial_cri=float(cri[0,z]),final_cri=float(cri[-1,z]),
            time_avg_cri=float(avg[z]),cri_first_hit_minutes=float(times[ch[0]]) if len(ch) else end,cri_censored=not len(ch),
            access_first_hit_minutes=float(times[ah[0]]) if len(ah) else end,access_censored=not len(ah)))
    o=v.get('optimization')
    if o:
        optimization.append(dict(**key,**{k:val for k,val in o.items() if k not in ['trace','reference']},
            reference=o.get('reference'),**v['metrics']))
        for t in o.get('trace',[]):sa.append(dict(**key,**t))
assert not unfinished and len(sa)==504
assert sum(x['permutations'] for x in construction)==384
duplicates=sum(by_stage[('a-main',seed,s)]==by_stage[('c-main',seed,s)] for seed in range(20264001,20264025) for s in ['CEN','JSH','IJSH'])
# Stage C attaches OD metadata to its controls. Compare the actual simulation content, excluding that annotation.
same_controls=0
for seed in range(20264001,20264025):
    for s in ['CEN','JSH','IJSH']:
        av=by_stage[('a-main',seed,s)];cv=by_stage[('c-main',seed,s)]
        assert all(av[k]==cv[k] for k in ['scenario_seed','events','dispatch','sequence','metrics','complete','remaining'])
        same_controls+=1
assert same_controls==72
for (group,seed),(sid,initial) in sorted(unique_scenarios.items()):
    for asset,factor in initial.items():damage.append(dict(scenario_group=group,scenario_seed=seed,scenario_id=sid,asset=asset,
        asset_type=asset.split(':')[0],residual_factor=factor,fully_failed=(factor==0)))
score_rows=[]
for n,t in sorted(tables.items(),key=lambda x:int(x[0])):
    for strategy,scores in t['tables'].items():
        for asset,score in scores.items():score_rows.append(dict(checkpoint=int(n),strategy=strategy,asset=asset,asset_type=asset.split(':')[0],
            score=score,observed_count=t['counts'][strategy].get(asset,0)))
coverage=[]
for n,t in sorted(tables.items(),key=lambda x:int(x[0])):
    for s,scores in t['tables'].items():
        for kind in ['power','road']:
            total=sum(a.startswith(kind+':') for a in cat['assets']);nscore=sum(a.startswith(kind+':') for a in scores)
            coverage.append(dict(checkpoint=int(n),strategy=s,asset_type=kind,total_assets=total,scored_assets=nscore,missing_scores=total-nscore,coverage_fraction=nscore/total))
fallback=Counter()
for row in rows:
    if row['stage']=='a-main' and row['strategy']=='JSH':
        for k in ['power','road']:fallback[k+'_occurrences']+=row['initial_'+k];fallback[k+'_unseen']+=row['fallback_'+k]
assert dict(fallback)=={'power_occurrences':280,'power_unseen':19,'road_occurrences':197,'road_unseen':191}
groups=defaultdict(list)
for r in rows:groups[(r['stage'],r['strategy'])].append(r)
summary=[]; full_aggregate=[]
strategy_order={'a-main':['CEN','JSH','IJSH'],'b-main':['baseline_reference','single_triangle','single_gini_restore','single_maximin_time_avg_cri','single_p90_access_restore','weighted_gini_restore_l050','weighted_maximin_time_avg_cri_l050','guardrail_gini_restore'],'c-main':['CEN','JSH','IJSH','OD_CEN','OD_JSH','OD_IJSH']}
for stage,strategies in strategy_order.items():
    for strategy in strategies:
        rr=groups[(stage,strategy)]
        summary.append(dict(stage=stage,strategy=strategy,n=len(rr),incomplete=sum(not r['complete'] for r in rr),recovery_median_minutes=statistics.median(r['completion_or_censor_time'] for r in rr),recovery_max_minutes=max(r['completion_or_censor_time'] for r in rr),over_24h=sum(r['completion_or_censor_time']>H for r in rr),
            **{m:statistics.mean(r[m] for r in rr) for m in metric_names}))
        for m in metric_names:
            vals=[r[m] for r in rr]
            full_aggregate.append(dict(stage=stage,strategy=strategy,metric=m,n=len(vals),mean=statistics.mean(vals),median=statistics.median(vals),
                p95=float(np.quantile(vals,.95)),minimum=min(vals),maximum=max(vals),sample_sd=statistics.stdev(vals),incomplete=sum(not r['complete'] for r in rr)))
assert len(summary)==17
for r in csvread(RUN/'analysis/aggregate.csv'):
    found=next(x for x in full_aggregate if (x['stage'],x['strategy'],x['metric'])==(r['stage'],r['strategy'],r['metric']))
    for k in ['mean','median','p95']:assert math.isclose(found[k],float(r[k]),rel_tol=1e-10,abs_tol=1e-8)

def colname(i):
    s=''
    while i:i,r=divmod(i-1,26);s=chr(65+r)+s
    return s

def xlsx(path,sheets):
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main';rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.writestr('[Content_Types].xml','<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'+''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(sheets)+1))+'</Types>')
        z.writestr('_rels/.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="'+rel+'/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>'+''.join(f'<sheet name="{html.escape(name,quote=True)}" sheetId="{i}" r:id="rId{i}"/>' for i,(name,_) in enumerate(sheets,1))+'</sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+''.join(f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheets)+1))+'</Relationships>')
        for i,(name,source) in enumerate(sheets,1):
            with source.open(encoding='utf-8-sig',newline='') as f, z.open(f'xl/worksheets/sheet{i}.xml','w') as stream:
                w=lambda s:stream.write(s.encode('utf-8'))
                w(f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><sheetData>')
                for rownum,rr in enumerate(csv.reader(f),1):
                    cells=[]
                    for cn,c in enumerate(rr,1):
                        addr=f'{colname(cn)}{rownum}'
                        try:n=float(c);numeric=rownum>1 and math.isfinite(n) and len(c)<24
                        except ValueError:numeric=False
                        if numeric:cells.append(f'<c r="{addr}"><v>{c}</v></c>')
                        else:
                            # Excel cells have a 32767-character limit; full values remain in CSV/JSON.
                            if len(c)>32767:c=c[:32700]+' [完整内容见 CSV/JSON]'
                            cells.append(f'<c r="{addr}" t="inlineStr"><is><t xml:space="preserve">{html.escape(c)}</t></is></c>')
                    w(f'<row r="{rownum}">'+''.join(cells)+'</row>')
                w(f'</sheetData><autoFilter ref="A1:{colname(len(rr))}{rownum}"/></worksheet>')
            print('Workbook sheet:',name,flush=True)


dictionary=[
('triangle_area','综合损失面积','∫[(1−P)+(1−R)]dt；模型分钟；越小越好'),
('weighted_triangle_area','重要负荷与避难所加权损失','∫[(1−R)+(1−加权P)+(1−避难所可达性)]dt；权重与普通面积不同'),
('power_loss','电力损失面积','P 相对于可行健康 AC 供电，而非名义负荷'),
('road_loss','交通损失面积','R=健康总出行成本/当前总出行成本，截断到[0,1]'),
('gini_restore','区域 CRI 首次恢复时间 Gini','包含初始已达90%阈值、恢复时间为0的区域；不是居民收入Gini'),
('var_restore','区域首次恢复时间方差','模型分钟²'),
('p90_restore / p95_restore','区域 CRI 首次恢复时间分位数','阈值0.9；未达到者填入结束时间，并另记删失数'),
('p90_access_restore / p95_access_restore','区域可达性首次恢复时间分位数','模型分钟；大量初始达标区会产生0'),
('gini_access_restore','可达性首次恢复时间 Gini','包含初始达标区及截止时间代入的删失区'),
('min_time_avg_cri','区域时间平均 CRI 的最小值','越大越好；平均区间为该条轨迹的完成/截止时间，不统一为24h'),
('mean_time_avg_cri','区域时间平均 CRI 的均值','同上；CRI=0.133×相对供电+0.867×相对可达性'),
('maximin_time_avg_cri_loss','最弱区域时间平均 CRI 损失','1−min_time_avg_cri；越小越好'),
('gini_final_cri / mean_final_cri / min_final_cri','终态区域 CRI 的 Gini/均值/最小值','本版终态为全部修复完成'),
('shelter_access_initial / shelter_access_final','避难所起始/终态可达性','从仓库出发，健康最短路成本/当前成本，截断到[0,1]'),
('completion_or_censor_time','本版全部修复完成时间','模型分钟；须与complete/stop_reason一起读取'),
('censored_cri_zones / censored_access_zones','从未达到恢复阈值的区域数','与remaining资产数量不同'),
('final_power_relative / final_power_absolute','终态相对健康/相对名义供电比例','健康基准为名义需求的26.26915%；相对健康100%不等于名义需求全供'),
('remaining_power / remaining_road','本版完成后尚未完成的修复数量，均为0','同时包含已派出未完成与尚未派出任务'),
('residual_factor / initial_residual_factor','初始剩余能力因子','0为完全故障；大于0为部分能力；不能解读为损失百分比'),
('mean_reduction_ci95','配对平均改善的95% bootstrap区间','与该指标同单位；不是百分比区间；正值代表候选更好'),
('exact_sign_p','配对符号检验 p 值','检验胜负次数；不同于均值差检验；未作多重比较校正'),
('scenario_seed','场景唯一依据','main_00001在构表、A/C和B中对应不同seed；请连同scenario_group/stage使用'),
('power in original scenario_rows','求解器资源统计对象','不是供电性能；实际供电指标见power_loss及事件power_func'),
]

def fmt(v,n=3):
    if isinstance(v,bool):return '是' if v else '否'
    if isinstance(v,(int,np.integer)):return str(v)
    if isinstance(v,(float,np.floating)):return f'{v:.{n}f}'
    return str(v)

def htable(records,columns,cls=''):
    return '<div class="scroll"><table class="'+cls+'"><thead><tr>'+''.join('<th>'+html.escape(label)+'</th>' for _,label in columns)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(r.get(k,'')))+'</td>' for k,_ in columns)+'</tr>' for r in records)+'</tbody></table></div>'

def picture(name,caption):
    b=base64.b64encode((FIG/(name+'.png')).read_bytes()).decode()
    return f'<figure><img alt="{html.escape(caption)}" src="data:image/png;base64,{b}"><figcaption>{caption}</figcaption></figure>'

# Revalidate the complete saved trajectories without calling audit(), which writes
# into the run directory. Only the pure trajectory checker is called here.
import importlib.util
spec=importlib.util.spec_from_file_location('offline_acceptance',RUN/'control/audit_completed_run.py')
auditor=importlib.util.module_from_spec(spec);spec.loader.exec_module(auditor)
acceptance=read(RUN/'control/completion_acceptance.json')
assert acceptance['passed'] and acceptance['source_sha256']==original_hashes
trajectory_checks=[dict(**prefix(d),**auditor.audit_trajectory(d['payload'],cfg)) for _,d in evaluations]
assert all(r['censored_cri_zones']==r['censored_access_zones']==0 for r in rows)
assert len(metric_names)==22
print('All 264 trajectories and 22 metrics rechecked; source hashes match final acceptance.',flush=True)

OLD=ROOT/'runtime/output/pilot-core-tapbfix-20261005'
old_docs={};old_hashes={}
for p in sorted((OLD/'results').glob('*/*.json')):
    d=read(p);old_hashes[str(p.relative_to(OLD))]=sha(p)
    if d['job']['stage']!='construct':
        k=(d['job']['stage'],d['payload']['scenario_seed'],d['job']['strategy']);old_docs[k]=d
assert len(old_hashes)==288 and len(old_docs)==264
assert old_hashes==read(ROOT/'reports/core_results_20261006/verification.json')['original_results_sha256']
comparisons=[];tails=[];over24=[];crew_rows=[]
for p,d in evaluations:
    v=d['payload'];key=prefix(d);end=v['events'][-1]['time']
    old=old_docs[(key['stage'],key['scenario_seed'],key['strategy'])]['payload']
    assert old['events'][0]['remaining']==v['events'][0]['remaining']
    before={m:0. for m in ['power_loss','road_loss','weighted_triangle_area']}
    tail={m:0. for m in before}
    for e,nxt in zip(v['events'],v['events'][1:]):
        rates={'power_loss':1-e['power_func'],'road_loss':1-e['road_func'],
               'weighted_triangle_area':3-e['road_func']-e['weighted_power_func']-e['shelter_access']}
        for m,rate in rates.items():
            before[m]+=rate*max(0.,min(H,nxt['time'])-e['time'])
            tail[m]+=rate*max(0.,nxt['time']-max(H,e['time']))
    before['triangle_area']=before['power_loss']+before['road_loss']
    tail['triangle_area']=tail['power_loss']+tail['road_loss']
    for m in before:assert math.isclose(before[m]+tail[m],v['metrics'][m],abs_tol=1e-8)
    at24=next((e for e in reversed(v['events']) if e['time']<=H),v['events'][0])
    tailrow=dict(**key,end_minutes=end,over_24h=end>H,extra_duration_minutes=max(0,end-H),
        remaining_at_former_horizon=at24['remaining'],remaining_count_at_former_horizon=len(at24['remaining']),
        **{'first_24h_'+m:x for m,x in before.items()},**{'after_24h_'+m:x for m,x in tail.items()})
    tails.append(tailrow)
    comp=dict(**key,old_complete=old['complete'],new_complete=v['complete'],old_remaining=len(old['remaining']),new_remaining=0,
        sequence_equal=old['sequence']==v['sequence'],old_end_minutes=old['events'][-1]['time'],new_end_minutes=end)
    for m in metric_names:
        comp['old_'+m]=old['metrics'][m];comp['new_'+m]=v['metrics'][m];comp['delta_'+m]=v['metrics'][m]-old['metrics'][m]
    comparisons.append(comp)
    crews={}
    for kind in ['power','road']:
        legs=[l for l in v['dispatch'] if l['crew']==kind+'-0']
        travel=sum(l['travel_minutes'] for l in legs)
        repair=sum(l['finish']-l['time']-l['travel_minutes'] for l in legs)
        finish=max(l['finish'] for l in legs)
        crews[kind]=dict(crew=kind+'-0',jobs=len(legs),travel_minutes=travel,repair_minutes=repair,
            final_finish_minutes=finish,idle_before_crew_finish_minutes=max(0,finish-travel-repair),
            idle_to_all_repairs_complete_minutes=max(0,end-travel-repair))
        crew_rows.append(dict(**key,**crews[kind]))
    if end>H:
        over24.append(dict(**key,end_minutes=end,end_hours=end/60,extra_minutes=end-H,remaining_final=0,
            remaining_at_24h=len(at24['remaining']),old_remaining=len(old['remaining']),
            power_jobs=crews['power']['jobs'],power_travel_minutes=crews['power']['travel_minutes'],
            power_repair_minutes=crews['power']['repair_minutes'],power_finish_minutes=crews['power']['final_finish_minutes'],
            road_finish_minutes=crews['road']['final_finish_minutes'],
            triangle_area=v['metrics']['triangle_area'],after_24h_triangle_area=tail['triangle_area']))
assert len(over24)==18 and {r['scenario_seed'] for r in over24}=={20264015,20264016}
assert sum(not r['old_complete'] for r in comparisons)==18

comparison_summary=[]
for s in summary:
    cc=[r for r in comparisons if (r['stage'],r['strategy'])==(s['stage'],s['strategy'])]
    tt=[r for r in tails if (r['stage'],r['strategy'])==(s['stage'],s['strategy'])]
    oldmean=statistics.mean(r['old_triangle_area'] for r in cc)
    old_incomplete=sum(not r['old_complete'] for r in cc)
    comparison_summary.append(dict(stage=s['stage'],strategy=s['strategy'],n=len(cc),old_incomplete=old_incomplete,new_incomplete=0,
        identical_sequences=sum(r['sequence_equal'] for r in cc),old_triangle_mean=oldmean,new_triangle_mean=s['triangle_area'],
        delta_triangle_mean=s['triangle_area']-oldmean,delta_percent=(s['triangle_area']/oldmean-1)*100,
        new_first_24h_triangle_mean=statistics.mean(r['first_24h_triangle_area'] for r in tt),
        new_after_24h_triangle_mean=statistics.mean(r['after_24h_triangle_area'] for r in tt),
        old_end_or_horizon_mean_minutes=statistics.mean(r['old_end_minutes'] for r in cc),
        new_recovery_mean_minutes=s['completion_or_censor_time']))

native_rows=[];native_directories=[]
for p in sorted(RUN.glob('scratch/construct/*/dss/*/recovery.json')):
    r=read(p);assert r['status']=='recovered_and_confirmed' and r['physical_differences']==[]
    native_rows.append(dict(job_id=p.parents[2].name,source=str(p.relative_to(RUN)),**r))
    native_directories.extend([p.parent,Path(r['recovery_directory']),Path(r['confirmation_directory'])])
assert len(native_rows)==2
guard=[r for r in optimization if r['strategy']=='guardrail_gini_restore']
guard_bad=sum(r['guardrail_violation']>1e-12 for r in guard)
main_pair=[r for r in paired if r['metric']=='triangle_area']
report_time=datetime.now(ZoneInfo('America/Chicago')).strftime('%Y-%m-%d %H:%M %Z')
dictionary.extend([
    ('recovery_mode / observation_horizon_minutes','运行停止条件','until_complete / null；仅全部修复后结束。本版没有观察窗。'),
    ('first_24h_* / after_24h_*','本版轨迹的面积分段','在本版已保存轨迹上以1440分钟分段积分；不重新优化或求解；不等于旧版结果。'),
    ('old_* / new_* / delta_*','新旧版本配对','old是8 worker、TAP-B4线程的限窗版；new是12 worker、TAP-B2线程的无窗版；delta=new−old。'),
    ('sequence_equal','新旧优先序是否一致','比较完整sequence字段；数值细微变化可影响分数、排序或SA路径；不意味着实际派工完全相同。'),
    ('end_minutes / end_hours','全部修复所需模型时间','机器计算用时单列为active_seconds；二者不可混淆。'),
    ('guardrail_violation','软约束违反量','0才达标；正值表示该场景最终解仍违反Gini目标，优化任务成功不等于约束达标。'),
])
exports={
    'scenario_results_all_264':rows,'strategy_summary_17':summary,'all_metric_aggregates':full_aggregate,
    'over_24h_18':over24,'dispatch_all':dispatch_rows,'events_all':event_rows,'crew_time_all':crew_rows,
    'zone_results_all':zone_rows,'construction_24':construction,'construction_scores':construction_scores,
    'damage_scenarios':damage,'score_tables_all_checkpoints':score_rows,'coverage_by_asset_type':coverage,
    'sa_iterations_504':sa,'optimization_48':optimization,'paired_statistics_90':paired,
    'old_vs_new_all_264':comparisons,'old_vs_new_strategy_17':comparison_summary,'new_trajectory_24h_decomposition':tails,
    'native_recoveries_2':native_rows,'trajectory_acceptance_264':trajectory_checks,'guardrail_6':guard,
    'metric_dictionary':[dict(field=k,label=l,definition=d) for k,l,d in dictionary],
}
for name,rr in exports.items():csvout(name,rr)
for p in sorted((RUN/'analysis').glob('*.csv')):shutil.copyfile(p,DATA/('original_'+p.name))
verification=dict(passed=True,report_time=report_time,run=str(RUN),fingerprint=fingerprint,
    source_result_sha256=original_hashes,previous_result_sha256=old_hashes,source_hashes_match_completion_acceptance=True,
    all_264_repairs_completed=True,metrics_recomputed_per_trajectory=22,
    regional_physical_checks=sum(r['regional_physical_checks'] for r in trajectory_checks),
    all_85_original_aggregate_rows_recomputed=True,identical_ac_controls=72,
    same_old_new_damage_scenarios=264,over_24h_trajectories=18,unique_over_24h_scenarios=2,
    unchanged_sequences=sum(r['sequence_equal'] for r in comparisons),
    native_recoveries_independently_confirmed=2,guardrail_violations=guard_bad,
    csv_rows={k:len(v) for k,v in exports.items()},solvers_started=0,
    scope='Offline result/program verification; not a claim of calibrated real-world recovery times.')
dump(OUT/'verification.json',verification)
print(json.dumps(dict(summary=summary,comparison_summary=comparison_summary,over24=over24,
    main_pairs=main_pair,guardrail_violations=guard_bad,export_counts=verification['csv_rows']),ensure_ascii=False),flush=True)

def short(s):
    return s.replace('baseline_reference','baseline').replace('maximin_time_avg_cri','maximin').replace('gini_restore','Gini').replace('p90_access_restore','P90')

def savefig(name,f):
    f.savefig(FIG/(name+'.png'),bbox_inches='tight');f.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(f)

fig,axes=plt.subplots(1,3,figsize=(14,5.6),gridspec_kw={'width_ratios':[.9,1.5,1.2]})
for ax,stage in zip(axes,strategy_order):
    ss=[r for r in summary if r['stage']==stage];ys=np.arange(len(ss))
    vals=[r['triangle_area'] for r in ss]
    ax.barh(ys,vals,color=['#c8752b' if r['triangle_area']==min(vals) else '#33778d' for r in ss])
    ax.set_yticks(ys,[short(r['strategy']) for r in ss]);ax.invert_yaxis();ax.set_xlim(0,max(vals)*1.25)
    for y,v in zip(ys,vals):ax.text(v+1,y,f'{v:.2f}',va='center',fontsize=9)
    ax.set_title(stage.upper()+f" · n={ss[0]['n']}");ax.set_xlabel('平均综合损失面积（模型分钟，越低越好）');ax.grid(axis='x',alpha=.2)
fig.suptitle('无观察窗：17组策略的全部修复结果',size=16);fig.tight_layout();savefig('strategy_loss',fig)

fig,ax=plt.subplots(figsize=(12,8))
for y,s in enumerate(summary):
    rr=groups[(s['stage'],s['strategy'])];times=[r['completion_or_censor_time']/60 for r in rr]
    jitter=np.linspace(-.15,.15,len(times))
    ax.scatter(times,y+jitter,s=18,alpha=.5,color='#33778d')
    ax.scatter([statistics.mean(times)],[y],marker='D',s=35,color='#c8752b',zorder=3)
ax.axvline(24,color='#a43c42',linestyle='--');ax.text(24.2,-.6,'原24小时观察窗',color='#a43c42')
ax.set_yticks(range(17),[r['stage'][0].upper()+' / '+short(r['strategy']) for r in summary]);ax.invert_yaxis()
ax.set_xlabel('完成全部修复所需模型小时；蓝点为场景，橙色菱形为均值');ax.grid(axis='x',alpha=.2)
ax.set_title('264条评估均已完成；18条超过原观察窗，来自相同的2场灾害');fig.tight_layout();savefig('recovery_duration',fig)

fig,axes=plt.subplots(2,1,figsize=(13,7.5),sharex=True)
for ax,seed in zip(axes,[20264015,20264016]):
    for y,s in enumerate(strategy_order['c-main']):
        v=by_stage[('c-main',seed,s)]
        for leg in v['dispatch']:
            if leg['crew']!='power-0':continue
            ax.broken_barh([(leg['time']/60,leg['travel_minutes']/60)],(y-.32,.64),facecolors='#67a1b7')
            repair=leg['finish']-leg['time']-leg['travel_minutes']
            ax.broken_barh([((leg['time']+leg['travel_minutes'])/60,repair/60)],(y-.32,.64),facecolors='#d37c2c')
        ax.text(v['events'][-1]['time']/60+.25,y,f"{v['events'][-1]['time']/60:.2f}h",va='center',fontsize=8)
    ax.axvline(24,color='#ad313b',linestyle='--',linewidth=1.5)
    ax.set_yticks(range(6),strategy_order['c-main']);ax.invert_yaxis();ax.set_xlim(0,43)
    ax.set_title(f'场景 {seed-20264000}：电力队出行（蓝）与修复（橙）；红线仅标原观察窗')
    ax.grid(axis='x',alpha=.2)
axes[-1].set_xlabel('灾后模型小时；全线来自本次实际模拟派工记录，直到全部修复完成')
fig.tight_layout();savefig('completed_timeline',fig)

fig,ax=plt.subplots(figsize=(11,5))
cc=[r for r in comparison_summary if r['stage']=='c-main'];x=np.arange(len(cc));w=.36
ax.bar(x-w/2,[r['old_triangle_mean'] for r in cc],w,color='#9aa9b0',label='旧版限窗结果')
ax.bar(x+w/2,[r['new_first_24h_triangle_mean'] for r in cc],w,color='#33778d',label='新版轨迹前24h')
ax.bar(x+w/2,[r['new_after_24h_triangle_mean'] for r in cc],w,bottom=[r['new_first_24h_triangle_mean'] for r in cc],color='#c8752b',label='新版轨迹24h后')
ax.set_xticks(x,[r['strategy'] for r in cc]);ax.set_ylabel('平均综合损失面积（模型分钟）');ax.legend()
ax.set_title('损失变化分解：新版完整积分 = 新版前24h + 24h后的恢复损失');ax.grid(axis='y',alpha=.2)
fig.tight_layout();savefig('horizon_decomposition',fig)

sheets=[('说明与指标',DATA/'metric_dictionary.csv'),('策略汇总17',DATA/'strategy_summary_17.csv'),
    ('全部评估264',DATA/'scenario_results_all_264.csv'),('超过24h已完成18',DATA/'over_24h_18.csv'),
    ('全部指标统计374',DATA/'all_metric_aggregates.csv'),('配对统计90',DATA/'paired_statistics_90.csv'),
    ('新旧版本逐项264',DATA/'old_vs_new_all_264.csv'),('新旧策略对照17',DATA/'old_vs_new_strategy_17.csv'),
    ('新版24h损失分解',DATA/'new_trajectory_24h_decomposition.csv'),('全部派工5221',DATA/'dispatch_all.csv'),
    ('全部事件5485',DATA/'events_all.csv'),('全部队伍耗时528',DATA/'crew_time_all.csv'),
    ('区域结果156816',DATA/'zone_results_all.csv'),('构表24',DATA/'construction_24.csv'),
    ('构表边际贡献2790',DATA/'construction_scores.csv'),('全部分数表28996',DATA/'score_tables_all_checkpoints.csv'),
    ('覆盖率36',DATA/'coverage_by_asset_type.csv'),('损伤场景1058',DATA/'damage_scenarios.csv'),
    ('SA迭代504',DATA/'sa_iterations_504.csv'),('优化48',DATA/'optimization_48.csv'),('公平约束6',DATA/'guardrail_6.csv'),
    ('原生故障恢复2',DATA/'native_recoveries_2.csv'),('逐轨迹验收264',DATA/'trajectory_acceptance_264.csv')]
for name,filename in [('原聚合85','aggregate.csv'),('原评估264','scenario_rows.csv'),('排序稳定性','table_stability.csv'),
    ('序列一致性','sequence_agreement.csv'),('B帕累托','task_b_pareto.csv'),('OD覆盖','od_coverage.csv'),('OD分数','od_asset_scores.csv')]:
    sheets.append((name,DATA/('original_'+filename)))
xlsx(OUT/'Austin_unbounded_all_results.xlsx',sheets)
with zipfile.ZipFile(OUT/'Austin_unbounded_all_results.xlsx') as z:
    assert z.testzip() is None
    for name in z.namelist():
        if name.endswith(('.xml','.rels')):
            with z.open(name) as f:
                for _,el in ET.iterparse(f,events=('end',)):el.clear()
print('Workbook verified:',len(sheets),'sheets',flush=True)

overview=[
    '本次无观察窗核心版已完整结束：288/288任务成功，0任务失败或缺失。264/264评估轨迹全部完成修复，remaining均为空，stop_reason均为completed；unfinished repairs=0。',
    '计算时间：2026-10-06 09:18:33至20:59:41（America/Chicago），记录的有效运行时间42067.292秒，约11小时41分。使用12个worker，每个TAP-B调用2线程。',
    '18条轨迹超过原1440分钟观察窗，来自场景15、16两场灾害。最长为场景15/OD_JSH，2357.196分钟，即约39小时17分。这是模拟的灾后恢复时间，不能当作机器计算耗时。',
    '24个构表任务包含384次排列；A=72、B=48、C=144，共264条评估；B共504次SA迭代。所有结果、原物理验收、独立终审、两次原生故障恢复证据及完整导出一并交付。',
    '本版有2次P3U原生SIGSEGV，各经一次受限恢复及第二次独立确认后成功；没有未处理任务错误。该结果不证明底层原生崩溃根因已解决。',
    '这是用户确认的核心288任务，不是完整研究矩阵。灾害规模、健康基准、交通时间与现实设施的校准仍未完成。报告过程只读取已保存结果，没有追加求解。',
]
methods=[
    '任务A：比较CEN、JSH、IJSH三种排序，在24个共同灾害场景上评价。CEN采用中心性；JSH采用抽样排列的联合功能边际贡献；IJSH在此基础上加入可达性边际贡献，本核心策略alpha=1。',
    '贡献表由另24个构表场景生成，每场景16次排列，保留8/16/24场景检查点。JSH与IJSH使用相同抽样联盟；未见资产按既定规则回退中心性并记录。',
    '任务B：6个独立于A/C的场景，各比较8个策略。基准为按场景种子打乱的顺序；7种SA目标每次12迭代，包括损失、Gini、最弱区域CRI、P90可达性及两种lambda=0.5加权目标。',
    '任务C：使用与A相同的24个场景，比较3个基础策略与3个OD道路排序策略。OD版本电力排序沿用对应基础策略，道路按所选OD路径分数排序；默认k=3、3个OD对、9条路径，75个道路组有非零分数。',
    '构表种子20263001–20263024；A/C种子20264001–20264024；B种子20260325–20260330。main_00001等编号跨任务重复，应以任务组和种子区分。A/C的72条基础对照数值完全一致，不是独立重复样本。',
    '本版stop_at_horizon=false、recovery_mode=until_complete、observation_horizon_minutes=null；配置保留的1440只供兼容，不再用于终止本次模拟。保留原AC/交通验收、资源保护及单次求解超时。',
]
notes={
    'a-main':[
        'CEN的平均综合损失面积最低，为54.507；JSH为82.582（比CEN高51.51%），IJSH为79.876（高46.54%）。JSH虽然17/24场景优于CEN，但少数较大的损失使均值更差。',
        'CEN→JSH的平均改善95% bootstrap区间为[-79.267, 8.555]；CEN→IJSH为[-76.004, 10.719]，均跨0。这里正值表示改善，不能据此断言贡献排序有总体优势。',
        '本次结果包含24小时以后的全部恢复损失；尤其场景15的尾部损失较大。所有策略均完成修复，未完成数量已不能作为策略优劣的区分。'],
    'b-main':[
        'single_triangle平均损失42.437，基准58.598，降低27.58%；6/6场景均改善，配对平均改善95%区间[3.805, 33.956]模型分钟。',
        '不同公平目标带来不同取舍：single_gini_restore的平均Gini为0.8628；single_maximin的最弱区平均CRI为0.8637；single_P90平均P90恢复16.846分钟，但平均综合损失升至67.972。',
        'guardrail仍有5/6场景违反Gini目标上限，仅1/6达标。当前实现是损失加违反量罚项的软约束，且每次仅12迭代；任务正常结束不能被写成公平约束已满足。',
        '同场景用于优化和评价，只有6个样本。这些改进描述本次有限搜索结果，不能当作在未见灾害上的泛化保证。'],
    'c-main':[
        'OD_CEN、OD_JSH、OD_IJSH平均损失分别为62.968、91.021、88.315；相对对应基础策略分别增加15.52%、10.22%、10.57%。',
        '三项配对平均改善95%区间均在0以下，表明本次所选OD设置的样本平均损失更高；未做多重比较校正，且覆盖仅3个OD对，不能概括所有OD配置。',
        'A/C的基础对照相同。若两处bootstrap区间略有不同，是原分析为不同阶段生成不同重采样序列；不是重复模拟出了不同基础轨迹。'],
}
why_complete=[
    '上一版18条unfinished repairs都在1440分钟观察窗按规则终止。本版取消该终止条件，按照下一次实际修复完成事件继续推进，直到每个初始受损资产都完成；并未删除剩余资产或把部分完成伪装为完成。',
    '原18条对应的本版轨迹均完成：场景15约35.78–39.29小时，场景16约26.66–26.73小时。全部道路更早修完，最后完成的是电力队的任务。',
    '这两场景各有15个受损电力资产，单支电力队总纯修复时间300分钟（15×20）。场景15的出行累计约1846.9–2057.2分钟，场景16约1299.5–1303.5分钟；模型恢复时间主要消耗在出行。',
    '此前旧版路线审计已发现健康交通基准中存在很长的拥堵行程。本版时间线直接来自新结果派工日志，证明了模型内的时间消耗；不能据此证明现实中必然需要这些小时，仍须单独审查交通和设施接入假设。',
    '解除观察窗使缺失的后续修复被算完，但不会改善模型中的交通速度或缩短既定维修时长。增加worker只影响计算调度。',
]
comparison_notes=[
    '新旧版逐项配对使用同一灾害种子和初始损伤；264条全部核实一致。旧版为8 worker、每TAP-B4线程、1440分钟窗口；本版12 worker、每TAP-B2线程、无观察窗。',
    f"{verification['unchanged_sequences']}/264条优先序完全一致，其余存在变化。TAP-B并行线程变化带来容差内的数值差异，可影响构表排序、SA接受路径和最终顺序，因此总版本差异不能全归于取消窗口。",
    '为单独识别尾部积分，本报告在本版已保存轨迹上以1440分钟分段，不重跑求解或优化。新版完整损失=本版前24h损失+本版24h后损失；前24h部分不必等于旧版结果。',
    'A中CEN的24h后损失平均仅0.0423；JSH为9.4013，IJSH为8.7158。JSH新旧平均损失增加9.0928，但本版尾部面积为9.4013，两者差额来自前24h轨迹/数值变化。',
    'B的所有轨迹都在24h以内，因此B的新旧差异不来自截掉的24h后损失。全部22个指标的新旧值、差额、顺序一致性和分段积分已逐项输出。',
]
correctness=[
    '物理验收于09:26:33通过；独立完成审计于21:00:38通过。报告生成时再次核对288份结果SHA256与终审记录一致，并逐条检查264份评估的全部22项指标。',
    f"独立复核包括：全部任务与计划、受损资产无遗漏/重复、派工时间/旅行/维修时长、队伍不重叠、事件remaining因果、事件严格推进、最终健康功能，以及{verification['regional_physical_checks']}项区域物理状态检查。",
    '原AC收敛、控制器稳定、电压、热限、功率平衡及交通gap/节点守恒门槛均保留；没有为完成报告放宽任何阈值。264条终态相对健康电力和道路功能均为1，CRI/可达性阈值删失区域均为0。',
    '构表场景21和24各发生一次P3U原生子进程SIGSEGV(-11)。两次均status=recovered_and_confirmed、maximum_recovery_attempts=1、physical_differences=[]；原失败日志、恢复与确认输入输出均随ZIP保存。',
    '启动前证据包括76项常规测试、5项原生测试、3项完成审计测试，以及12路并发TAP-B预检；道路18710条原始记录的列数及长度/速度/自由流时间单位关系已检查。测试通过说明执行与已检查规则一致，不能替代现实校准。',
    '与旧版健康和故障/降额/修复状态的AC物理字段对照一致。12路交通预检全部通过原gap=1e-4，和4线程基准的最大TSTT相对差约0.0122%；这些预检值不是全部研究状态的误差上界。',
]
limits=[
    '健康电力基准已有降载：实际可行供电814153.921 kW，名义需求3099277.740 kW，仅26.26915%。本报告相对供电P=1表示恢复到该可行基准，不能写成全部名义需求都得到满足。',
    '贡献表覆盖不足：JSH/IJSH在24构表场景后有分数的电力资产119/128、道路组193/8414。A/C独立场景共有280次电力损伤与197次道路损伤，其中19次电力、191次道路未见贡献分数，按规则回退中心性；道路回退率96.954%。',
    '灾害规模仍未校准：保留交接中的main分布，8–15个变电站组、3–13个道路组，含完全故障及部分剩余能力。这是现有算法情景，不是Austin实测灾损率；没有恢复已撤回的损坏百分比建议。',
    '静态交通、OD需求/容量时间口径、站点接入、单支电力队和单支道路队、固定修复时长及健康供电基准都需要现实校准。输入单位自洽并不证明上述物理假设已充分代表现实。',
    '594个模型区域，CRI=0.133×相对供电+0.867×相对可达性，首次达到0.9计为恢复；初始达标者恢复时间为0，首次达标不保证之后始终达标。恢复时间Gini不是收入Gini。时间平均CRI以每条轨迹实际完成时刻为分母。',
    'A/C各24个场景且共用对照，B只有6个；SA每次12迭代，guardrail为软约束。配对bootstrap和符号检验未做多重比较校正；没有恢复全研究矩阵中的广泛敏感性、迁移、鲁棒性与扩展比较。',
]

shortcols=[('strategy','策略'),('n','n'),('triangle_area','损失面积↓'),('weighted_triangle_area','加权损失↓'),
    ('gini_restore','恢复Gini↓'),('min_time_avg_cri','最弱区平均CRI↑'),('p90_access_restore','P90可达/min↓'),
    ('completion_or_censor_time','平均修复/min'),('incomplete','未完成')]
overcols=[('stage','任务'),('scenario_seed','种子'),('strategy','策略'),('end_minutes','全部修复/min'),('end_hours','小时'),
    ('remaining_at_24h','24h时剩余'),('remaining_final','最终剩余'),('power_travel_minutes','电力队出行/min'),('after_24h_triangle_area','24h后损失')]
compcols=[('stage','任务'),('strategy','策略'),('old_triangle_mean','旧平均损失'),('new_triangle_mean','新平均损失'),
    ('delta_triangle_mean','差额'),('new_after_24h_triangle_mean','本版24h后损失'),('identical_sequences','顺序一致数'),('n','n')]
paircols=[('stage','任务'),('reference','参考'),('candidate','候选'),('n','n'),('mean_reduction','平均改善'),
    ('mean_reduction_ci95','95%区间'),('wins','胜'),('losses','负'),('ties','平'),('exact_sign_p','符号p')]
def paras(ps):return ''.join('<p>'+html.escape(p)+'</p>' for p in ps)
parts=['<header><div class="eyebrow">AUSTIN · 12 WORKERS · UNTIL COMPLETE</div><h1>无观察窗核心版<br>完整结果与验收报告</h1><p>288个任务 · 264条完整恢复轨迹 · 未完成修复为0</p><p class="small">生成时间：'+report_time+'；数据来自2026-10-06 20:59:41完成的新版运行</p></header>',
    '<nav><a href="#summary">完成情况</a><a href="#strategies">全部策略</a><a href="#completion">18条长恢复</a><a href="#compare">新旧对照</a><a href="#audit">验收</a><a href="#all">全部264结果</a><a href="#download">下载</a></nav>',
    '<section id="summary"><h2>1. 完成情况</h2>'+paras(overview)+'</section>',
    '<section><h2>2. 本次到底完成了什么</h2>'+paras(methods)+'</section>',
    '<section id="strategies"><h2>3. 全部17组策略结果</h2><p>表格为各场景均值；损失/时间单位为模型分钟，↑越大越好、↓越小越好。表中保留3位小数，Excel/CSV保留完整数值。</p>']
for stage in strategy_order:
    parts.append('<h3>'+stage.upper()+'</h3>'+htable([r for r in summary if r['stage']==stage],shortcols)+paras(notes[stage]))
parts.append(picture('strategy_loss','全部17组平均损失；A/C相同基础对照不能当作独立样本。')+'</section>')
parts.append('<section id="completion"><h2>4. 原先的unfinished repairs现在如何</h2>'+paras(why_complete)+htable(over24,overcols)+picture('completed_timeline','本次完整派工时间线；A基础策略与C对应对照相同，此图只画C的6种策略。')+picture('recovery_duration','全部264条轨迹的完成时间分布。')+'</section>')
parts.append('<section id="compare"><h2>5. 与上一版的全部对照</h2>'+paras(comparison_notes)+htable(comparison_summary,compcols)+picture('horizon_decomposition','仅在本版轨迹上分段，明确区分新旧版本总差异与本版尾部损失。')+'</section>')
parts.append('<section id="audit"><h2>6. 程序验收与两次原生故障</h2>'+paras(correctness)+htable(native_rows,[('job_id','构表任务'),('region','区域'),('initial_returncode','首次退出码'),('status','恢复状态'),('wall_seconds','恢复流程秒'),('physical_differences','独立确认差异')])+'<p>本次没有剩余修复，也没有未处理任务失败；这与“从未发生原生崩溃”是不同的事实。</p></section>')
parts.append('<section><h2>7. 必须保留的解释边界</h2>'+paras(limits)+'<h3>公平性guardrail逐场景</h3>'+htable(guard,[('scenario_seed','种子'),('gini_restore','最终Gini'),('guardrail_limit','Gini上限'),('guardrail_violation','违反量'),('triangle_area','综合损失')])+'</section>')
allcols=[('stage','任务'),('scenario_seed','种子'),('strategy','策略'),('triangle_area','损失面积'),('weighted_triangle_area','加权损失'),('gini_restore','Gini'),('min_time_avg_cri','最弱区CRI'),('p90_access_restore','P90/min'),('completion_or_censor_time','修复/min'),('remaining_power','剩余电力'),('remaining_road','剩余道路')]
parts.append('<section id="all"><h2>8. 全部264条评估结果</h2><p>输入任务、策略或种子筛选。全部22个原始指标及新旧差额见Excel/CSV；每区域、每事件完整向量见原始JSON。</p><input id="filter" placeholder="例如 OD_JSH、20264015、b-main" aria-label="筛选结果"><span id="count">264条</span>'+htable(rows,allcols,'allresults')+'</section>')
parts.append('<section><h2>9. 构表、配对统计与指标词典</h2><details><summary>24个构表任务 / 384次排列</summary>'+htable(construction,[('scenario_id','编号'),('scenario_seed','种子'),('damaged_power','损伤电力'),('damaged_road','损伤道路'),('permutations','排列'),('wall_seconds','计算秒')])+'</details>')
parts.append('<details><summary>全部90条配对统计</summary><p>平均改善为正表示候选更好；区间单位与指标相同，非百分比。未作多重比较校正；符号p检验胜负，不直接检验均值差。</p>'+htable(paired,[('metric','指标')]+paircols)+'</details>')
parts.append('<details><summary>指标与导出字段词典</summary>'+htable(exports['metric_dictionary'],[('field','字段'),('label','中文含义'),('definition','计算与边界')])+'</details></section>')
parts.append('<section id="download"><h2>10. 全部结果下载</h2><p><a class="button" href="Austin_unbounded_report.pdf">中文PDF报告</a><a class="button" href="Austin_unbounded_all_results.xlsx">全部结果Excel</a><a class="button" href="Austin_unbounded_complete_results_20261006.zip">完整结果ZIP</a></p>')
parts.append(paras([
    f'Excel共{len(sheets)}张表：264评估、17策略、374项指标汇总、18条超过24h的完整轨迹、新旧逐项对照、5221次派工、5485个事件、156816条区域结果、24构表、28996条分数、504次SA、90项配对统计、验收与恢复证据等。',
    'ZIP含全部288份原始任务JSON、原始分析/计划/检查点、冻结的运行代码和配置、原验收与终审、两次原生恢复的失败/重算/确认证据、报告/表格/CSV/图表。区域逐事件完整电力、可达性和AC物理字段在原始任务JSON中。',
    '6754个新版缓存及全部旧结果/缓存留在项目原位置，没有删除或跨指纹导入；大体量原始模型与全缓存不重复打包。因此这是全部任务结果包，不是可脱离原环境直接重跑的完整环境镜像。',
    '旧版桌面报告文件夹仍保留，日期相同但对应限窗版本；本版使用独立“无观察窗完整结果报告”文件夹。ZIP内SHA256SUMS.txt校验每个成员，verification.json保存逐原结果哈希与检查计数。',
]))
parts.append('<p>原运行：<code>'+str(RUN)+'</code><br>数值指纹：<code>'+fingerprint+'</code></p><details><summary>全部CSV文件</summary><ul>'+''.join('<li><a href="data/'+p.name+'">'+p.name+'</a></li>' for p in sorted(DATA.glob('*.csv')))+'</ul></details></section>')
style='''*{box-sizing:border-box}body{margin:0;background:#f1f5f7;color:#243642;font:16px/1.8 "Microsoft YaHei",sans-serif}main{max-width:1400px;margin:auto;padding:32px}header{background:#173b4b;color:white;padding:42px;border-radius:12px}h1{font-size:36px;line-height:1.45}h2{font-size:25px;color:#173b4b}h3{color:#245c75}section{background:white;margin:24px 0;padding:30px;border-radius:10px}p{margin:14px 0}.eyebrow{color:#b8d3dd;letter-spacing:2px}.small,figcaption{font-size:13px;color:#607987}header .small{color:#bfd5df}nav{display:flex;gap:20px;flex-wrap:wrap;padding:20px 0}a{color:#1b647e}.scroll{overflow:auto}table{width:100%;border-collapse:collapse;font-size:13px;line-height:1.65;margin:16px 0}th{background:#e7f0f4;white-space:nowrap}td,th{padding:8px;text-align:left;border-bottom:1px solid #dce5e9;vertical-align:top}tbody tr:nth-child(even){background:#f7f9fa}td{overflow-wrap:anywhere}code{font-size:13px;overflow-wrap:anywhere}figure{margin:25px 0}img{max-width:100%;height:auto}details{margin:16px 0;padding:12px;border:1px solid #dce5e9;border-radius:6px}summary{cursor:pointer;font-weight:bold}.button{display:inline-block;background:#245c75;color:white;text-decoration:none;padding:10px 15px;border-radius:5px;margin:5px}input{padding:12px;font-size:16px;width:75%;border:1px solid #8ba6b2;border-radius:5px}#count{margin-left:15px}@media(max-width:800px){main{padding:12px}header,section{padding:20px}h1{font-size:26px}}@media print{body{background:white;font-size:11pt}main{padding:0}section,header{padding:18px;break-inside:auto}nav,input,.button{display:none}table{font-size:8pt}tr,figure{break-inside:avoid}.scroll{overflow:visible}}'''
script="document.getElementById('filter').addEventListener('input',function(){let q=this.value.toLowerCase(),n=0;document.querySelectorAll('table.allresults tbody tr').forEach(r=>{let ok=r.textContent.toLowerCase().includes(q);r.hidden=!ok;if(ok)n++;});document.getElementById('count').textContent=n+'条';});"
(OUT/'Austin_unbounded_report.html').write_text('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Austin无观察窗完整结果报告</title><style>'+style+'</style></head><body><main>'+''.join(parts)+'</main><script>'+script+'</script></body></html>')

with PdfPages(OUT/'Austin_unbounded_report.pdf') as pdf:
    page=0;qa=OUT/'qa';qa.mkdir(exist_ok=True)
    def newpage(title):
        global page
        page+=1;f=plt.figure(figsize=(11.69,8.27))
        f.text(.055,.94,title,size=20,weight='bold',color='#173b4b')
        f.text(.055,.035,'Austin 无观察窗核心版 · 2026-10-06 · 全部原始数值保持不变',size=9,color='#607987')
        f.text(.93,.035,str(page),size=9,color='#607987');return f
    def paragraphs(f,ps,y=.86,size=12):
        f.canvas.draw();renderer=f.canvas.get_renderer();prop=FontProperties(fname=font_path,size=size)
        max_width=f.bbox.width*.88;glyph_widths={}
        for p in ps:
            ls=[];line='';used=0.
            for c in p:
                if c not in glyph_widths:glyph_widths[c]=renderer.get_text_width_height_descent(c,prop,ismath=False)[0]
                w=glyph_widths[c]
                if line and used+w>max_width:ls.append(line);line='';used=0.
                line+=c;used+=w
            if line:ls.append(line)
            t=f.text(.06,y,'\n'.join(ls),size=size,va='top',linespacing=1.48)
            y-=t.get_window_extent(renderer).height/f.bbox.height+.024
        assert y>.052,('PDF text overflow',page,y)
    def ptable(f,body,labels,rect,widths=None,size=9,scale=1.65):
        ax=f.add_axes(rect);ax.axis('off')
        t=ax.table(cellText=body,colLabels=labels,colWidths=widths,loc='center',cellLoc='center')
        t.auto_set_font_size(False);t.set_fontsize(size);t.scale(1,scale)
        for (r,c),cell in t.get_celld().items():
            cell.set_edgecolor('#dce5e9');cell.set_facecolor('#e7f0f4' if r==0 else '#ffffff' if r%2 else '#f4f7f8')
        f.canvas.draw();renderer=f.canvas.get_renderer()
        for cell in t.get_celld().values():
            tb=cell.get_text().get_window_extent(renderer);cb=cell.get_window_extent(renderer)
            assert tb.width<=cb.width+1,('PDF table horizontal overflow',page,cell.get_text().get_text())
        return t
    def savepage(f):
        f.canvas.draw();renderer=f.canvas.get_renderer()
        for t in f.texts:
            b=t.get_window_extent(renderer)
            assert b.x0>=0 and b.x1<=f.bbox.width and b.y0>=0 and b.y1<=f.bbox.height,('PDF clipped text',page,t.get_text())
        f.savefig(qa/f'pdf_page_{page:02}.png',dpi=100);pdf.savefig(f);plt.close(f)
    for title,ps in [('全部作业已完成：288任务，未完成修复为0',overview),('范围与方法：本次核心版包含哪些比较',methods)]:
        f=newpage(title);paragraphs(f,ps,size=11.8);savepage(f)
    for stage in strategy_order:
        f=newpage('全部策略结果 · '+stage.upper())
        ss=[r for r in summary if r['stage']==stage]
        body=[[short(r['strategy']),r['n'],f"{r['triangle_area']:.3f}",f"{r['weighted_triangle_area']:.3f}",f"{r['gini_restore']:.4f}",f"{r['min_time_avg_cri']:.4f}",f"{r['p90_access_restore']:.3f}",f"{r['completion_or_censor_time']/60:.2f}"] for r in ss]
        ptable(f,body,['策略','n','损失面积↓','加权损失↓','Gini↓','最弱区CRI↑','P90/min↓','恢复/h'],[.04,.49,.92,.37],
            [.29,.04,.115,.115,.085,.095,.09,.09],size=8.7,scale=1.8)
        paragraphs(f,notes[stage]+['表中为场景均值，全部未完成数为0。损失单位模型分钟；恢复/h为全部修复所需时间均值。'],y=.45,size=10.8)
        savepage(f)
    f=newpage('超过原24小时窗口的18条轨迹：全部完成')
    body=[[r['stage'],str(r['scenario_seed'])[-2:],r['strategy'],f"{r['end_minutes']:.3f}",f"{r['end_hours']:.3f}",r['remaining_at_24h'],0,f"{r['power_travel_minutes']:.2f}"] for r in over24]
    ptable(f,body,['任务','场景','策略','完成/min','完成/h','24h剩余','最终剩余','电力出行/min'],[.05,.22,.90,.65],
        [.1,.065,.145,.15,.105,.115,.1,.22],size=9,scale=1.50)
    paragraphs(f,['只有2个独立灾害场景，分别在A的3策略与C的6策略中出现。每场景15项电力修复共300分钟，单支电力队出行消耗了其余主要时间。','这些小时是模型结果；交通和维修假设的现实校准仍未完成。'],y=.18,size=10.5);savepage(f)
    f=newpage('完整派工时间线：长时间主要花在路上')
    ax=f.add_axes([.04,.19,.92,.68]);ax.imshow(plt.imread(FIG/'completed_timeline.png'));ax.axis('off')
    paragraphs(f,[why_complete[0]],y=.16,size=10.5);savepage(f)
    f=newpage('全部策略的新旧损失对照')
    body=[[r['stage'][0].upper()+'/'+short(r['strategy']),f"{r['old_triangle_mean']:.3f}",f"{r['new_triangle_mean']:.3f}",f"{r['delta_triangle_mean']:+.3f}",f"{r['new_after_24h_triangle_mean']:.3f}",str(r['identical_sequences'])+'/'+str(r['n'])] for r in comparison_summary]
    ptable(f,body,['策略','旧平均损失','新平均损失','差额','新版24h后','顺序一致'],[.04,.26,.92,.6],
        [.37,.13,.13,.105,.14,.125],size=9,scale=1.55)
    paragraphs(f,[comparison_notes[1],comparison_notes[2]],y=.22,size=10.5);savepage(f)
    f=newpage('新版轨迹的24小时分段积分')
    ax=f.add_axes([.055,.36,.89,.50]);ax.imshow(plt.imread(FIG/'horizon_decomposition.png'));ax.axis('off')
    paragraphs(f,comparison_notes[3:],y=.31,size=11);savepage(f)
    f=newpage('264条完整恢复轨迹的耗时分布')
    ax=f.add_axes([.04,.13,.92,.73]);ax.imshow(plt.imread(FIG/'recovery_duration.png'));ax.axis('off')
    f.text(.06,.085,'A/C共用灾害和基础对照；橙色为各组均值，24h竖线仅标旧版窗口。',size=11);savepage(f)
    f=newpage('主要配对统计：综合损失面积')
    body=[[r['stage'][0].upper(),short(r['reference'])+' → '+short(r['candidate']),f"{r['mean_reduction']:+.3f}",f"[{r['mean_reduction_ci95'][0]:.3f}, {r['mean_reduction_ci95'][1]:.3f}]",f"{r['wins']}/{r['losses']}/{r['ties']}",f"{r['exact_sign_p']:.4f}"] for r in main_pair]
    ptable(f,body,['任务','参考 → 候选','平均改善','95% bootstrap区间','胜/负/平','符号p'],[.04,.24,.92,.62],
        [.05,.48,.115,.19,.085,.08],size=8.2,scale=1.6)
    paragraphs(f,['改善为正表示候选损失更低。区间单位为模型分钟；符号检验与均值差检验不同。未做多重比较校正，全部90项配对统计见Excel。','C阶段基础对照与A相同；两阶段重采样随机序列不同，可产生略不同的区间。'],y=.19,size=10.5);savepage(f)
    f=newpage('正确性验收与两次P3U原生故障恢复');paragraphs(f,correctness,size=11.2);savepage(f)
    f=newpage('科学解释边界：哪些结论还不能下');paragraphs(f,limits,size=11.2);savepage(f)
    f=newpage('指标、文件与阅读顺序')
    paragraphs(f,[
        '综合损失面积 = ∫[(1−P)+(1−R)]dt。P为相对可行健康AC的供电，R为健康/当前交通总成本比并截断到[0,1]。越低越好；本版积分到每条轨迹全部修复完成。',
        '加权损失另含重要负荷与避难所可达性，不能与普通综合损失面积混为同一个目标。最弱区时间平均CRI越高越好；恢复时间Gini越低越好，但它不是收入公平指标。',
        '先读本PDF概览，再用HTML筛选264条结果。Excel含30张工作表；data/提供30份CSV，含全部22指标、逐派工、逐事件、每区域结果、新旧比较、SA过程和全部配对统计。',
        '完整ZIP收录288份原始任务JSON、分析/计划/检查点、冻结代码和配置、验收/终审及两次原生恢复证据。逐事件区域向量和完整AC物理字段保存在JSON里，未在扁平Excel中重复展开。',
        '原模型输入、6754个新版缓存与所有旧结果/缓存继续保留在原项目位置。结果包不是全部环境镜像，离线build_report.py需要原项目依赖与路径；不启动AC或TAP-B。',
        'SHA256SUMS.txt逐项校验ZIP成员，package_verification.json记录ZIP自身哈希；verification.json记录原始结果哈希和本次重新核验项目。完整研究矩阵与现实灾害校准不在本次完成范围内。',
    ],size=11.8);savepage(f)
print('PDF generated:',page,'pages',flush=True)

readme=f'''Austin 无观察窗核心版完整结果报告
生成：{report_time}
运行：{RUN}
指纹：{fingerprint}

本版：12 worker / TAP-B 2线程 / 无观察窗 / 288任务已完成 / 264条评估全部修复。
与旧版限窗报告使用独立文件夹，不要混淆。两次P3U原生故障经恢复和独立确认通过。

1. Austin_unbounded_report.pdf：{page}页中文报告，包含全部17组结果及18条超过24h的完成时间。
2. Austin_unbounded_report.html：完整中文报告、264结果筛选、全部90配对统计及下载链接。
3. Austin_unbounded_all_results.xlsx：{len(sheets)}张表，完整22指标及逐场景、派工、事件、区域等数据。
4. data/：UTF-8 BOM CSV。嵌套对象以JSON保存。Excel超过32767字符的单元格提示参看CSV/原JSON。
5. figures/：四张PNG及矢量PDF图。qa/：PDF每页预览。
6. Austin_unbounded_complete_results_20261006.zip：以上文件及288份原始任务JSON、原分析、计划、
   检查点、代码/配置、终审/物理验收、两次原生恢复的失败和独立确认记录。

ZIP original_run/results/保留完整逐事件区域电力/可达性向量与AC物理字段。
6754个新版缓存、旧缓存/结果和大体量原模型输入留在原项目，未删除或重复归档。
这是一份全部任务结果包，不是完整运行环境镜像，也不是完整研究矩阵结果。

old_vs_new_*为新旧版本比较；new_trajectory_24h_decomposition是本版轨迹事后分段积分。
旧/新TAP-B线程数不同；新旧差异不能全部归于去掉观察窗。
SHA256SUMS.txt位于ZIP内，校验每个成员；package_verification.json校验ZIP本身。
verification.json记录原始结果与前版结果SHA256、22指标重算及原物理门槛检查。

build_report.py可在原项目中用.venv-runtime/bin/python重新离线导出，不启动求解器。
程序验收通过不代表真实灾害、健康降载、交通基准和修复时间假设已获现实校准。
'''
(OUT/'README.txt').write_text(readme)
for p,d in docs:assert sha(p)==original_hashes[str(p.relative_to(RUN))]
for name,h in old_hashes.items():assert sha(OLD/name)==h
archive=OUT/'Austin_unbounded_complete_results_20261006.zip';members={}
for p in OUT.rglob('*'):
    if p.is_file() and p.name not in [archive.name,'package_verification.json','SHA256SUMS.txt'] and '__pycache__' not in p.parts:
        members[str(p.relative_to(OUT))]=p
for directory in ['results','analysis','checkpoints','plans','control']:
    for p in (RUN/directory).rglob('*'):
        if p.is_file() and p.suffix not in ['.lock','.tmp']:members['original_run/'+str(p.relative_to(RUN))]=p
for name in ['tables.json','validation.json','run_manifest.json','RUN_SCOPE.txt']:
    p=RUN/name
    if p.is_file():members['original_run/'+name]=p
for directory in native_directories:
    for p in directory.iterdir():
        if p.is_file() and p.suffix in ['.json','.log']:members['original_run/'+str(p.relative_to(RUN))]=p
for p in (RUN/'snapshot/runtime').rglob('*.py'):members['code_snapshot/runtime/'+str(p.relative_to(RUN/'snapshot/runtime'))]=p
for p in (RUN/'snapshot/runtime').rglob('*.toml'):members['code_snapshot/runtime/'+str(p.relative_to(RUN/'snapshot/runtime'))]=p
for name in ['runtime/prepared/catalog.json','data/processed/road/links.csv','data/processed/road/nodes.csv']:
    p=RUN/'snapshot'/name
    if p.is_file():members['code_snapshot/'+name]=p
for p in sorted((ROOT/'reports').glob('core_unbounded_*20261006.*')):
    if p.suffix in ['.json','.log']:members['provenance/'+p.name]=p
for name in ['runtime/ASSUMPTIONS.md','runtime/configs/core-unbounded-12.toml','runtime/output/preflight-unbounded-12-20261006/summary.json']:
    p=ROOT/name;members['provenance/'+('preflight_' if p.name=='summary.json' else '')+p.name]=p
for name in ['analysis/aggregate.csv','analysis/scenario_rows.csv','analysis/paired_statistics.json','control/status.json','control/launch.json','run_manifest.json']:
    members['previous_24h_run/'+name]=OLD/name
hashes={name:sha(p) for name,p in sorted(members.items())}
sums=''.join(f'{h}  {name}\n' for name,h in hashes.items())
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name,p in sorted(members.items()):z.write(p,name)
    z.writestr('SHA256SUMS.txt',sums)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for name,h in hashes.items():assert hashlib.sha256(z.read(name)).hexdigest()==h
assert sum(name.startswith('original_run/results/') and name.endswith('.json') for name in members)==288
# Validate every relative link in the HTML against concrete delivered files.
from html.parser import HTMLParser
class LinkCheck(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ['href','src'] and not v.startswith(('#','data:','https:','http:')):
                assert (OUT/v).is_file(),('Missing HTML link',v)
LinkCheck().feed((OUT/'Austin_unbounded_report.html').read_text())
dump(OUT/'package_verification.json',dict(archive=archive.name,sha256=sha(archive),bytes=archive.stat().st_size,
    members=len(members)+1,original_result_jsons=288,csv_files=len(list(DATA.glob('*.csv'))),excel_sheets=len(sheets),pdf_pages=page,
    all_archive_member_hashes_verified=True,workbook_xml_verified=True,html_local_links_verified=True,
    original_results_unchanged=True,previous_results_unchanged=True,solvers_started=0,
    verified_at=datetime.now(ZoneInfo('America/Chicago')).isoformat()))
print('DONE',archive,archive.stat().st_size,'bytes',len(members)+1,'members',flush=True)
