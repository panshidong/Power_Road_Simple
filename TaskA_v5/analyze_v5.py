from pathlib import Path
import json,csv,collections,statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent;OLD=Path('/home/workenv/TaskA_v4')
OUT=ROOT/'analysis';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
COLORS={'CEN':'#687582','JSH':'#2672a4','IJSH':'#c96b21'}
LABELS={'CEN':'Structural','JSH':'Service Shapley','IJSH':'Access-augmented Shapley'}

def read_rows(folder):
    rows=collections.defaultdict(dict)
    for p in folder.glob('worker_*/rows_*.csv'):
        for r in csv.DictReader(p.open()):
            variant=r.get('variant','main');rows[(variant,r['scenario_id'])][r['strategy_id']]=r
    return rows

def paired_percent(ref,cand,seed):
    rng=np.random.default_rng(seed);results=[];area=[]
    for _ in range(20):
        ix=rng.integers(0,len(ref),size=(1000,len(ref)))
        rm=ref[ix].mean(axis=1);cm=cand[ix].mean(axis=1)
        results.extend(100*(rm-cm)/rm);area.extend(rm-cm)
    return np.quantile(results,[.025,.975]).tolist(),np.quantile(area,[.025,.975]).tolist()

rows=read_rows(OLD/'results/evaluation');rows.update(read_rows(OLD/'results/robustness'))
records=[]
for vindex,variant in enumerate(['main','small','large','light','clustered']):
    chosen=sorted(k for k in rows if k[0]==variant)
    for h,(ref,cand) in enumerate([('CEN','JSH'),('JSH','IJSH')],1):
        a=np.array([float(rows[k][ref]['triangle_area']) for k in chosen]);b=np.array([float(rows[k][cand]['triangle_area']) for k in chosen]);d=a-b
        ci,aci=paired_percent(a,b,2026090600+10*vindex+h)
        records.append(dict(variant=variant,hypothesis=h,reference=ref,candidate=cand,n=len(a),mean_reduction=float(d.mean()),percent_reduction=float(100*d.mean()/a.mean()),percent_ci=ci,area_ci=aci,wins=int((d>1e-12).sum()),losses=int((d< -1e-12).sum()),ties=int((abs(d)<=1e-12).sum())))
(OUT/'corrected_percentage_intervals.json').write_text(json.dumps(records,indent=2))
# Figures use descriptive method names rather than manuscript-specific acronyms.
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig,ax=plt.subplots(figsize=(10,2.4),layout='constrained');ax.set_xlim(0,12);ax.set_ylim(0,2.8);ax.axis('off')
boxes=['Sample plausible\ndamage cases','Calculate component\nShapley values','Average scores\nacross damage cases','Prepare standing\npriority tables','Identify damaged\ncomponents','Dispatch eligible\nrepair crews']
for i,label in enumerate(boxes):
    x=.1+2*i;ax.add_patch(FancyBboxPatch((x,.7),1.7,1.25,boxstyle='round,pad=.03',fc='#e0edf5' if i<4 else '#fae7d6',ec='#577080'));ax.text(x+.85,1.325,label,ha='center',va='center',fontsize=8)
    if i<5:ax.add_patch(FancyArrowPatch((x+1.72,1.325),(x+1.95,1.325),arrowstyle='-|>',mutation_scale=11,color='#577080'))
ax.text(4,2.35,'Before a disruption',ha='center',weight='bold');ax.text(10,2.35,'After a disruption',ha='center',weight='bold')
fig.savefig(OUT/'figure_1_workflow.png',dpi=240);plt.close(fig)
conv=list(csv.DictReader((OLD/'analysis/table_convergence.csv').open()))
fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
for strategy in ['JSH','IJSH']:
    for trade,style in [('power','-'),('road','--')]:
        subset=sorted([r for r in conv if r['strategy']==strategy and r['asset_type']==trade],key=lambda r:int(r['checkpoint']))
        for ax,metric,title in zip(axes,['spearman_mean','top3_overlap'],['Agreement in repair order','Overlap among first three priorities']):
            ax.plot([int(r['checkpoint']) for r in subset],[float(r[metric]) for r in subset],marker='o',ls=style,color=COLORS[strategy],label=LABELS[strategy]+' / '+trade)
            ax.set_ylim(.75,1.01);ax.set_xlabel('Damage cases used to prepare priorities');ax.set_title(title,fontsize=10);ax.grid(alpha=.2)
axes[0].set_ylabel('Mean Spearman rank correlation');axes[1].set_ylabel('Mean overlap fraction');axes[1].legend(fontsize=6.6,loc='lower right')
fig.savefig(OUT/'figure_2_priority_stability.png',dpi=240);plt.close(fig)
main_ids=sorted(k for k in rows if k[0]=='main');data=[np.array([float(rows[k][s]['triangle_area']) for k in main_ids]) for s in ['CEN','JSH','IJSH']]
fig,ax=plt.subplots(figsize=(7.5,4.2),layout='constrained');bp=ax.boxplot(data,patch_artist=True,showmeans=True,meanprops={'marker':'D','markerfacecolor':'white','markeredgecolor':'black'},tick_labels=[LABELS[s] for s in ['CEN','JSH','IJSH']])
for patch,color in zip(bp['boxes'],COLORS.values()):patch.set_facecolor(color);patch.set_alpha(.7)
ax.set_yscale('log');ax.set_ylabel('System resilience loss (model time units)');ax.grid(axis='y',alpha=.2)
fig.savefig(OUT/'figure_3_loss_distributions.png',dpi=240);plt.close(fig)
original_effects=list(csv.DictReader((OLD/'analysis/paired_effects.csv').open()))
fig,axes=plt.subplots(1,2,figsize=(9,4.2),sharey=True,layout='constrained');rng=np.random.default_rng(121)
for i,(ax,ref,cand,title) in enumerate(zip(axes,data[:2],data[1:],['Service Shapley vs structural','Access-augmented vs service Shapley'])):
    diff=ref-cand;rec=[r for r in original_effects if r['variant']=='main'][i];lo=float(rec['mean_reduction_ci95_lower']);hi=float(rec['mean_reduction_ci95_upper']);m=diff.mean()
    ax.scatter(rng.normal(0,.07,len(diff)),diff,s=12,alpha=.35,color=list(COLORS.values())[i+1]);ax.errorbar(0,m,yerr=[[m-lo],[hi-m]],fmt='D',mec='black',color='black',mfc=list(COLORS.values())[i+1],capsize=4);ax.scatter(.3,np.median(diff),marker='s',facecolors='white',edgecolors='black')
    ax.set_yscale('symlog',linthresh=10);ax.set_xticks([0,.3],['Mean and 95% interval','Median']);ax.tick_params(axis='x',labelsize=8);ax.axhline(0,ls='--',color='#666666');ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.2)
axes[0].set_ylabel('Paired reduction in system loss (positive favors candidate)')
fig.savefig(OUT/'figure_4_paired_effects.png',dpi=240);plt.close(fig)
fig,ax=plt.subplots(figsize=(8.5,4.2),layout='constrained')
variants=['main','small','large','light','clustered']
for h,color,label in [(1,COLORS['JSH'],'Service Shapley vs structural'),(2,COLORS['IJSH'],'Access-augmented vs service Shapley')]:
    selected=[r for r in records if r['hypothesis']==h]
    y=np.arange(5)+(-.11 if h==1 else .11);x=np.array([r['percent_reduction'] for r in selected]);ci=np.array([r['percent_ci'] for r in selected])
    ax.errorbar(x,y,xerr=np.stack([x-ci[:,0],ci[:,1]-x]),fmt='o',capsize=3,color=color,label=label)
ax.set_yticks(range(5),['Main damage assumptions','Smaller disruptions','Larger disruptions','Lighter road damage','Spatially clustered damage']);ax.invert_yaxis();ax.axvline(0,color='#666666',ls='--',lw=1)
ax.set_xlabel('Reduction in mean system loss (%) with paired bootstrap 95% interval');ax.grid(axis='x',alpha=.2);ax.legend(fontsize=8,loc='upper right')
fig.savefig(OUT/'figure_5_damage_patterns.png',dpi=240);plt.close(fig)

selection=json.loads((ROOT/'case_selection.json').read_text());case_records=[]
fig,axes=plt.subplots(3,2,figsize=(9.2,8),layout='constrained',sharey=True)
for row,(label,case) in enumerate(selection['cases'].items()):
    for strategy in ['CEN','JSH','IJSH']:
        p=ROOT/'case_results'/f'{case}_{strategy}_penalty_9999.json'
        if not p.exists():raise RuntimeError(f'Missing case: {p}')
        r=json.loads(p.read_text());events=r['event_log'];t=[e['time'] for e in events]
        axes[row,0].step(t,[e['power_func'] for e in events],where='post',color=COLORS[strategy],label=LABELS[strategy],lw=1.5)
        axes[row,1].step(t,[e['state']['road_func'] for e in events],where='post',color=COLORS[strategy],lw=1.5)
        case_records.append({k:r[k] for k in ['case','label','strategy','triangle_area','power_loss','road_loss','power_sequence','road_sequence']})
    for col,title in enumerate(['Power service','Road service']):
        axes[row,col].set_title(f"{case.replace('scenario_','Scenario ')}: {label.replace('_',' ')} — {title}",fontsize=10)
        axes[row,col].set_ylim(-.025,1.035);axes[row,col].grid(alpha=.2)
        axes[row,col].set_xlabel('Elapsed model time');axes[row,col].set_ylabel('Functionality')
axes[0,0].legend(fontsize=7,loc='lower right')
fig.savefig(OUT/'figure_6_recovery_cases.png',dpi=240);plt.close(fig)
(OUT/'case_summary.json').write_text(json.dumps(case_records,indent=2))

# Supplementary case sensitivity; the selection deliberately spans outcome extremes.
case_sensitivity=[]
for penalty in [1000,9999,100000]:
    for label,case in selection['cases'].items():
        outcomes={s:json.loads((ROOT/'case_results'/f'{case}_{s}_penalty_{penalty}.json').read_text())['triangle_area'] for s in ['CEN','JSH','IJSH']}
        case_sensitivity.append(dict(penalty=penalty,case=case,label=label,**outcomes,H1=outcomes['CEN']-outcomes['JSH'],H2=outcomes['JSH']-outcomes['IJSH']))
(OUT/'case_closure_sensitivity.json').write_text(json.dumps(case_sensitivity,indent=2))

closure_files=list((ROOT/'closure_check').glob('scenario_*.json'))
if len(closure_files)==300:
    cr={json.loads(p.read_text())[0]['scenario_id']:{r['strategy_id']:r for r in json.loads(p.read_text())} for p in closure_files}
    ids=sorted(cr);closure=[]
    for h,(ref,cand) in enumerate([('CEN','JSH'),('JSH','IJSH')],1):
        a=np.array([cr[k][ref]['triangle_area'] for k in ids]);b=np.array([cr[k][cand]['triangle_area'] for k in ids]);d=a-b
        pci,aci=paired_percent(a,b,2026090670+h)
        closure.append(dict(hypothesis=h,n=300,reference_mean=float(a.mean()),candidate_mean=float(b.mean()),mean_reduction=float(d.mean()),percent_reduction=float(100*d.mean()/a.mean()),percent_ci=pci,area_ci=aci,wins=int((d>1e-12).sum()),losses=int((d< -1e-12).sum()),ties=int((abs(d)<=1e-12).sum())))
    (OUT/'closure_check_summary.json').write_text(json.dumps(closure,indent=2))
    print('CLOSURE',json.dumps(closure))
print('PERCENTAGE INTERVALS',json.dumps(records[:2]))
print('CASE SUMMARY',json.dumps(case_records))
