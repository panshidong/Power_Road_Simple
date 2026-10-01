"""Publication-exportable plots built only from completed native experiment records."""
from __future__ import annotations
from collections import defaultdict
import numpy as np
from .common import digest, read_json


def supplemental_figures(output,dest,rows,stability,fp):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder=dest/"figures";folder.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout()
        for suffix in ("png","svg"):fig.savefig(folder/(name+"."+suffix),dpi=180)
        plt.close(fig)
    facets=defaultdict(list)
    for row in rows:
        if row["stage"] in ("a-shift","a-alpha","c-k","c-shift","c-alternatives","c-decomp","b-sensitivity"):
            facets[(row["stage"],row["variant"],row["tag"],row["k"],row["options_digest"])].append(row)
    for key,values in sorted(facets.items(),key=lambda item:str(item[0])):
        names=sorted({r["strategy"] for r in values})
        fig,axes=plt.subplots(1,2,figsize=(max(10,len(names)*2),4))
        for ax,metric in zip(axes,("triangle_area","weighted_triangle_area")):
            ax.bar(names,[np.mean([r[metric] for r in values if r["strategy"]==name]) for name in names])
            ax.set_ylabel(metric+" (mean)");ax.tick_params(axis="x",rotation=35)
        fig.suptitle(" / ".join(map(str,key)))
        save(fig,key[0]+"_"+digest(key)[:12])
    if stability:
        fig,axes=plt.subplots(1,2,figsize=(11,4))
        for ax,kind in zip(axes,("power:","road:")):
            for strategy in sorted({r["strategy"] for r in stability}):
                subset=sorted((r for r in stability if r["asset_type"]==kind and r["strategy"]==strategy and r["spearman"] is not None),key=lambda r:r["checkpoint"])
                if subset:ax.plot([r["checkpoint"] for r in subset],[r["spearman"] for r in subset],marker="o",label=strategy)
            ax.set(xlabel="Construction scenarios",ylabel="Spearman vs final table",title=kind)
            ax.legend(fontsize=7)
        save(fig,"task_a_table_stability")
    # Read one trajectory at a time. No whole-study event log is held in memory.
    cases=defaultdict(list)
    for path in sorted((output/"results/a-cases").glob("*.json")):
        if path.name.endswith(".error.json"):continue
        record=read_json(path)
        if record["fingerprint"]!=fp:raise ValueError("Mixed case fingerprints")
        job=record["job"];cases[(job["index"],job["options"]["case_label"],job["options"]["penalty"])].append(path)
    for key,paths in cases.items():
        fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
        for path in paths:
            record=read_json(path);events=record["payload"]["events"]
            for ax,metric in zip(axes,("power_func","road_func")):
                ax.step([e["time"] for e in events],[e[metric] for e in events],where="post",label=record["job"]["strategy"])
                ax.set_ylabel(metric);ax.legend()
        axes[-1].set_xlabel("Time (assumed minutes)");fig.suptitle(f"Case {key}")
        save(fig,"task_a_case_"+digest(key)[:12])
    manifest=read_json(output/"run_manifest.json")
    representative=manifest["config"]["task_b"]["representative_scenario"]-1
    eq=manifest["config"]["equity"]
    for path in sorted((output/"results/b-main").glob("*.json")):
        if path.name.endswith(".error.json"):continue
        record=read_json(path);job=record["job"]
        if job["index"]!=representative:continue
        events=record["payload"]["events"];times=[e["time"] for e in events]
        if len(times)<2 or times[-1]<=0:continue
        cri=np.asarray([np.asarray(e["electric"])*eq["electric_weight"]+np.asarray(e["access"])*eq["access_weight"] for e in events])
        fig,ax=plt.subplots(figsize=(10,6))
        mesh=ax.pcolormesh(times,np.arange(cri.shape[1]+1),cri[:-1].T,vmin=0,vmax=1,cmap="viridis",shading="flat")
        fig.colorbar(mesh,ax=ax,label="CRI")
        ax.set(xlabel="Time (assumed minutes)",ylabel="Zones ordered by TAZ ID",title=job["strategy"])
        save(fig,"task_b_cri_"+job["strategy"])
