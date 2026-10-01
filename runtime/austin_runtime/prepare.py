"""Explicit preparation on the destination host; never executed on import."""
from __future__ import annotations
import csv, gzip, io, math, re, shutil
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile
from .common import AUSTIN, atomic_json, digest, file_lock, local_path, read_json, rows, sha_file, write_rows


def commands(text):
    current = ""
    for line in text.splitlines():
        line = line.split("!",1)[0].strip()
        if not line: continue
        if line.startswith("~"): current += " "+line[1:].strip()
        else:
            if current: yield current
            current = line
    if current: yield current


def sphere(points):
    import numpy as np
    a = np.radians(np.asarray(points,dtype=float))
    return np.column_stack([np.cos(a[:,1])*np.cos(a[:,0]),np.cos(a[:,1])*np.sin(a[:,0]),np.sin(a[:,1])])


def prepare(cfg):
    import networkx as nx
    import numpy as np
    from scipy.spatial import cKDTree
    dst = local_path(cfg,"prepared"); dst.mkdir(parents=True,exist_ok=True)
    source = AUSTIN/"data/processed"
    with file_lock(dst/"prepare.lock",blocking=False):
        locked=read_json(AUSTIN/"sources.lock.json")
        for item in locked["sources"]:
            path=AUSTIN/item["path"]
            if not path.exists(): raise FileNotFoundError(f"{path}; run the documented fetch/build steps first")
            if sha_file(path)!=item["sha256"]: raise ValueError(f"Source hash mismatch: {path}")
        nodes={int(r["node_id"]):r for r in rows(source/"road/nodes.csv")}
        links=list(rows(source/"road/links.csv"))
        graph=nx.MultiDiGraph(); physical=[]; roads=defaultdict(list)
        for r in links:
            if r["centroid_connector"]=="1": continue
            u,v,i=int(r["from_node"]),int(r["to_node"]),int(r["link_id"])
            graph.add_edge(u,v,key=i,weight=float(r["free_flow_time_source"]))
            roads[f"road:{min(u,v)}-{max(u,v)}"].append(i)
        main_component=max(nx.strongly_connected_components(graph),key=len)
        physical=sorted(main_component)
        pts=[(nodes[i]["longitude"],nodes[i]["latitude"]) for i in physical]
        tree=cKDTree(sphere(pts))
        def nearest(point):
            chord,index=tree.query(sphere([point])[0]); node=physical[int(index)]
            return node, float(2*6371008.8*np.arcsin(min(1,float(chord)/2)))
        depot,_=nearest(cfg["recovery"]["depot_lon_lat"])
        shelter,_=nearest(cfg["recovery"]["shelter_lon_lat"])
        subs=list(rows(source/"power/distribution_v03/substations.csv"))
        feeders=list(rows(source/"power/distribution_v03/feeders.csv"))
        sub_region={r["substation_id"]:r["region"] for r in feeders}
        assets={}
        for key,ids in sorted(roads.items()):
            first=links[ids[0]-1]; u,v=int(first["from_node"]),int(first["to_node"])
            if not ({u,v}&main_component): continue
            assets[key]=dict(kind="road",link_ids=ids,targets=[u,v],offroad_m=0,
                lon=(float(nodes[u]["longitude"])+float(nodes[v]["longitude"]))/2,
                lat=(float(nodes[u]["latitude"])+float(nodes[v]["latitude"]))/2)
        for sub in subs:
            sid=sub["Name"]; node,distance=nearest([sub["Longitude"],sub["Latitude"]])
            assets["power:"+sid]=dict(kind="power",substation=sid,region=sub_region[sid],targets=[node],
                offroad_m=distance,lon=float(sub["Longitude"]),lat=float(sub["Latitude"]),
                nominal_kw=float(sub["Total Real Load"])/1000)
        # Materialize the complete original regional networks, including cross-feeder controls.
        with ZipFile(AUSTIN/"data/raw/power/syn-Austin-TDgrid-v03.zip") as z:
            for region in sorted(set(sub_region.values())):
                prefix=f"syn-austin-D_only-v03/{region}/base/opendss/"
                for name in z.namelist():
                    if not name.startswith(prefix) or not name.lower().endswith(".dss"): continue
                    rel=Path(name[len(prefix):])
                    if rel.is_absolute() or ".." in rel.parts: raise ValueError("Unsafe archive member")
                    target=dst/"models"/region/rel; target.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(name) as inp, target.open("wb") as out: shutil.copyfileobj(inp,out)
                master=dst/"models"/region/"Master.dss"
                # Keep native setup/voltage bases; remove the automatic load-flow command.
                # Solve only after runtime fault isolation, loading and DataPath have been set.
                text="\n".join(line for line in master.read_text().splitlines() if not re.match(r"(?i)^\s*solve\b",line))
                (master.parent/"Master_runtime.dss").write_text(text+"\n")
        for key,a in assets.items():
            if a["kind"]!="power": continue
            folder=dst/"models"/a["region"]/a["substation"]
            transformers=[]
            for filename in ("Transformers.dss","Regulators.dss"):
                path=folder/filename
                if not path.exists(): continue
                for cmd in commands(path.read_text()):
                    m=re.match(r"(?i)new\s+(transformer\.\S+)",cmd)
                    if not m: continue
                    kv=[float(x) for x in re.findall(r"(?i)\bkv\s*=\s*([\d.eE+-]+)",cmd)]
                    array=re.search(r"(?i)\bkvs\s*=\s*\[([^]]+)\]",cmd)
                    if array: kv += [float(x) for x in re.findall(r"[\d.eE+-]+",array.group(1))]
                    if kv and max(kv)>20 and min(kv)<20: transformers.append(m.group(1).lower())
            if not transformers: raise ValueError(f"Cannot identify HV/MV transformer for {key}; inspect {folder}")
            a["isolation_transformers"]=sorted(set(transformers))
        zone_ids=sorted(i for i,r in nodes.items() if r["node_type"]=="centroid")
        ztree=cKDTree(sphere([(nodes[i]["longitude"],nodes[i]["latitude"]) for i in zone_ids]))
        zone_kw=defaultdict(float); sub_kw=defaultdict(float); load_count=0; unique=set()
        catalog_path=dst/"load_catalog.csv.gz"
        fields=["load_id","feeder_id","substation_id","region","zone","bus_id","kw","kvar"]
        with catalog_path.open("wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as compressed, io.TextIOWrapper(compressed,encoding="utf-8",newline="") as handle:
            w=csv.DictWriter(handle,fieldnames=fields);w.writeheader()
            batch=[]
            def flush():
                nonlocal load_count
                if not batch:return
                _,indices=ztree.query(sphere([(r["longitude"],r["latitude"]) for r in batch]))
                for r,index in zip(batch,indices):
                    region=sub_region[r["substation_id"]]; identity=(region,r["load_id"].lower())
                    if identity in unique: raise ValueError(f"Duplicate regional DSS load identity: {identity}")
                    unique.add(identity); zone=zone_ids[int(index)]; kw=float(r["kw"])
                    zone_kw[zone]+=kw; sub_kw[r["substation_id"]]+=kw
                    w.writerow(dict(load_id=r["load_id"].lower(),feeder_id=r["feeder_id"],substation_id=r["substation_id"],
                        region=region,zone=zone,bus_id=r["bus_id"],kw=kw,kvar=r["kvar"]))
                    load_count+=1
                batch.clear()
            for r in rows(source/"power/distribution_v03/loads.csv.gz"):
                batch.append(r)
                if len(batch)>=8192:flush()
            flush()
        # Sampled road betweenness is explicit and reproducible; 0 selects exact calculation.
        k=cfg["traffic"]["centrality_sources"]
        bc=nx.edge_betweenness_centrality(graph,k=min(k,len(graph)) if k else None,normalized=True,
                                         weight="weight",seed=cfg["traffic"]["centrality_seed"])
        cen={key:0.0 for key in assets}
        arc_asset={i:key for key,a in assets.items() if a["kind"]=="road" for i in a["link_ids"]}
        for (_,_,i),score in bc.items():
            if i in arc_asset: cen[arc_asset[i]]+=score
        total=sum(sub_kw.values())
        for sid,value in sub_kw.items():cen["power:"+sid]=value/total
        critical=[]
        for lon,lat in cfg["recovery"]["critical_sites_lon_lat"]:
            choices=[sid for sid in sub_kw if sid not in critical]
            if not choices:raise ValueError("More critical facilities than available substations")
            critical.append(min(choices,key=lambda sid:((assets["power:"+sid]["lon"]-lon)*math.cos(math.radians(lat)))**2+(assets["power:"+sid]["lat"]-lat)**2))
        signals=list(rows(source/"coupling/signal_road_power_provisional.csv"))
        signal_ids={r["signal_id"] for r in signals}
        approaches=[r for r in rows(source/"coupling/signal_approach_links.csv") if r["signal_id"] in signal_ids]
        inputs={str(p.relative_to(AUSTIN)):sha_file(p) for p in sorted(source.rglob("*")) if p.is_file()}
        inputs.update({r["path"]:r["sha256"] for r in locked["sources"]})
        inputs["runtime/prepared/load_catalog.csv.gz"]=sha_file(catalog_path)
        model_hashes={str(p.relative_to(dst)):sha_file(p) for p in sorted((dst/"models").rglob("*.dss"))}
        payload=dict(schema=1,input_hashes=inputs,derived_model_hashes=model_hashes,assets=assets,centrality=cen,regions=sorted(set(sub_region.values())),
            depot=depot,shelter=shelter,critical_substations=critical,zone_nominal_kw=dict(zone_kw),
            substation_nominal_kw=dict(sub_kw),substation_region=sub_region,signals=signals,approaches=approaches,
            load_count=load_count,road_nodes=len(nodes),road_links=len(links),
            provenance=dict(centroid_through_traffic=False,parallel_arcs_retained=True,datum_declared=False,
                power_scope="six original regional AC networks with ideal 230kV boundary sources; not full transmission co-simulation",
                zone_mapping="nearest TAZ centroid, not a polygon overlay",power_damage_unit="HV/MV substation transformer group",
                road_damage_unit="all parallel/directional arcs sharing an unordered endpoint pair",
                damage_road_groups_excluded_outside_main_scc=sum(key not in assets for key in roads),
                power_access="nearest physical node in main SCC; last-mile distance retained, actual driveway unverified",
                centrality_sample_sources=k,prepared_config=digest(dict(traffic=cfg["traffic"],recovery=cfg["recovery"])) ))
        atomic_json(dst/"catalog.json",payload)
        return {"catalog":str(dst/"catalog.json"),"power_assets":len(subs),"road_assets":sum(a["kind"]=="road" for a in assets.values()),"loads":load_count}
