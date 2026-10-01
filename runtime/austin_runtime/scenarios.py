from __future__ import annotations
import math, random
from .common import digest

VARIANTS={
    "main":([8,15],[3,13],[.5,1.0],False,20264001),
    "small":([3,7],[1,5],[.5,1.0],False,20265001),
    "large":([15,20],[12,18],[.5,1.0],False,20266001),
    "light":([8,15],[3,13],[.2,.6],False,20267001),
    "clustered":([8,15],[3,13],[.5,1.0],True,20268001),
}


def generate(catalog,cfg,variant,index,seed0=None):
    pr,rr,drops,clustered,seed=VARIANTS[variant]
    seed=(seed0 if seed0 is not None else seed)+index
    rng=random.Random(seed);npower=rng.randint(*pr);nroad=rng.randint(*rr)
    power=sorted(k for k,a in catalog["assets"].items() if a["kind"]=="power")
    roads=sorted(k for k,a in catalog["assets"].items() if a["kind"]=="road")
    if npower>len(power) or nroad>len(roads):raise ValueError("Insufficient damage candidates")
    if clustered:
        centre=catalog["assets"][rng.choice(power)]
        def distance(key):
            a=catalog["assets"][key]
            return ((a["lon"]-centre["lon"])*math.cos(math.radians(centre["lat"])))**2+(a["lat"]-centre["lat"])**2
        # A shared expanding geographical radius supplies both trades.
        p=sorted(power,key=lambda k:(distance(k),k));r=sorted(roads,key=lambda k:(distance(k),k))
        radius=max(distance(p[npower-1]),distance(r[nroad-1]))*1.5+1e-12
        power=[k for k in p if distance(k)<=radius];roads=[k for k in r if distance(k)<=radius]
    selected_power=rng.sample(power,npower); selected_roads=rng.sample(roads,nroad)
    damage={}
    for key in selected_roads:
        remaining=1-rng.uniform(*drops)
        # Finite closure probability; light ensemble explicitly has no closures.
        damage[key]=0.0 if remaining<.1 else remaining
    for key in selected_power:
        local=random.Random(int(digest([seed,key,"partial-transformer-rating"])[:16],16))
        damage[key]=local.uniform(*cfg["power"]["partial_rating_range"]) if local.random()<cfg["power"]["partial_failure_probability"] else 0.0
    return dict(id=f"{variant}_{index+1:05d}",variant=variant,index=index,seed=seed,damage=dict(sorted(damage.items())))
