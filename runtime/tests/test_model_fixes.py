"""2026-10-07 model fixes: crews, repair times, crew routing, OD order, antithetic
permutations and local-shedding topology. No native solver is started."""
from __future__ import annotations
import copy, unittest
from collections import defaultdict
from pathlib import Path
from austin_runtime.common import RUNTIME, config
from austin_runtime.power import downstream_loads
from austin_runtime.rankings import construct, critical_set, order
from austin_runtime.simulation import Coupled, crew_counts, repair_minutes
from austin_runtime.traffic import TrafficState


V2=RUNTIME/"configs/core-v2.toml"


class CrewsAndRepairs(unittest.TestCase):
    def test_crews_scale_with_damaged_units(self):
        cfg=config(V2)
        damage={**{f"power:p{i}":0.0 for i in range(8)},**{f"road:{i}-{i+1}":.3 for i in range(3)}}
        self.assertEqual(crew_counts(cfg,damage),{"power":3,"road":1})
        damage.update({f"power:q{i}":.5 for i in range(30)})
        self.assertEqual(crew_counts(cfg,damage)["power"],cfg["recovery"]["max_crews_per_trade"])
        self.assertEqual(crew_counts(cfg,{"power:p":0.0})["road"],0)
        self.assertEqual(crew_counts(config(),damage),{"power":1,"road":1}) # research.toml keeps the 2026-10-06 rules

    def test_repair_time_follows_severity(self):
        cfg=config(V2);r=cfg["recovery"]
        self.assertEqual(repair_minutes(cfg,"power",0.0),r["power_repair_minutes_full"])
        least,most=r["power_repair_minutes_partial"]
        self.assertAlmostEqual(repair_minutes(cfg,"power",0.8),least+(most-least)*0.2)
        self.assertAlmostEqual(repair_minutes(cfg,"road",0.25),r["road_repair_minutes_partial"][0]+
                               (r["road_repair_minutes_partial"][1]-r["road_repair_minutes_partial"][0])*0.75)
        self.assertLess(repair_minutes(cfg,"power",0.8),repair_minutes(cfg,"power",0.2))
        self.assertEqual(repair_minutes(config(),"power",0.0),20.0)

    def test_simulation_dispatches_parallel_crews_with_severity_durations(self):
        cfg=config(V2)
        damage={"power:a":0.0,"power:b":0.5,"power:c":0.0,"power:d":0.7,"road:1-2":0.0}
        engine=Coupled.__new__(Coupled)
        engine.cfg=cfg;engine.catalog={"assets":{k:{"kind":k.split(":")[0]} for k in damage},"critical_substations":[]}
        engine.depot=0;engine.zones=[1];engine.excluded_zones=[];engine.healthy_power={"served_kw":1.,"nominal_kw":1.}
        engine.state=lambda remaining,penalty=None:dict(remaining=dict(remaining),traffic=None)
        engine.site_travel=lambda traffic,origin,asset,roundtrip_lastmile=False:(30.0,7)
        def record(time,remaining,state,critical):
            f=1-len(remaining)/len(damage)
            return dict(time=time,remaining=dict(remaining),power_func=f,road_func=f,weighted_power_func=f,shelter_access=f,electric=[f],access=[f])
        engine.record=record
        run=engine.simulate({"id":"t","seed":1,"damage":damage},sorted(damage))
        self.assertEqual(run["crews"],{"power":2,"road":1})
        first=[x for x in run["dispatch"] if x["time"]==0]
        self.assertEqual(sorted(x["crew"] for x in first),["power-0","power-1","road-0"])
        for x in run["dispatch"]:
            self.assertAlmostEqual(x["finish"]-x["time"],30.0+repair_minutes(cfg,x["asset"].split(":")[0],damage[x["asset"]]))
        self.assertTrue(run["complete"])


class CrewRouting(unittest.TestCase):
    def test_crews_use_capped_costs_while_access_keeps_equilibrium_costs(self):
        links=[dict(from_node=u,to_node=v,link_id=i+1) for i,(u,v) in enumerate([(3,4),(4,5),(3,5)])]
        congested=[100.,1.,50.];capped=[6.,1.,50.]
        state=TrafficState(links,congested,[0]*3,[1]*3,{},first_thru=1,crew_costs=capped)
        self.assertEqual(state.travel(3,[5]),(50.,5))
        self.assertEqual(state.travel(3,[5],crew=True),(7.,5))
        plain=TrafficState(links,congested,[0]*3,[1]*3,{},first_thru=1)
        self.assertEqual(plain.travel(3,[5],crew=True),(50.,5))


class Rankings(unittest.TestCase):
    def test_od_ties_follow_paired_table_not_asset_name(self):
        cfg=config(V2)
        assets=["power:p1","road:1-2","road:3-4","road:5-6"]
        catalog={"centrality":{"power:p1":.1,"road:1-2":.01,"road:3-4":.5,"road:5-6":.2}}
        tables={"CEN":catalog["centrality"]}
        sequence,_=order(assets,"OD_CEN",tables,catalog,cfg,{"road:5-6":0.0})
        self.assertEqual([a for a in sequence if a.startswith("road")],["road:3-4","road:5-6","road:1-2"])
        sequence,_=order(assets,"OD_CEN",tables,catalog,cfg,{"road:1-2":2.0})
        self.assertEqual([a for a in sequence if a.startswith("road")],["road:1-2","road:3-4","road:5-6"])
        legacy=copy.deepcopy(cfg);legacy["task_c"]["od_tie_break"]="asset_id"
        sequence,_=order(assets,"OD_CEN",tables,catalog,legacy,{})
        self.assertEqual([a for a in sequence if a.startswith("road")],["road:1-2","road:3-4","road:5-6"])

    def test_zone_pairs_cover_loaded_zones(self):
        catalog={"assets":{"power:s1":{"kind":"power","targets":[2000]}},"critical_substations":["s1"],"shelter":1500,"depot":1200,
                 "centrality":{"power:s1":1.},"zone_nominal_kw":{"5":10.,"6":0.,"7":3.}}
        selection=critical_set(catalog,"default",config(V2))
        self.assertEqual(selection["pair_mode"],"zones_to_essential");self.assertEqual(selection["zones"],[5,7])
        self.assertEqual(selection["destinations"],[1500,2000])
        legacy=critical_set(catalog,"default",config())
        self.assertEqual(legacy["pair_mode"],"depot_to_critical");self.assertNotIn("zones",legacy)

    def test_antithetic_permutations_reverse_their_partner(self):
        cfg=config(V2);cfg["task_a"]["shapley_permutations"]=4
        orders=[]
        class Engine:
            catalog={"centrality":defaultdict(float)}
            def coalition(self,scenario,repaired):
                orders.append(tuple(repaired));return float(len(repaired)),0.
        scenario={"seed":11,"damage":{f"power:p{i}":0.0 for i in range(5)}}
        with self.subTest("antithetic"):
            import tempfile
            with tempfile.TemporaryDirectory() as folder:
                construct(Engine(),scenario,cfg,Path(folder)/"c.json","sig")
        # Each permutation evaluates the empty coalition, then adds one asset at a time.
        perms=[];current=[]
        for coalition in orders:
            if not coalition:
                if current:perms.append(current)
                current=[];previous=set()
            else:
                current.append(next(iter(set(coalition)-previous)));previous=set(coalition)
        perms.append(current)
        self.assertEqual(len(perms),4)
        self.assertEqual(perms[1],perms[0][::-1]);self.assertEqual(perms[3],perms[2][::-1])
        self.assertNotEqual(perms[0],perms[2])


class LocalSheddingTopology(unittest.TestCase):
    def build(self,elements):
        # elements: (buses, terminal powers) with positive power entering the element.
        terminals=[b for b,_ in elements];powers=[p for _,ps in elements for p in ps]
        term_first=[];count=0
        for buses in terminals:term_first.append(count);count+=len(buses)
        incidence=defaultdict(list)
        for i,buses in enumerate(terminals):
            for t,bus in enumerate(buses):incidence[bus].append((i,t))
        return dict(incidence),powers,term_first,terminals

    def test_radial_subtree_only(self):
        # src -0-> a -1-> b (load L1); a -2-> c (load L2); src -3-> d (load L3)
        inc,p,first,terms=self.build([(["src","a"],[10,-10]),(["a","b"],[4,-4]),(["a","c"],[6,-6]),(["src","d"],[5,-5])])
        loads={"b":["L1"],"c":["L2"],"d":["L3"]}
        self.assertEqual(sorted(downstream_loads(["a"],inc,p,first,terms,loads,0)),["L1","L2"])
        self.assertEqual(downstream_loads(["b"],inc,p,first,terms,loads,1),["L1"])

    def test_flow_direction_is_respected_and_dead_branches_skipped(self):
        # a feeds b; e (dead, zero flow) joins b to x which hosts L9
        inc,p,first,terms=self.build([(["a","b"],[3,-3]),(["b","x"],[0,0])])
        self.assertEqual(downstream_loads(["b"],inc,p,first,terms,{"b":["L1"],"x":["L9"]},0),["L1"])


class Configurations(unittest.TestCase):
    def test_new_and_legacy_configurations(self):
        cfg=config(RUNTIME/"configs/core-v2.toml")
        self.assertEqual((cfg["power"]["operator"],cfg["recovery"]["crew_model"],cfg["task_a"]["shapley_permutations"]),("local","per_damage",4))
        for legacy in (config(),config(RUNTIME/"configs/core-unbounded-12.toml")):
            self.assertEqual((legacy["power"].get("operator","uniform_grid"),legacy["recovery"].get("repair_time_model","fixed"),
                              legacy["task_c"].get("od_tie_break","asset_id")),("uniform_grid","fixed","asset_id"))
        self.assertNotIn("extends",cfg)

    def test_invalid_settings_are_rejected(self):
        import tempfile
        # Mode selectors are checked everywhere; mode parameters only when their mode is on (model-v2).
        cases=[('[power]\noperator = "lp"\n',False),('[task_c]\nod_tie_break = "random"\n',False),
               ('[recovery]\ncrew_model = "all"\n',False),('[power]\nemergency_factor = 0.5\n',True),
               ('[recovery]\npower_units_per_crew = 0\n',True),('[recovery]\npower_repair_minutes_partial = [300.0, 100.0]\n',True),
               ('[task_a]\nshapley_permutations = 3\n',True),('[recovery]\ncrew_congestion_cap = 0.5\n',True)]
        for text,needs_v2 in cases:
            prefix='extends = ["%s"]\n'%(RUNTIME/"configs/model-v2.toml") if needs_v2 else ""
            with tempfile.NamedTemporaryFile("w",suffix=".toml",delete=False) as f:f.write(prefix+text)
            with self.subTest(text=text),self.assertRaises(ValueError):config(f.name)
            Path(f.name).unlink()

if __name__=="__main__":
    unittest.main()
