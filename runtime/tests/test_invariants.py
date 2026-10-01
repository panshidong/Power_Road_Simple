"""Destination-host tests. These have NOT been run on the source host."""
from __future__ import annotations
import contextlib, copy, io, json, math, tempfile, unittest
from pathlib import Path
from unittest.mock import MagicMock, patch
import numpy as np
from austin_runtime.cli import main
from austin_runtime.common import config, read_json
from austin_runtime.experiments import objective
from austin_runtime.plan import B_RULES, C_RULES, jobs
from austin_runtime.power import PowerEngine
from austin_runtime.rankings import construct, order
from austin_runtime.scenarios import generate
from austin_runtime.simulation import Coupled, measure
from austin_runtime.traffic import TrafficState


class Invariants(unittest.TestCase):
    def test_cli_plan_starts_no_engine(self):
        stream=io.StringIO()
        with patch("austin_runtime.power.PowerEngine.__init__",side_effect=AssertionError("Native engine started")), contextlib.redirect_stdout(stream):
            main(["plan"])
        plan=json.loads(stream.getvalue())
        self.assertGreater(plan["total_jobs"],0)
        self.assertEqual(plan["stages"]["a-main"],900)

    def test_matrix_has_all_strategies_and_unique_ids(self):
        cfg=config()
        self.assertEqual({j["strategy"] for j in jobs(cfg,"b-main")},set(B_RULES))
        self.assertEqual({j["strategy"] for j in jobs(cfg,"c-main")},set(C_RULES))
        self.assertEqual(len(jobs(cfg,"b-sensitivity")),40)
        for stage in ("a-greedy","a-alpha","c-alternatives","b-sensitivity"):
            plan=jobs(cfg,stage)
            self.assertEqual(len(plan),len({j["id"] for j in plan}))

    def test_parallel_arcs_centroids_and_physical_closures(self):
        links=[dict(from_node=u,to_node=v,link_id=i+1) for i,(u,v) in enumerate([(3,4),(3,4),(3,1),(1,5),(4,5)])]
        state=TrafficState(links,[7,2,1,1,3],[0]*5,[1]*5,{},first_thru=3)
        self.assertEqual(state.travel(3,[5]),(5.,5)) # no shortcut through centroid 1
        self.assertEqual(state.travel(1,[5]),(1.,5)) # centroid is a valid origin
        closed=TrafficState(links,[7,9999,1,1,3],[0]*5,[1,0,1,1,1],{},first_thru=3)
        self.assertEqual(closed.travel(3,[5]),(10.,5))
        self.assertEqual(state.distances([5],reverse=True)[3],5.)

    def test_metrics_use_left_steps_and_report_censoring(self):
        def e(t,p,r,x,a):return dict(time=t,power_func=p,road_func=r,weighted_power_func=p,shelter_access=r,electric=x,access=a)
        events=[e(0,.5,.8,[0,0],[0,0]),e(10,1,.9,[1,0],[1,0]),e(20,1,1,[1,1],[1,1])]
        eq=dict(electric_weight=.5,access_weight=.5,cri_threshold=.9,access_threshold=.9)
        result=measure(events,[1,2],eq)
        self.assertAlmostEqual(result["triangle_area"],8.)
        self.assertAlmostEqual(result["p90_restore"],19.)
        self.assertAlmostEqual(result["gini_restore"],1/6)
        self.assertEqual(result["censored_cri_zones"],0)
        censored=measure(events[:-1],[1,2],eq)
        self.assertEqual(censored["censored_cri_zones"],1)

    def test_weighted_objective_and_guardrail(self):
        cfg=config();ref=dict(triangle_area=10.,gini_restore=.5);m=dict(triangle_area=8.,gini_restore=.25)
        self.assertAlmostEqual(objective(m,{"kind":"weighted","metric":"gini_restore","lambda":.5},ref,cfg),1.05)
        self.assertEqual(objective(m,{"kind":"guardrail","metric":"gini_restore"},ref,cfg),8.)
        m["gini_restore"]=.6
        self.assertGreater(objective(m,{"kind":"guardrail","metric":"gini_restore"},ref,cfg),100000.)

    def test_fault_types_do_not_resample_assets(self):
        catalog={"assets":{**{f"power:p{i}":dict(kind="power") for i in range(30)},**{f"road:{i}-{i+1}":dict(kind="road") for i in range(40)}}}
        cfg=config();partial=copy.deepcopy(cfg);partial["power"]["partial_failure_probability"]=1.
        full=copy.deepcopy(cfg);full["power"]["partial_failure_probability"]=0.
        a=generate(catalog,partial,"main",0);b=generate(catalog,full,"main",0)
        self.assertEqual(set(a["damage"]),set(b["damage"]))
        self.assertEqual({k:v for k,v in a["damage"].items() if k.startswith("road:")},{k:v for k,v in b["damage"].items() if k.startswith("road:")})
        self.assertTrue(all(v==0 for k,v in b["damage"].items() if k.startswith("power:")))

    def test_missing_scores_are_recorded_or_rejected(self):
        cfg=config();assets=["power:p1","road:3-4"]
        catalog={"centrality":dict(zip(assets,[.1,.2]))};tables={"JSH":{"power:p1":-.1}}
        sequence,missing=order(assets,"JSH",tables,catalog,cfg)
        self.assertEqual(missing,["road:3-4"]);self.assertEqual(set(sequence),set(assets))
        cfg["task_a"]["unseen_asset_policy"]="error"
        with self.assertRaises(ValueError):order(assets,"JSH",tables,catalog,cfg)

    def test_negative_marginals_and_checkpoint_resume(self):
        cfg=config();cfg["task_a"]["shapley_permutations"]=3
        scenario=dict(seed=10,damage={"power:p1":0.,"road:3-4":0.})
        engine=MagicMock();engine.catalog={"centrality":{"power:p1":.1,"road:3-4":.2}}
        engine.coalition.side_effect=lambda s,repaired:(.5+sum({"power:p1":-.1,"road:3-4":.3}[a] for a in repaired),.4+(.2 if "power:p1" in repaired else 0.))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"checkpoint.json"
            result=construct(engine,scenario,cfg,path,"sig")
            self.assertAlmostEqual(result["scores"]["JSH"]["power:p1"],-.1)
            self.assertAlmostEqual(sum(result["scores"]["JSH"].values()),.2)
            engine.coalition.reset_mock()
            self.assertEqual(result,construct(engine,scenario,cfg,path,"sig"))
            engine.coalition.assert_not_called()

    def test_simultaneous_repairs_and_unreachable_censoring(self):
        engine=Coupled.__new__(Coupled);engine.cfg=config()
        engine.cfg["recovery"].update(power_repair_minutes=10.,road_repair_minutes=10.)
        engine.catalog={"assets":{"power:s":{"kind":"power"},"road:3-4":{"kind":"road"}},"critical_substations":["s"]}
        engine.depot=1;engine.zones=[1];engine.excluded_zones=[]
        engine.healthy_power={"served_kw":1.,"nominal_kw":1.}
        seen=[]
        def state(remaining,penalty=None):
            seen.append(set(remaining))
            return dict(traffic=None,power_func=float("power:s" not in remaining),road_func=float("road:3-4" not in remaining))
        engine.state=state
        engine.site_travel=lambda traffic,origin,asset,roundtrip_lastmile=False:(0.,3)
        engine.record=lambda t,remaining,s,critical:dict(time=t,remaining=dict(remaining),power_func=s["power_func"],road_func=s["road_func"],weighted_power_func=s["power_func"],shelter_access=s["road_func"],electric=[s["power_func"]],access=[s["road_func"]])
        scenario=dict(id="example",seed=1,damage={"power:s":0.,"road:3-4":0.})
        run=engine.simulate(scenario,list(scenario["damage"]))
        self.assertTrue(run["complete"])
        self.assertEqual(seen,[set(scenario["damage"]),set()])
        self.assertEqual([e["time"] for e in run["events"]],[0.,10.])
        self.assertEqual(len(run["dispatch"]),2)
        engine.site_travel=lambda traffic,origin,asset,roundtrip_lastmile=False:(math.inf,origin)
        run=engine.simulate(scenario,list(scenario["damage"]))
        self.assertFalse(run["complete"])
        self.assertEqual(run["stop_reason"],"unreachable_pending_jobs")
        self.assertEqual(run["metrics"]["completion_or_censor_time"],1440.)

    def test_native_fault_isolation_does_not_edit_impedance(self):
        engine=PowerEngine.__new__(PowerEngine);engine.dss=MagicMock();d=engine.dss
        engine.scratch=Path("/tmp/austin-unit-scratch");engine.prepared=Path("/tmp/austin-unit-prepared")
        engine.loads={"R":{}};engine.catalog={"assets":{"power:p1":{"isolation_transformers":["transformer.t1"]},"power:p2":{"isolation_transformers":["transformer.t2"]}}}
        d.Loads.Count.return_value=0;d.RegControls.AllNames.return_value=["r1","r2"]
        active={}
        d.RegControls.Name.side_effect=lambda name:active.update(name=name)
        d.RegControls.Transformer.side_effect=lambda:{"r1":"Transformer.t1","r2":"t2"}[active["name"]]
        ratings,disabled=engine._compile("R",{"power:p1":0.,"power:p2":.4},.8,[])
        self.assertEqual(ratings,{"transformer.t2":.4});self.assertEqual(disabled,{"t1"})
        commands=[c.args[0] for c in d.call_args_list]
        self.assertIn("Disable transformer.t1",commands);self.assertIn("Disable RegControl.r1",commands)
        self.assertNotIn("Disable RegControl.r2",commands)
        self.assertFalse(any("edit transformer" in command.lower() for command in commands))


if __name__=="__main__":unittest.main()
