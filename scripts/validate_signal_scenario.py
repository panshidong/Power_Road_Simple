#!/usr/bin/env python3
"""Exercise provisional signal-power-road joins with a static capacity-loss scenario."""
import csv, hashlib, json, math
from collections import Counter
from pathlib import Path
from build_dataset import ROOT, OUT, dump, write_csv
from validate_solvers import tapb


def rows(path):
    with path.open() as f:return list(csv.DictReader(f))


def main():
    matches=rows(OUT/'coupling/signal_road_power_provisional.csv')
    counts=Counter(r['substation_id'] for r in matches)
    substation=sorted(counts,key=lambda s:(-counts[s],s))[0]
    selected=[r for r in matches if r['substation_id']==substation]
    signal_ids={r['signal_id'] for r in selected}
    approaches=[r for r in rows(OUT/'coupling/signal_approach_links.csv') if r['signal_id'] in signal_ids]
    link_ids={int(r['link_id']) for r in approaches}
    assert selected and link_ids and len(link_ids)==len(approaches)
    folder=ROOT/'data/scenarios/signal_outage_demo';folder.mkdir(parents=True,exist_ok=True)
    write_csv(folder/'signals.csv',list(selected[0]),selected)
    write_csv(folder/'affected_links.csv',list(approaches[0]),approaches)
    source=(OUT/'road/Austin_net.tntp').read_text().splitlines()
    lines=source[:8]
    assert len(source[8:])==18710
    for link_id,line in enumerate(source[8:],1):
        parts=line.split(';')[0].split()
        if link_id in link_ids:parts[2]=format(float(parts[2])*.5,'.12g')
        lines.append('\t'.join(parts)+'\t;')
    network=folder/'Austin_signal_outage_net.tntp'
    network.write_text('\n'.join(lines)+'\n')
    for i,(a,b) in enumerate(zip(source[8:],lines[8:]),1):
        x=list(map(float,a.split(';')[0].split()));y=list(map(float,b.split(';')[0].split()))
        assert x[:2]+x[3:]==y[:2]+y[3:]
        assert y[2]==x[2]*(.5 if i in link_ids else 1)
    baseline=json.loads((ROOT/'reports/tapb_validation.json').read_text())
    assert baseline['passed'] and baseline['network_sha256']==hashlib.sha256((OUT/'road/Austin_net.tntp').read_bytes()).hexdigest()
    result=tapb(ROOT/'vendor/tap-b/bin/tap',network,'tapb_signal_outage','tapb_signal_outage_validation.json')
    assert result['trips_sha256']==baseline['trips_sha256']
    result.update(scope='Hypothetical signals-off capacity assignment; the substation choice uses unverified proximity candidates. No power-flow fault simulation.',
        provisional_substation_id=substation,signals_affected=len(selected),directed_approach_links_affected=len(link_ids),
        capacity_factor=.5,baseline_tstt_source_units=baseline['tstt_source_units'],
        tstt_change_percent=100*(result['tstt_source_units']/baseline['tstt_source_units']-1),
        only_selected_capacity_fields_changed=True,real_event_or_calibrated_outage=False)
    dump(ROOT/'reports/tapb_signal_outage_validation.json',result)
    dump(folder/'scenario.json',dict(description=result['scope'],substation_id=substation,
        selection_rule='Substation with most joint provisional road-power matches; lexicographic tie break',
        signals_affected=len(selected),links_affected=len(link_ids),capacity_factor=.5,
        signal_power_mapping='synthetic nearest load within 200m',
        approach_mapping='all incoming physical links at uniquely selected candidate node',
        restoration_or_power_contingency_solved=False))
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
