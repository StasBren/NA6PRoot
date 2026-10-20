#!/usr/bin/env python3
"""Build electron-ONLY two-anode PI viewer with TRUE Geant4 local u-v tracks.

Consumes the already-generated 3+3-mm TWO-ANODE viewer (selected real native
events and Garfield electron arrival times) and JOINs TGeo-local gas entrance
and exit coordinates from the existing first-chunk Geant4 CSV exports.

No Garfield/Heed resampling, analog pulse synthesis, ROOT export or Geant4.
Rejects missing/mismatched track keys. Does not infer paths from entrance
angles and does not represent arrival-time distributions as ASIC signals.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

from mwpc_checkerboard_projection import build_chambers
from pi_native_centered_timing_pilot import BASE, SOURCES, EXPORT
from pi_centered_electron_sampler import GEOMETRY_ID

HERE=Path(__file__).resolve().parent
PAGES=BASE/"pi_fixed_wire_maps"
OLD=PAGES/"pi_two_anode_signal_explorer.html"
NEW=PAGES/"pi_two_anode_electron_explorer.html"
TEMPLATE=HERE/"pi_two_anode_electron_explorer_template.html"
START='<script id="twoWireData" type="application/json">'
END='</script>'
MARKER='__PI_TWO_ANODE_ELECTRON_DATA__'
REQUIRED={'event','hit_index','track_id','det_id','station','pdg',
          'u_in_mm','v_in_mm','u_out_mm','v_out_mm',
          'w_in_mm','w_out_mm','theta_u_in_deg','theta_w_in_deg'}
GEOM={c.det_id:c for group in build_chambers() for c in group}


def load_old(path):
    text=path.read_text(encoding='utf-8')
    if text.count(START)!=1:
        raise ValueError("Original 3+3mm two-anode HTML data not found")
    raw=text.split(START,1)[1].split(END,1)[0]
    data=json.loads(raw)
    if data['meta'].get('geometry_id')!=GEOMETRY_ID:
        raise ValueError("Original viewer did not use the correct centered field")
    if not data.get('examples'):
        raise ValueError("No two-anode native events available")
    return data,hashlib.sha256(raw.encode('utf-8')).hexdigest()


def simplify_event(original):
    wires=original['wires']
    if len(wires)!=2 or wires[0]==wires[1] or len(original['anodes'])!=2:
        raise ValueError("Requested event is not a two-distinct-anode event")
    anodes=[]
    for wire,part in zip(wires,original['anodes']):
        if int(part['wire'])!=int(wire):
            raise ValueError("Original wire/channel mismatch")
        arrivals=[float(t) for t in part['arrivalTimes_ns']]
        if any(t<0 for t in arrivals):
            raise ValueError("Negative collected-electron arrival time")
        if not arrivals:
            raise ValueError("Both physical anodes must have an arrival")
        if len(arrivals)!=int(part['nTimed']):
            raise ValueError("Original stored collected arrival times truncated")
        anodes.append({
            'wire':int(wire),'times_ns':arrivals,
            'first_ns':min(arrivals),
            'n_timed':len(arrivals),
            'n_collected':int(part['nPrimaryCollected']),
        })
    if sum(x['n_collected'] for x in anodes)!=int(original['nCollected']):
        raise ValueError("Collected primary electron count not conserved")
    st=int(original['st'])
    n=77 if st==0 else 128
    if not all(0<=int(w)<n for w in wires):
        raise ValueError("Physical wire number outside nominal finite grid")
    return {
        'src':original['src'],'st':st,'det':int(original['det']),
        'layer':original['layer'],'event':int(original['event']),
        'track':int(original['track']),'hitIndex':int(original['hitIndex']),
        'pdg':int(original['pdg']),
        'thetaU':float(original['thetaU']),'thetaW':float(original['thetaW']),
        'edep_eV':float(original['edep_eV']),
        'path_mm':float(original['path_mm']),
        'nPrimary':int(original['nPrimary']),
        'nCollected':int(original['nCollected']),
        'N_t100':int(original['N_t100']),
        'wires':[int(w) for w in wires],
        'anodes':anodes,
        'deltaFirst_ns':(anodes[1]['first_ns']-anodes[0]['first_ns']),
    }


def build(args):
    original,sha=load_old(args.original)
    events={}
    selected={src:set() for src in SOURCES}
    total=0
    for group,items in original['examples'].items():
        new=[]
        for item in items:
            e=simplify_event(item)
            source=e['src']
            if source not in selected or e['det'] not in GEOM or (
                GEOM[e['det']].station!=e['st']
            ):
                raise ValueError("Physical chamber/source label invalid")
            if group!=f"{source}|{e['det']}":
                raise ValueError("Original group/chamber inconsistent")
            identity=(e['event'],e['hitIndex'],e['track'])
            if identity in selected[source]:
                raise ValueError("Duplicate event/hit/track in curated two-wire selection")
            selected[source].add(identity)
            new.append(e)
            total+=1
        events[group]=new
    if not total:
        raise ValueError("No valid examples")

    # Lookup ONLY selected events, and attach native u,v,w entrance and
    # exit from the real Geant4 TGeo-local export (never angle extrapolation).
    found={src:set() for src in SOURCES}
    index={
        src:{(e['event'],e['hitIndex'],e['track']):e
             for key,ls in events.items() for e in ls if e['src']==src}
        for src in SOURCES
    }
    for src,channel in SOURCES.items():
        path=args.export_dir/f"{channel}_geant_gas_crossings.csv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing original native TGeo geometry export {path}")
        with path.open(newline='',encoding='utf-8') as stream:
            reader=csv.DictReader(stream)
            if not REQUIRED.issubset(reader.fieldnames or []):
                raise ValueError(f"Native Geant4 exporter lacks local u/v/w: {path}")
            for row in reader:
                identity=tuple(int(row[k]) for k in ('event','hit_index','track_id'))
                item=index[src].get(identity)
                if item is None:
                    continue
                if identity in found[src]:
                    raise ValueError("Repeated native hit identity "+str(identity))
                found[src].add(identity)
                if (abs(int(row['pdg']))!=13 or
                    int(row['det_id'])!=item['det'] or
                    int(row['station'])!=item['st'] or
                    abs(float(row['theta_u_in_deg'])-item['thetaU'])>0.0002 or
                    abs(float(row['theta_w_in_deg'])-item['thetaW'])>0.0002):
                    raise ValueError("Original two-anode record does not match native Geant4 hit")
                segment={}
                for direction in ('in','out'):
                    for coord in 'uvw':
                        val=float(row[f"{coord}_{direction}_mm"])
                        if not (-1.e6<val<1.e6):
                            raise ValueError("Nonfinite or enormous native local coordinate")
                        segment[f"{coord}_{direction}_mm"]=val
                    if abs(segment[f"v_{direction}_mm"])>3.05:
                        raise ValueError("Geant4 local hit endpoint outside symmetric gas gap")
                item['segment']=segment
        if found[src]!=selected[src]:
            raise ValueError(f"Cannot construct true Geant4 u-v tracks for {src}: "
                             f"missing {len(selected[src]-found[src])} events")
        print(f"PASS {src}: real TGeo-local muon entrance/exit joined for "
              f"{len(found[src]):,} selected 2-wire events",flush=True)

    geometry=[{'id':c.det_id,'st':c.station,'x':c.x,'y':c.y,
               'hx':c.half_x,'hy':c.half_y,
               'row':c.row,'col':c.col,'layer':c.layer}
              for c in GEOM.values()]
    compact={
      'sources':list(SOURCES),
      'layout':geometry,'examples':events,
      'meta':{
        'geometry_id':GEOMETRY_ID,
        'sample':'CURATED, FIRST 3000 GENERATED EVENTS PER SOURCE: true '
                 'native Geant4 gas-hit muon segments with exactly two '
                 'distinct physical anodes collecting primary electrons',
        'physics':'Real TGeo-local muon u-v paths; independently sampled '
                  'Heed primaries; VERIFIED symmetric Garfield B0 endpoint '
                  'and drift-time data only',
        'electronSignal':'Collected PRIMARY ELECTRON COUNTS per time bin, '
                         'NOT avalanche charge, voltage, ASIC readout or '
                         'a physical summed-electronics channel',
        'waveScope':'A+B means only arithmetic sum of the two electron '
                    'count histograms. No analog voltage anywhere.',
        'wirePitch_mm':4.,'gap_half_mm':3.,
        'preserved_original_payload_sha256':sha,
        'unique_selected_native_events':total,
      }
    }
    template=args.template.read_text(encoding='utf-8')
    if template.count(MARKER)!=1:
        raise ValueError("Electron viewer template marker must occur exactly once")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(template.replace(MARKER,json.dumps(
        compact,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')),
        encoding='utf-8')
    result={
      'geometry_id':GEOMETRY_ID,'curated_events':total,
      'counts_per_source':{src:len(found[src]) for src in SOURCES},
      'data_source':'Original 3+3mm two-anode electron arrival-time records '
                    'plus native Geant4 local u/v/w endpoint coordinates',
      'old_viewer_payload_sha256':sha,
      'new_electronic_voltage_model_used':False,
      'native_geant4_recomputed':False,
      'heed_or_garfield_recomputed':False,
      'cathode_strips_used':False,
    }
    report=args.out.with_name('pi_two_anode_electron_explorer_manifest.json')
    report.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',
                      encoding='utf-8')
    print(f"PASS: {total:,} real two-wire events with explicit TGeo local u-v track")
    print("Electron-only independent Garfield arrival-time histograms ready:",args.out)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=OLD)
    p.add_argument('--export-dir',type=Path,default=EXPORT)
    p.add_argument('--template',type=Path,default=TEMPLATE)
    p.add_argument('--out',type=Path,default=NEW)
    args=p.parse_args()
    build(args)


if __name__=='__main__':
    main()
