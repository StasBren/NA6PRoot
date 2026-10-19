#!/usr/bin/env python3
"""Fast synthetic PI anode drift-time histogram / probability regression.

No ROOT, Geant4, Heed or Garfield input file required.
"""
import numpy as np
from pi_native_drift_times import add_hit, blank, EDGES, TIMED_WINDOWS


def main():
    assert TIMED_WINDOWS == (50,100,200,400,800)
    assert EDGES[0]==0 and EDGES[-1]==np.inf
    b=blank()
    # Hit 1: two physical anodes, one collected e at 40ns,
    # and the other with a 150ns arrival.
    add_hit(b,[7,8],[40.,150.],0)
    # Hit 2: one anode, two electrons and one unknown time.
    add_hit(b,[7,7],[70.,np.nan],0)
    # Hit 3: no collection, which MUST remain in gas-hit denominator.
    add_hit(b,[],[],2)
    # Hit 4: a collected electron beyond the last finite histogram bin.
    add_hit(b,[9],[2500.],0)

    assert b["n"]==4
    assert b["collected"]==5
    assert b["validTimes"]==4
    assert b["unknownTimes"]==1
    assert b["edgeLoss"]==2
    assert b["zeroCollectedHits"]==1
    assert b["hist"].sum()==4
    assert b["hist"][-1]==1  # a real 2,500 ns tail event
    assert b["multiplicity"][0].tolist()==[3,1,0,0,0,0]
    assert b["multiplicity"][1].tolist()==[2,2,0,0,0,0]
    assert b["multiplicity"][2].tolist()==[2,1,1,0,0,0]
    assert b["multiplicity"][-1].tolist()==[1,2,1,0,0,0]
    assert np.all(b["multiplicity"].sum(axis=1)==b["n"])

    # Correctly preserve the individual anode definition: two electrons
    # on wire 7 still make ONE active anode.
    b2=blank()
    add_hit(b2,[7,7,7],[15.,20.,45.],0)
    assert b2["multiplicity"][-1].tolist()==[0,1,0,0,0,0]

    for bad in [
        ([1],[np.nan,np.nan],0),
        ([1],[8.],-1),
    ]:
        try:
            add_hit(blank(),*bad)
        except ValueError:
            pass
        else:
            raise AssertionError("Malformed timing physics input accepted")

    print("PASS: primary-electron Garfield time histogram with >2000ns overflow")
    print("PASS: zero-electron gas hits retained; unknown times never treated as t=0")
    print("PASS: timing-gated physical-wire multiplicity vs all-collected reference")
    print("PASS: distinct 77/128-wire identities, not electron number, determine N")


if __name__=="__main__":
    main()
