#!/usr/bin/env python3
"""Pure/fast smoke checks for rare native zero-chord Heed resume policy.
No ROOT, no Garfield++, no Geant4, no archived data modifications.
"""
from __future__ import annotations
import hashlib
import numpy as np
import pi_six_geometry_native_timing as physics
import pi_2m_zero_path_compatible_resume as guard


def fake(length="0",out=(0.1,0.,0.2),edep="0.0465400322447"):
    return {
      "event":"4040","hit_index":"3","track_id":"2",
      "pdg":"13","det_id":"110","station":"3",
      "path_length_mm":str(length),"edep_eV":str(edep),
      "u_in_mm":"0.1","v_in_mm":"0","w_in_mm":"0.2",
      "u_out_mm":str(out[0]),"v_out_mm":str(out[1]),
      "w_out_mm":str(out[2]),
    }


def main():
    row=fake()
    assert guard.coincident_zero_chord(row)
    before=len(guard.seen)
    pos,q=guard.guarded_birth_sample(row,lam=3.0,marks=np.array([1,2]),seed=123)
    assert pos.shape==(0,3) and q.shape==(0,)
    assert len(guard.seen)==before+1
    assert guard.seen[-1]["event"]==4040 and guard.seen[-1]["det_id"]==110
    assert not guard.coincident_zero_chord(fake(length="0.1"))
    assert not guard.coincident_zero_chord(fake(edep="15"))
    assert not guard.coincident_zero_chord(fake(out=(0.2,0.,0.2)))
    for wrong in (fake(edep="15"),fake(out=(0.2,0.,0.2))):
        try:
            guard.guarded_birth_sample(wrong,3.,np.array([1]),123)
        except ValueError:
            pass
        else:
            raise AssertionError("Bad zero-length geometry was not rejected")
    valid=fake(length="2.",out=(0.1,2.0,0.2),edep="42")
    marks=np.array([1,1,2,3],dtype=np.int64)
    a,b=physics.birth_sample(valid,3.,marks,12345)
    c,d=guard.guarded_birth_sample(valid,3.,marks,12345)
    assert np.array_equal(a,c) and np.array_equal(b,d)
    assert guard.ORIGINAL is physics.birth_sample
    assert physics.aggregate.__globals__["birth_sample"] is physics.birth_sample
    # The wrapper is a distinct file. Original SHA-salted frozen
    # physical fold source inputs MUST NOT be edited by this mitigation.
    assert guard.HERE.name=="pi_2m_zero_path_compatible_resume.py"
    print("PASS: strict Edep<1eV and geometrically identical zero chord "
          "-> zero Heed primaries; native gas exposure retained; "
          "nonzero Heed samples and random draws bitwise unchanged; "
          "ambiguous/bad zero chords rejected.")


if __name__=="__main__":
    main()
