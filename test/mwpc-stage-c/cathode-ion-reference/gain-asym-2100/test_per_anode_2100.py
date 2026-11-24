#!/usr/bin/env python3
"""Synthetic unit tests of exact per-anode aggregation (not detector data)."""
from __future__ import annotations
import math
import tempfile
import unittest
from pathlib import Path
import pandas as pd

from per_anode_census_2100 import read_per_anode


class PerAnodeCensusTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix="test_per_anode_")
        self.addCleanup(self.tmp.cleanup)
        self.home=Path(self.tmp.name)
        self.run=self.home/"results"
        self.run.mkdir()
        self.event=self.home/"source"/"muon_0000"
        self.event.mkdir(parents=True)
        (self.run/"rootfiles.txt").write_text(str(self.event/"ion_births.root")+"\n")
        self.mu=pd.DataFrame([dict(muon_slot=0,channel="jpsi",station=0,det_id=1,
            theta_u_deg=3.0,n_primaries=2,
            ion_births_exact=30,collected_electrons=20)])
        self.pr=pd.DataFrame([
            dict(muon_slot=0,primary=0,t_primary_ns=10.0,
                 ion_births_exact=10,collected_electrons=5),
            dict(muon_slot=0,primary=1,t_primary_ns=20.0,
                 ion_births_exact=20,collected_electrons=15)])
        self.seed=pd.DataFrame([
           dict(primary_index=0,wire_index_global=7,ion_births=10,collected_e=5,
                ion_up=6,ion_down=4,ion_boundary=0,A_ion_hemi=0.2,
                mean_t_birth_ns=12.0,sigma_t_birth_ns=.5,
                mean_t_collected_ns=13.0,sigma_t_collected_ns=.3,
                sigma_w_birth_mm=.02),
           dict(primary_index=1,wire_index_global=7,ion_births=20,collected_e=15,
                ion_up=5,ion_down=15,ion_boundary=0,A_ion_hemi=-.5,
                mean_t_birth_ns=26.0,sigma_t_birth_ns=1.,
                mean_t_collected_ns=27.0,sigma_t_collected_ns=.4,
                sigma_w_birth_mm=.04)])
        self.save()

    def save(self):
        self.seed.to_csv(self.event/"microscopic_wire_seed.csv",index=False)

    def test_exact_hemisphere_and_parent_relative_times(self):
        frag,anode=read_per_anode(self.run,self.pr,self.mu)
        self.assertEqual(len(frag),2)
        self.assertEqual(len(anode),1)
        self.assertAlmostEqual(frag.iloc[0]["ion_birth_delay_ns"],2.)
        self.assertAlmostEqual(frag.iloc[1]["ion_birth_delay_ns"],6.)
        self.assertAlmostEqual(anode.iloc[0]["A_UD_exact"],(11-19)/30)
        self.assertAlmostEqual(anode.iloc[0]["ion_birth_delay_ns"],(10*2+20*6)/30)
        second=(10*(.5*.5+2*2)+20*(1+6*6))/30
        avg=(10*2+20*6)/30
        self.assertAlmostEqual(anode.iloc[0]["sigma_t_birth_relative_ns"],
                               math.sqrt(second-avg*avg))
        self.assertAlmostEqual(anode.iloc[0]["electron_collection_delay_ns"],
                               (5*3+15*7)/20)

    def test_hemispheric_boundary_excluded_from_asymmetry(self):
        # The first primary has four off-boundary births and six on the axis.
        self.seed.loc[0,"ion_up"]=1
        self.seed.loc[0,"ion_down"]=3
        self.seed.loc[0,"ion_boundary"]=6
        self.seed.loc[0,"A_ion_hemi"]=-.5
        self.save()
        frag,anode=read_per_anode(self.run,self.pr,self.mu)
        self.assertAlmostEqual(frag.iloc[0]["A_UD_exact"],-.5)
        self.assertAlmostEqual(anode.iloc[0]["A_UD_exact"],(1+5-(3+15))/(4+20))
        self.assertEqual(int(anode.iloc[0]["ion_births"]),30)

    def test_original_census_mismatch_must_stop(self):
        self.seed.loc[0,"ion_births"]=11
        self.save()
        with self.assertRaisesRegex(ValueError,"exact closure|birth"):
            read_per_anode(self.run,self.pr,self.mu)

    def test_original_counts_are_not_gains(self):
        frag,anode=read_per_anode(self.run,self.pr,self.mu)
        self.assertEqual(int(frag.iloc[0]["G_ion_primary"]),10)
        self.assertEqual(int(frag.iloc[1]["G_ion_primary"]),20)
        self.assertAlmostEqual(anode.iloc[0]["G_ion_anode_contribution"],30/2)


if __name__=="__main__":
    unittest.main()
