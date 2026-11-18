import unittest
import numpy as np
from wire_core import (nearest_wire,n_anodes,build_lut_kernel,compare_one,
                       observables,stable_hash,cell_ids,physical_wire_charges)

class TestNA60MultiAnode(unittest.TestCase):
    def setUp(self):
        ii=np.repeat(np.arange(30),30)
        vv=np.tile(np.arange(30),30)
        self.central=build_lut_kernel(ii,vv,np.ones(900,dtype=int),np.zeros(900,dtype=int))
        # A synthetic balanced {-1,0,1} endpoint kernel, all cells present.
        self.triple=build_lut_kernel(np.repeat(ii,3),np.repeat(vv,3),
                         np.ones(2700,dtype=int),np.tile([-1,0,1],900))

    def test_native_station_lattice(self):
        self.assertEqual(n_anodes(0),102)
        self.assertEqual(n_anodes(3),171)
        self.assertEqual(n_anodes(5),171)
        u=np.array([-1.5,1.5,4.5,7.5])
        k,local=nearest_wire(u,0)
        self.assertTrue(np.allclose(local,0.))
        self.assertEqual(len(set(k)),4)
        k3,loc3=nearest_wire(np.array([-3.,0,3.]),3)
        self.assertTrue(np.allclose(loc3,0))

    def test_three_plus_anodes_from_single_trajectory(self):
        # Six well-separated source cells; each electron has the SAME local
        # LUT offset 0, but the complete muon touches SIX physical anodes.
        u=np.array([-7.5,-4.5,-1.5,1.5,4.5,7.5])
        v=np.zeros_like(u);gain=np.array([20,30,40,50,60,70.])
        g,l,ge,le=compare_one(u,v,gain,0,self.central,rng=np.random.default_rng(42))
        self.assertEqual(len(g),6)
        self.assertEqual(len(l),6)
        self.assertEqual(observables(l)['n_gt5pct'],6)
        self.assertAlmostEqual(sum(l.values()),sum(gain))
        self.assertEqual((ge,le),(0,0))

    def test_random_adjacent_offsets_but_unbounded_total_multiplicity(self):
        rng=np.random.default_rng(57)
        n=2500
        u=np.linspace(-30,30,n)
        v=np.zeros(n)
        g,l,_,_=compare_one(u,v,np.ones(n),3,self.triple,rng=rng)
        self.assertGreaterEqual(len(g),20)
        self.assertGreaterEqual(len(l),20)
        self.assertGreaterEqual(observables(l)['n_gt1pct'],18)

    def test_same_zero_gain_and_edge_accounting(self):
        u=np.array([0.,3.,6.]);v=np.zeros(3);weight=np.array([0,9,1.])
        g,l,ge,le=compare_one(u,v,weight,3,self.central,rng=np.random.default_rng(1))
        self.assertEqual(g,l)
        self.assertEqual(observables(g)['n_positive'],2)
        self.assertEqual(observables({})['n_gt5pct'],0)
        self.assertEqual((ge,le),(0,0))

    def test_lut_count_and_missing_cells(self):
        with self.assertRaises(ValueError):
            build_lut_kernel(np.array([0]),np.array([0]),np.array([1]),np.array([0]))
        with self.assertRaises(ValueError):
            build_lut_kernel(np.repeat(np.arange(30),30),np.tile(np.arange(30),30),
                             np.ones(900),np.full(900,2))

    def test_stable_seed(self):
        self.assertEqual(stable_hash('Jpsi',1,3),stable_hash('Jpsi',1,3))
        self.assertNotEqual(stable_hash('Jpsi',1,3),stable_hash('Jpsi',1,4))

if __name__=='__main__':unittest.main()