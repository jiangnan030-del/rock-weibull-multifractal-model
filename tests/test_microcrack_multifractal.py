import unittest
import numpy as np
from microcrack_multifractal import box_count_statistics,compute_multifractal_spectrum,generate_microcrack_network,select_scaling_window
class Tests(unittest.TestCase):
 def test_network_deterministic_and_inside(self):
  a=generate_microcrack_network(5,1,30,seed=7); b=generate_microcrack_network(5,1,30,seed=7)
  np.testing.assert_allclose(a.points,b.points); self.assertTrue(np.all(a.points>=0)); self.assertTrue(np.all(a.points<=5)); self.assertEqual(len(a.centers),10)
 def test_probabilities(self):
  n=generate_microcrack_network(5,1,40,seed=8)
  for s in box_count_statistics(n.points,5,[3,4,6,8]): self.assertAlmostEqual(float(s.probabilities.sum()),1,12)
 def test_uniform_cube_dimension(self):
  points=np.random.default_rng(10).uniform(0,1,(150000,3)); scales=box_count_statistics(points,1,[2,3,4,5,6,8,10]); spectrum=compute_multifractal_spectrum(scales,np.array([0.,1.]),.95,4)
  self.assertTrue(spectrum[0].threshold_passed); self.assertGreater(spectrum[0].f,2.8); self.assertLessEqual(spectrum[0].f,3.01); self.assertAlmostEqual(spectrum[1].alpha,spectrum[1].f,10)
 def test_window(self):
  x=np.arange(6.); start,end,a,f,passed=select_scaling_window(x,2*x+1,1.5*x-2,.99,3)
  self.assertTrue(passed); self.assertEqual((start,end),(0,6)); self.assertAlmostEqual(a.slope,2); self.assertAlmostEqual(f.slope,1.5)
if __name__=='__main__': unittest.main()
