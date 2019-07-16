#!/usr/bin/env python
import numpy as np
from forlib.gencell import nsupercell, generate_all_tilematrices
from forlib.gencell import generate_all_supercells

def test_nsupercell():
  # test up to nprim=63
  answers = np.array([1 , 7 , 13 , 35 , 31 , 91 , 57 , 155 , 130 , 217 , 133 , 455 , 183 , 399 , 403 , 651 , 307 , 910 , 381 , 1085 , 741 , 931 , 553 , 2015 , 806 , 1281 , 1210 , 1995 , 871 , 2821 , 993 , 2667 , 1729 , 2149 , 1767 , 4550 , 1407 , 2667 , 2379 , 4805 , 1723 , 5187 , 1893 , 4655 , 4030 , 3871 , 2257 , 8463 , 2850 , 5642 , 3991 , 6405 , 2863 , 8470 , 4123 , 8835 , 4953 , 6097 , 3541 , 14105 , 3783 , 6951 , 7410], dtype=int)
  for i, ans in enumerate(answers):
    nprim=i+1
    assert nsupercell(nprim) == ans

def test_tilematrix():
  # test a few random ones
  nprim=3
  nhnf = nsupercell(nprim)
  hnfs = generate_all_tilematrices(nprim, nhnf)
  ihnf = 1
  ans = np.array([[1, 0, 0], [0, 1, 1], [0, 0, 3]], dtype=int)
  assert np.allclose(hnfs[:, :, ihnf], ans)
  # N96
  nprim=4
  nhnf = nsupercell(nprim)
  hnfs = generate_all_tilematrices(nprim, nhnf)
  ihnf = 20
  ans = np.array([[1, 1, 0], [0, 2, 0], [0, 0, 2]])
  assert np.allclose( hnfs[:, :, ihnf], ans)

def test_opt_tile():
  from structure import optimal_tilematrix_hnf
  nprim=4  # 24 to 96
  axes=np.array([
    [2.6812469112636386, -1.5259544931042469e-4, -4.6349129093427388],
    [-3.039595757365180e-2, 9.4905815242660267, 0.0],
    [2.6812469112636386, -1.5259544931042469e-4, 4.6349129093427388]
  ])
  topt, ropt = optimal_tilematrix_hnf(axes, nprim)
  assert np.isclose(ropt, 5.354577733)
  assert np.allclose(topt, np.array([[1, 1, 0,], [2, 0, 0], [0, 0, 2]]))
  nprim=32  # 24 to 768
  topt, ropt = optimal_tilematrix_hnf(axes, nprim)
  assert np.isclose(ropt, 10.725, atol=1e-4)
  assert np.allclose(topt, np.array([[3, 1, -1], [-1, 1, 3], [4, 0, 4]]))

if __name__ == '__main__':
  #test_nsupercell()
  #test_tilematrix()
  test_opt_tile()
# end __main__
