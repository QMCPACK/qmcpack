#!/usr/bin/env python
import numpy as np

def example_structure_h4():
  # hydrogen at rs=1.31
  from structure import Structure
  natom = 4
  alat = 3.3521298178767225
  axes = alat*np.eye(3)
  elem = ['H']*natom
  pos = np.array([
    [0, 0, 0], [alat/2., 0, 0], [0, alat/2, 0], [0, 0, alat/2]
  ])
  s1 = Structure(axes=axes, elem=elem, pos=pos, units='B')
  return s1
#end def example_structure_h4

def test_kf_rpa():
  from physical_system import generate_physical_system
  s1 = example_structure_h4()
  ps = generate_physical_system(
    structure = s1,
    net_charge = 1,
    net_spin = 1,
    H = 1
  )
  kfs = ps.kf_rpa()
  assert np.isclose(kfs[0], 1.465, atol=1e-3)
  assert np.isclose(kfs[1], 1.465/2**(1./3), atol=1e-3)
#end def test_kf_rpa

if __name__ == '__main__':
  test_kf_rpa()
#end __main__
