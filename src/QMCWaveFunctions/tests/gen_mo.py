
# Compute values for Gaussian Type Orbitals
# Used in test_MO.cpp

import read_qmcpack
import gaussian_orbitals
import numpy as np

def gen_He():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("he_sto3g.wfj.xml")
  gto = gaussian_orbitals.GTO(basis_sets['He'])
  for pos in ([0.1, 0.0, 0.0], [1.0, 0.0, 0.0]):
    atomic_orbs = gto.eval_v(*pos)

    print '  // Generated from gen_mo.py for position %s'%str(pos)
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(0,atomic_orbs[0])
    print ''

    v,g,l = gto.eval_vgl(*pos)
    print '  // Generated from gen_mo.py for position %s'%str(pos)
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(0,v[0])
    print '  REQUIRE(dpsi[%d][0] == Approx(%15.10g));'%(0,g[0][0])
    print '  REQUIRE(dpsi[%d][1] == Approx(%15.10g));'%(0,g[0][1])
    print '  REQUIRE(dpsi[%d][2] == Approx(%15.10g));'%(0,g[0][2])
    print '  REQUIRE(d2psi[%d] == Approx(%15.10g));'%(0,l[0])
    print ''

if __name__ == '__main__':
  gen_He()
