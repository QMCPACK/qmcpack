
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

def gen_Ne():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("ne_def2_svp.wfnoj.xml")
  gto = gaussian_orbitals.GTO(basis_sets['Ne'])
  for pos in ([0.00001, 0.0, 0.0], [1.0, 0.0, 0.0]):
    atomic_orbs = gto.eval_v(*pos)
    mol_orbs =  np.dot(MO_matrix, atomic_orbs)

    print '  // Generated from gen_mo.py for position %s'%str(pos)
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(0,mol_orbs[0])

    v,g,l = gto.eval_vgl(*pos)
    mo_v = np.dot(MO_matrix, v)
    mo_g = np.dot(MO_matrix, g)
    mo_l = np.dot(MO_matrix, l)
    print '  // Generated from gen_mo.py for position %s'%str(pos)
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(0,mo_v[0])
    print '  REQUIRE(dpsi[%d][0] == Approx(%15.10g));'%(0,mo_g[0][0])
    print '  REQUIRE(dpsi[%d][1] == Approx(%15.10g));'%(0,mo_g[0][1])
    print '  REQUIRE(dpsi[%d][2] == Approx(%15.10g));'%(0,mo_g[0][2])
    print '  REQUIRE(d2psi[%d] == Approx(%15.10g));'%(0,mo_l[0])


def gen_HCN():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("hcn.wfnoj.xml")
  pos_list, elements = read_qmcpack.read_structure_file("hcn.structure.xml")

  gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)
  pos = [0.0, 0.0, 0.0]
  atomic_orbs = gtos.eval_v(*pos)
  #print 'first MO',MO_matrix[0,:]
  #print 'atomic_orbs',atomic_orbs
  mol_orbs =  np.dot(MO_matrix, atomic_orbs)
  print '  // Generated from gen_mo.py for position %s'%str(pos)
  for i in range(7):
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(i,mol_orbs[i])

  v,g,l = gtos.eval_vgl(*pos)
  mo_v = np.dot(MO_matrix, v)
  mo_g = np.dot(MO_matrix, g)
  mo_l = np.dot(MO_matrix, l)
  print '  // Generated from gen_mo.py for position %s'%str(pos)
  for i in range(7):
    print '  REQUIRE(values[%d] == Approx(%15.10g));'%(i,mo_v[i])
    print '  REQUIRE(dpsi[%d][0] == Approx(%15.10g));'%(i,mo_g[i][0])
    print '  REQUIRE(dpsi[%d][1] == Approx(%15.10g));'%(i,mo_g[i][1])
    print '  REQUIRE(dpsi[%d][2] == Approx(%15.10g));'%(i,mo_g[i][2])
    print '  REQUIRE(d2psi[%d] == Approx(%15.10g));'%(i,mo_l[i])
    print ''



if __name__ == '__main__':
  #gen_He()
  #gen_Ne()
  gen_HCN()
