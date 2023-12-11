
# Compute values for Gaussian Type Orbitals
# Used in test_MO.cpp

import read_qmcpack
import gaussian_orbitals
import numpy as np

def gen_He():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("he_sto3g.wfj.xml",['He'])
  gto = gaussian_orbitals.GTO(basis_sets['He'])
  for pos in ([0.1, 0.0, 0.0], [1.0, 0.0, 0.0]):
    atomic_orbs = gto.eval_v(*pos)

    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,atomic_orbs[0]))
    print('')

    v,g,l = gto.eval_vgl(*pos)
    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,v[0]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(0,g[0][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(0,g[0][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(0,g[0][2]))
    print('  CHECK(d2psi[%d] == Approx(%15.10g));'%(0,l[0]))
    print('')
   
    v,g,h = gto.eval_vgh(*pos)
    gh    = gto.eval_gradhess(*pos)
    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,v[0]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(0,g[0][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(0,g[0][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(0,g[0][2]))
    print('  //Hessian (xx,xy,xz,yy,yz,zz) ')
    print('  CHECK(dhpsi[%d][0] == Approx(%15.10g));'%(0,h[0][0]))
    print('  CHECK(dhpsi[%d][1] == Approx(%15.10g));'%(0,h[0][1]))
    print('  CHECK(dhpsi[%d][2] == Approx(%15.10g));'%(0,h[0][2]))
    print('  CHECK(dhpsi[%d][3] == Approx(%15.10g));'%(0,h[0][3]))
    print('  CHECK(dhpsi[%d][4] == Approx(%15.10g));'%(0,h[0][4]))
    print('  CHECK(dhpsi[%d][5] == Approx(%15.10g));'%(0,h[0][5]))
    print('  //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz) ')
    print('  CHECK(dghpsi[%d][0] == Approx(%15.10g));'%(0,gh[0][0]))
    print('  CHECK(dghpsi[%d][1] == Approx(%15.10g));'%(0,gh[0][1]))
    print('  CHECK(dghpsi[%d][2] == Approx(%15.10g));'%(0,gh[0][2]))
    print('  CHECK(dghpsi[%d][3] == Approx(%15.10g));'%(0,gh[0][3]))
    print('  CHECK(dghpsi[%d][4] == Approx(%15.10g));'%(0,gh[0][4]))
    print('  CHECK(dghpsi[%d][5] == Approx(%15.10g));'%(0,gh[0][5]))
    print('  CHECK(dghpsi[%d][6] == Approx(%15.10g));'%(0,gh[0][6]))
    print('  CHECK(dghpsi[%d][7] == Approx(%15.10g));'%(0,gh[0][7]))
    print('  CHECK(dghpsi[%d][8] == Approx(%15.10g));'%(0,gh[0][8]))
    print('  CHECK(dghpsi[%d][9] == Approx(%15.10g));'%(0,gh[0][9]))
    print('')
  
       

def gen_Ne():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("ne_def2_svp.wfnoj.xml",['Ne'])
  gto = gaussian_orbitals.GTO(basis_sets['Ne'])
  for pos in ([0.1, 0.0, 0.0], [1.0, 0.0, 0.0]):
    atomic_orbs = gto.eval_v(*pos)
    mol_orbs =  np.dot(MO_matrix, atomic_orbs)

    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,mol_orbs[0]))

    v,g,l = gto.eval_vgl(*pos)
    mo_v = np.dot(MO_matrix, v)
    mo_g = np.dot(MO_matrix, g)
    mo_l = np.dot(MO_matrix, l)
    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,mo_v[0]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(0,mo_g[0][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(0,mo_g[0][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(0,mo_g[0][2]))
    print('  CHECK(d2psi[%d] == Approx(%15.10g));'%(0,mo_l[0]))

  
    v,g,h = gto.eval_vgh(*pos)
    gh    = gto.eval_gradhess(*pos)
    mo_v = np.dot(MO_matrix,v)
    mo_g = np.dot(MO_matrix,g)
    mo_h = np.dot(MO_matrix,h)
    mo_gh = np.dot(MO_matrix,gh)
    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(0,mo_v[0]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(0,mo_g[0][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(0,mo_g[0][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(0,mo_g[0][2]))
    print('  //Hessian (xx,xy,xz,yy,yz,zz) ')
    print('  CHECK(dhpsi[%d][0] == Approx(%15.10g));'%(0,mo_h[0][0]))
    print('  CHECK(dhpsi[%d][1] == Approx(%15.10g));'%(0,mo_h[0][1]))
    print('  CHECK(dhpsi[%d][2] == Approx(%15.10g));'%(0,mo_h[0][2]))
    print('  CHECK(dhpsi[%d][3] == Approx(%15.10g));'%(0,mo_h[0][3]))
    print('  CHECK(dhpsi[%d][4] == Approx(%15.10g));'%(0,mo_h[0][4]))
    print('  CHECK(dhpsi[%d][5] == Approx(%15.10g));'%(0,mo_h[0][5]))
    print('  //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz) ')
    print('  CHECK(dghpsi[%d][0] == Approx(%15.10g));'%(0,mo_gh[0][0]))
    print('  CHECK(dghpsi[%d][1] == Approx(%15.10g));'%(0,mo_gh[0][1]))
    print('  CHECK(dghpsi[%d][2] == Approx(%15.10g));'%(0,mo_gh[0][2]))
    print('  CHECK(dghpsi[%d][3] == Approx(%15.10g));'%(0,mo_gh[0][3]))
    print('  CHECK(dghpsi[%d][4] == Approx(%15.10g));'%(0,mo_gh[0][4]))
    print('  CHECK(dghpsi[%d][5] == Approx(%15.10g));'%(0,mo_gh[0][5]))
    print('  CHECK(dghpsi[%d][6] == Approx(%15.10g));'%(0,mo_gh[0][6]))
    print('  CHECK(dghpsi[%d][7] == Approx(%15.10g));'%(0,mo_gh[0][7]))
    print('  CHECK(dghpsi[%d][8] == Approx(%15.10g));'%(0,mo_gh[0][8]))
    print('  CHECK(dghpsi[%d][9] == Approx(%15.10g));'%(0,mo_gh[0][9]))
    print('')
   

def gen_HCN():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("hcn.wfnoj.xml",['N','C','H'])
  pos_list, elements = read_qmcpack.read_structure_file("hcn.structure.xml")

  gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)
  pos = [0.0, 0.0, 0.0]
  atomic_orbs = gtos.eval_v(*pos)
  #print 'first MO',MO_matrix[0,:]
  #print 'atomic_orbs',atomic_orbs
  mol_orbs =  np.dot(MO_matrix, atomic_orbs)
  print('  // Generated from gen_mo.py for position %s'%str(pos))
  for i in range(7):
    print('  CHECK(values[%d] == Approx(%15.10g));'%(i,mol_orbs[i]))

  v,g,l = gtos.eval_vgl(*pos)
  mo_v = np.dot(MO_matrix, v)
  mo_g = np.dot(MO_matrix, g)
  mo_l = np.dot(MO_matrix, l)
  print('  // Generated from gen_mo.py for position %s'%str(pos))
  for i in range(7):
    print('  CHECK(values[%d] == Approx(%15.10g));'%(i,mo_v[i]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(i,mo_g[i][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(i,mo_g[i][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(i,mo_g[i][2]))
    print('  CHECK(d2psi[%d] == Approx(%15.10g));'%(i,mo_l[i]))
    print('')

  v,g,h = gtos.eval_vgh(*pos)
  gh    = gtos.eval_gradhess(*pos)
  mo_v = np.dot(MO_matrix,v)
  mo_g = np.dot(MO_matrix,g)
  mo_h = np.dot(MO_matrix,h)
  mo_gh = np.dot(MO_matrix,gh)
  for i in range(7):
    print('  // Generated from gen_mo.py for position %s'%str(pos))
    print('  CHECK(values[%d] == Approx(%15.10g));'%(i,mo_v[i]))
    print('  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(i,mo_g[i][0]))
    print('  CHECK(dpsi[%d][1] == Approx(%15.10g));'%(i,mo_g[i][1]))
    print('  CHECK(dpsi[%d][2] == Approx(%15.10g));'%(i,mo_g[i][2]))
    print('  //Hessian (xx,xy,xz,yy,yz,zz) ')
    print('  CHECK(dhpsi[%d][0] == Approx(%15.10g));'%(i,mo_h[i][0]))
    print('  CHECK(dhpsi[%d][1] == Approx(%15.10g));'%(i,mo_h[i][1]))
    print('  CHECK(dhpsi[%d][2] == Approx(%15.10g));'%(i,mo_h[i][2]))
    print('  CHECK(dhpsi[%d][3] == Approx(%15.10g));'%(i,mo_h[i][3]))
    print('  CHECK(dhpsi[%d][4] == Approx(%15.10g));'%(i,mo_h[i][4]))
    print('  CHECK(dhpsi[%d][5] == Approx(%15.10g));'%(i,mo_h[i][5]))
    print('  //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz) ')
    print('  CHECK(dghpsi[%d][0] == Approx(%15.10g));'%(i,mo_gh[i][0]))
    print('  CHECK(dghpsi[%d][1] == Approx(%15.10g));'%(i,mo_gh[i][1]))
    print('  CHECK(dghpsi[%d][2] == Approx(%15.10g));'%(i,mo_gh[i][2]))
    print('  CHECK(dghpsi[%d][3] == Approx(%15.10g));'%(i,mo_gh[i][3]))
    print('  CHECK(dghpsi[%d][4] == Approx(%15.10g));'%(i,mo_gh[i][4]))
    print('  CHECK(dghpsi[%d][5] == Approx(%15.10g));'%(i,mo_gh[i][5]))
    print('  CHECK(dghpsi[%d][6] == Approx(%15.10g));'%(i,mo_gh[i][6]))
    print('  CHECK(dghpsi[%d][7] == Approx(%15.10g));'%(i,mo_gh[i][7]))
    print('  CHECK(dghpsi[%d][8] == Approx(%15.10g));'%(i,mo_gh[i][8]))
    print('  CHECK(dghpsi[%d][9] == Approx(%15.10g));'%(i,mo_gh[i][9]))
    print('')
  

def gen_HCN_force():
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("hcn.wfnoj.xml",['N','C','H'])
  ionpos, elements = read_qmcpack.read_structure_file("hcn.structure.xml")

  pos = [0.02, -0.1, 0.05]

  delta=1.0e-6
  deltainv=1.0/delta

  Natom=ionpos.shape[0];
  norb=7

  print('  // Generated from gen_mo.py for position %s'%str(pos))

  for iat in range(0,Natom):
    for idim in range(0,3):
      ionpos_p=np.array(ionpos)
      ionpos_m=np.array(ionpos)
      
      ionpos_p[iat][idim]+=delta
      ionpos_m[iat][idim]-=delta

      gtos_p = gaussian_orbitals.GTO_centers(ionpos_p, elements, basis_sets)
      gtos_m = gaussian_orbitals.GTO_centers(ionpos_m, elements, basis_sets)
      
      atomic_orbs_p = gtos_p.eval_v(*pos)
      atomic_orbs_m = gtos_m.eval_v(*pos)
      
      mol_orbs_p =  np.dot(MO_matrix, atomic_orbs_p)
      mol_orbs_m =  np.dot(MO_matrix, atomic_orbs_m)

      v_p,g_p,l_p = gtos_p.eval_vgl(*pos)
      mo_v_p = np.dot(MO_matrix, v_p)
      mo_g_p = np.dot(MO_matrix, g_p)
      mo_l_p = np.dot(MO_matrix, l_p)

      v_m,g_m,l_m = gtos_m.eval_vgl(*pos)
      mo_v_m = np.dot(MO_matrix, v_m)
      mo_g_m = np.dot(MO_matrix, g_m)
      mo_l_m = np.dot(MO_matrix, l_m)

      dmo_v = 0.5*deltainv*(mo_v_p-mo_v_m)
      dmo_g = 0.5*deltainv*(mo_g_p-mo_g_m)
      dmo_l = 0.5*deltainv*(mo_l_p-mo_l_m)
      print("//============== Ion ",iat," Component ",idim,"===================")
      for iorb in range(0,norb):
        print('  CHECK( dionpsi[0][%d][%d]       == Approx(%15.10g) );  '%(iorb,idim,dmo_v[iorb]))
        print('  CHECK( diongradpsi[0][%d](%d,0) == Approx(%15.10g) );  '%(iorb,idim,dmo_g[iorb][0]))
        print('  CHECK( diongradpsi[0][%d](%d,1) == Approx(%15.10g) );  '%(iorb,idim,dmo_g[iorb][1]))
        print('  CHECK( diongradpsi[0][%d](%d,2) == Approx(%15.10g) );  '%(iorb,idim,dmo_g[iorb][2]))
        print('  CHECK( dionlaplpsi[0][%d][%d]  == Approx(%15.10g) );  '%(iorb,idim,dmo_l[iorb]))
    
  print('  // Generated from gen_mo.py for position %s'%str(pos))
 # for i in range(7):
 #   print '  CHECK(values[%d] == Approx(%15.10g));'%(i,mol_orbs[i])

#  v,g,l = gtos.eval_vgl(*pos)
#  mo_v = np.dot(MO_matrix, v)
#  mo_g = np.dot(MO_matrix, g)
#  mo_l = np.dot(MO_matrix, l)
#  print '  // Generated from gen_mo.py for position %s'%str(pos)
#  for i in range(7):
#    print '  CHECK(values[%d] == Approx(%15.10g));'%(i,mo_v[i])
#    print '  CHECK(dpsi[%d][0] == Approx(%15.10g));'%(i,mo_g[i][0])
#    print '  CHECK(d2psi[%d] == Approx(%15.10g));'%(i,mo_l[i])
#    print ''

#  v,g,h = gtos.eval_vgh(*pos)
#  gh    = gtos.eval_gradhess(*pos)
#  mo_v = np.dot(MO_matrix,v)
#  mo_g = np.dot(MO_matrix,g)
#  mo_h = np.dot(MO_matrix,h)
#  mo_gh = np.dot(MO_matrix,gh)
 
if __name__ == '__main__':
  #gen_He()
  #gen_Ne()
  #gen_HCN()
  gen_HCN_force()
