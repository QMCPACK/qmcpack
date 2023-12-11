from __future__ import print_function

# Cusp corrections for gaussian orbitals

# From "Scheme for adding electron-nucleus cusps to Gaussian orbitals" A. Ma, D. Towler, N. D. Drummond, R. J. Needs, Journal of Chemical Physics 122, 224322(2005) https://doi.org/10.1063/1.1940588

# Also qmc_algorithms/Wavefunctions/CuspCorrection.ipynb


from sympy import *
import gaussian_orbitals
import read_qmcpack
import math
import numpy as np



alpha = IndexedBase('alpha')
rc = Symbol('r_c')
X1,X2,X3,X4,X5 = symbols('X_1 X_2 X_3 X_4 X_5')
r = Symbol('r')
Zeff = Symbol('Z_eff')
sgn = Symbol('sgn')
C = Symbol('C')




def solve_for_alpha():
  p = alpha[0] + alpha[1]*r + alpha[2]*r**2 + alpha[3]*r**3 + alpha[4]*r**4


  # Constraint equations
  # Value matches at r_c
  eq1 = Eq(p.subs(r,rc), X1)
  #  1st derivative matches at r_c
  eq2 = Eq(diff(p,r).subs(r,rc), X2)
  #  2nd derivative matches at r_c
  eq3 = Eq((diff(p,r,2)+diff(p,r)**2).subs(r,rc),X3)
  # Cusp condition
  eq4 = Eq(diff(p,r).subs(r,0),X4)
  # Value of phi tilde at r=0
  eq5 = Eq(p.subs(r,0),X5)

  sln = solve([eq1, eq2, eq3, eq4, eq5],[alpha[0], alpha[1], alpha[2], alpha[3], alpha[4]])

  print('Symbolic solution for alpha:')
  for i,s in enumerate(sln[0]):
    print(alpha[i],' = ',s)
  print()

  return sln[0]

alpha_sln = solve_for_alpha()


def solve_for_alpha(Xvals, rc_val):
  svals = {rc: rc_val, X1:Xvals[1], X2:Xvals[2], X3:Xvals[3], X4:Xvals[4], X5:Xvals[5]}
  alpha_vals = []
  for i,s in enumerate(alpha_sln):
    alpha_vals.append(s.subs(svals))
  return alpha_vals

def simple_X_vals():
  rc_val = 0.1
  Xvals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] # first entry is unused

  print('For X: ',Xvals[1:])
  alpha_vals = solve_for_alpha(Xvals, rc_val)
  for i,a in enumerate(alpha_vals):
    print(alpha[i], ' = ', a)
  #svals = {rc: rc_val, X1:Xvals[1], X2:Xvals[2], X3:Xvals[3], X4:Xvals[4], X5:Xvals[5]}
  #for i,s in enumerate(alpha_sln):
  #  print(alpha[i], ' = ', s.subs(svals))

def evalX(gto, Z_val, rc_val):
  Xvals = [0.0]*6
  val_rc, grad_rc, lap_rc  = [v[0] for v in gto.eval_vgl(rc_val, 0.0, 0.0)]
  val_zero, grad_zero, lap_zero  = [v[0] for v in gto.eval_vgl(0.0, 0.0, 0.0)]
  Xvals[1] = log(val_rc)
  Xvals[2] = grad_rc[0]/val_rc
  Xvals[3] = (lap_rc - 2.0*grad_rc[0]/rc_val)/val_rc
  Xvals[4] = -Z_val
  Xvals[5] = log(abs(val_zero)) # initially use phi at 0
  return Xvals


def output_required_xvals_and_alphas(Xvals, alpha_vals):
  print('Xvals = ',Xvals)
  print('  // From gen_cusp_corr.py')
  for i in range(5):
    print('  CHECK(X[%d] == Approx(%.15f));'%(i,Xvals[i+1]))

  print()
  print('  // From gen_cusp_corr.py')
  for i in range(5):
    print('  CHECK(cusp.alpha[%d] == Approx(%.15f));'%(i,alpha_vals[i]))

def output_array(v, name):
  print('  // From gen_cusp_corr.py')
  for i,a in enumerate(v):
    print('  CHECK(%s[%d] == Approx(%.15f));'%(name, i, a))
  print()


def del_spherical(e, r):
    """Compute Laplacian for expression e with respect to symbol r.
        Currently works only with radial dependence"""
    t1 = r*r*diff(e, r)
    t2 = diff(t1, r)/(r*r)
    return simplify(t2)

def get_symbolic_effective_local_energy():
  p = alpha[0] + alpha[1]*r + alpha[2]*r**2 + alpha[3]*r**3 + alpha[4]*r**4
  R_sym = exp(p)
  phi_tilde = R_sym
  effEl_sym = -S.Half * del_spherical(phi_tilde, r)/phi_tilde - Zeff/r
  return effEl_sym


def eval_El(El_sym, r_val, Zeff_val, alpha_vals):
  slist = {alpha[0]:alpha_vals[0], alpha[1]:alpha_vals[1], alpha[2]: alpha_vals[2], alpha[3]:alpha_vals[3],
           alpha[4]:alpha_vals[4], Zeff:Zeff_val, r:r_val}
  val = El_sym.subs(slist).evalf()
  return val


def get_grid():
  pos = []
  for i in range(10):
    rval = .012*(i+1)
    pos.append(rval)
  return pos

def get_original_local_energy(pos, gto, El_sym, alpha_vals, rc_val, Zeff_val):
  vals = []
  for rval in pos:
    val, grad, lap  = gto.eval_vgl(rval, 0.0, 0.0)
    real_el = -.5*lap[0]/val[0] - Zeff_val/rval
    vals.append(real_el)

  return vals


def get_current_local_energy(pos, gto, El_sym, alpha_vals, rc_val, dE, Zeff_val):
  vals = []
  for rval in pos:
    el = eval_El(El_sym, rval, 2.0, alpha_vals)
    if rval < rc_val:
      vals.append(el + dE)
    else:
      val, grad, lap  = gto.eval_vgl(rval, 0.0, 0.0)
      real_el = -.5*lap[0]/val[0] - Zeff_val/rval
      vals.append(real_el + dE)

  return vals



def get_ideal_local_energy(pos, rc_val, beta0_val=0.0):
  ideal_EL = []
  beta0 = Symbol('beta_0')
  beta_vals = [beta0, 3.25819, -15.0126, 33.7308, -42.8705, 31.2276, -12.1316, 1.94692]
  El_terms = [beta_vals[n]*r**(n+1) for n in range(1,8)]
  Z = Symbol('Z')
  Z_val = 2.0
  EL_ideal_sym = Z*Z*(beta0 + sum(El_terms))

  slist = {beta0: 0.0, Z:Z_val, r: rc_val}

  idealEL_at_rc = EL_ideal_sym.subs(slist).evalf()
  #print('idealEL at rc = ',idealEL_at_rc)
  beta0_val = (-idealEL_at_rc)/Z_val/Z_val
  #print('beta0_val = ',beta0_val)

  for rval in pos:
    slist = {beta0: beta0_val, Z:Z_val, r: rval}
    v = EL_ideal_sym.subs(slist).evalf()
    ideal_EL.append(v)

  return ideal_EL


def get_phi_tilde_and_derivatives():
  p = alpha[0] + alpha[1]*r + alpha[2]*r**2 + alpha[3]*r**3 + alpha[4]*r**4
  pt = C + sgn*exp(p)

  pt_dr = diff(pt, r)
  pt_d2r = del_spherical(pt, r)

  # See the generated python code
  #pt_func_str = lambdastr([C,sgn,alpha,r],pt)
  #pt_dr_func_str = lambdastr([C,sgn,alpha,r],pt_dr)
  #pt_d2r_func_str = lambdastr([C,sgn,alpha,r],pt_d2r)
  #print('phi tilde',pt_func_str)
  #print('dr phi tilde',pt_dr_func_str)
  #print('d2r phi tilde',pt_d2r_func_str)

  pt_func = lambdify([C,sgn,alpha,r],pt)
  pt_dr_func = lambdify([C,sgn,alpha,r],pt_dr)
  pt_d2r_func = lambdify([C,sgn,alpha,r],pt_d2r)
  return pt_func, pt_dr_func, pt_d2r_func



def values_for_He():
  rc_val = 0.1
  Z_val = 2

  basis_set,MO = read_qmcpack.parse_qmc_wf('he_sto3g.wfj.xml',['He'])
  he_gto = gaussian_orbitals.GTO(basis_set['He'])

  Xvals = evalX(he_gto, Z_val, rc_val)

  alpha_vals = solve_for_alpha(Xvals, rc_val)
  output_required_xvals_and_alphas(Xvals, alpha_vals)

  El_sym = get_symbolic_effective_local_energy()
  print("El_sym = ",El_sym)

  Zeff_val = 2.0
  el_at_rc = -eval_El(El_sym, rc_val, Zeff_val, alpha_vals)
  dE = el_at_rc
  print('el at rc_val = ',el_at_rc)


  pos = get_grid()

  current_EL = get_current_local_energy(pos, he_gto, El_sym, alpha_vals, rc_val, dE, Zeff_val)
  original_EL = get_original_local_energy(pos, he_gto, El_sym, alpha_vals, rc_val, Zeff_val)
  #print('Current effective local energy')
  #for (p,v) in zip(pos,current_EL):
  #  print(p,v)
  print("  // Grid for local energy evaluations")
  output_array(pos, "cusp.pos");

  print("  // Original local energy")
  output_array(original_EL, "cusp.ELorig")

  print("  // Current local energy")
  output_array(current_EL, "cusp.ELcurr")

  print("  // Ideal local energy")
  ideal_EL = get_ideal_local_energy(pos,rc_val)
  output_array(ideal_EL, "cusp.ELideal")

  chi2 = 0.0
  for rval, ideal, curr in zip(pos, ideal_EL, current_EL):
    chi2 += (ideal - curr)**2
  print('  CHECK(chi2 == Approx(%.10f)); '%chi2)


def split_by_angular_momentum(pos_list, elements, basis_sets, MO_matrix):
  phi_list = []
  eta_list = []
  basis_by_index = gaussian_orbitals.get_center_and_ijk_by_index(pos_list, elements, basis_sets)
  C_no_S = MO_matrix.copy()
  for pos_idx in range(len(pos_list)):
    C_phi = MO_matrix.copy()
    C_eta = MO_matrix.copy()
    #print(' pos idx ',pos_idx)
    for ao_idx in range(MO_matrix.shape[1]):
      #print('     ao idx ',ao_idx)
      center_idx, basis_set, angular_info = basis_by_index[ao_idx]
      #print('       orbtype  ',basis_set.orbtype)


      if basis_set.orbtype == 0 and center_idx == pos_idx:
        #print('Setting eta to zero',ao_idx,C_phi[0:7, ao_idx])
        # s-type, part of phi, but not eta
        C_eta[:, ao_idx] = 0.0
      else:
        # not s-type, part of eta, but not phi
        C_phi[:, ao_idx] = 0.0

      # No s orbitals on any site - used during evaluation
      if basis_set.orbtype == 0:
        C_no_S[:, ao_idx] = 0.0

    phi_list.append(C_phi)
    eta_list.append(C_eta)

  return phi_list, eta_list, C_no_S



class CuspCorrection:
  def __init__(self, cusp_data, pos_list, gtos, phi_list):
    self.cusp_data = cusp_data
    # List of positions of centers
    self.pos_list = pos_list
    self.gtos = gtos
    # MO matrix for the phi wavefunction for each center
    self.phi_list = phi_list
    # generate these functions for faster evaluation
    self.phi_tilde_func, self.phi_tilde_dr_func, self.phi_tilde_d2r_func = get_phi_tilde_and_derivatives()

  def eval_phiBar_vgl(self, pos, mo_idx, pos_center, gto_mo, cusp_orb_data):
    rel_pos = [pos[0]-pos_center[0], pos[1]-pos_center[1], pos[2]-pos_center[2]]
    r = math.sqrt(rel_pos[0]**2 + rel_pos[1]**2 + rel_pos[2]**2)
    if r <= cusp_orb_data.Rc:
      a = cusp_orb_data.alpha
      #pr = a[0] + a[1]*r + a[2]*r*r + a[3]*r*r*r + a[4]*r*r*r*r
      #phiB = cusp_orb_data.C + cusp_orb_data.sg * exp(pr)
      phiB = self.phi_tilde_func(cusp_orb_data.C, cusp_orb_data.sg, a, r)
      grad = self.phi_tilde_dr_func(cusp_orb_data.C, cusp_orb_data.sg, a, r) * np.array(rel_pos)/r
      lap = self.phi_tilde_d2r_func(cusp_orb_data.C, cusp_orb_data.sg, a, r)
    else:
      phiB,grad,lap = gto_mo.eval_vgl_one_MO(pos[0], pos[1], pos[2], mo_idx)
      grad = [float(g) for g in grad]

    return phiB,grad,lap

  def compute_cusp_correction(self, center_idx, mo_idx, epos):
    C_phi = self.phi_list[center_idx]
    gto_mo = gaussian_orbitals.MolecularOrbital(self.gtos, C_phi)
    cusp_orb_data = self.cusp_data[center_idx]
    v,g,l = self.eval_phiBar_vgl(epos, mo_idx, self.pos_list[center_idx], gto_mo, cusp_orb_data[mo_idx])
    return v,g,l


def gen_correction_for_hcn():
  pos_list, elements = read_qmcpack.read_structure_file("hcn.structure.xml")
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("hcn.wfnoj.xml", elements)

  gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)

  #print("MO matrix")
  #for i in range(MO_matrix.shape[0]):
  #  print(i, MO_matrix[0:7, i])

  phi_list, eta_list, C_no_S = split_by_angular_momentum(pos_list, elements, basis_sets, MO_matrix)

  spo_name, cusp_data = read_qmcpack.read_cusp_correction_file("hcn_downdet.cuspInfo.xml")

  cusp = CuspCorrection(cusp_data, pos_list, gtos, phi_list)

  xgrid = get_grid()
  for center_idx in range(3):
    for mo_idx in range(7):
      #val_data = vals[(center_idx, mo_idx)][0]
      print('  //  Center ',center_idx, ' MO',mo_idx, 'rc = ',cusp_data[center_idx][mo_idx].Rc)
      for i,x in enumerate(xgrid):
        epos = np.array([x, 0.0, 0.0]) + pos_list[center_idx]
        v,g,l = cusp.compute_cusp_correction(center_idx, mo_idx, epos)
        print('  CHECK(rad_orb[%d] == Approx(%.10f)); // x = %g'%(i,v,x))
      print()

def gen_wavefunction_plus_correction_for_hcn():
  pos_list, elements = read_qmcpack.read_structure_file("hcn.structure.xml")
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("hcn.wfnoj.xml", elements)

  gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)

  phi_list, eta_list, C_no_S = split_by_angular_momentum(pos_list, elements, basis_sets, MO_matrix)

  spo_name, cusp_data = read_qmcpack.read_cusp_correction_file("hcn_downdet.cuspInfo.xml")

  cusp = CuspCorrection(cusp_data, pos_list, gtos, phi_list)

  # Bulk of the MO's - the non-cusp corrected parts
  other_mo = gaussian_orbitals.MolecularOrbital(gtos, C_no_S)

  # set electron position so it falls within rc for one center
  xyzgrid = [(-1.09, 0.0, 0.0)]
  #xyzgrid = [(0.0, 0.0, 0.0)]
  iat = 0
  for i,epos in enumerate(xyzgrid):
    mo_v, mo_g, mo_l = other_mo.eval_vgl(epos[0], epos[1], epos[2])
    for mo_idx in range(7):
      final_v = mo_v[mo_idx]
      final_g = mo_g[mo_idx, :]
      final_l = mo_l[mo_idx]
      for center_idx in range(3):
        cv,cg,cl = cusp.compute_cusp_correction(center_idx, mo_idx, epos)
        final_v += cv
        final_g += cg
        final_l += cl

      print('  # MO %d'%mo_idx)
      print('  CHECK(values[%d] == Approx(%.10f));'%(mo_idx,final_v))
      print('  CHECK(dpsi[%d][0] == Approx(%.10f));'%(mo_idx, final_g[0]))
      print('  CHECK(dpsi[%d][1] == Approx(%.10f));'%(mo_idx, final_g[1]))
      print('  CHECK(dpsi[%d][2] == Approx(%.10f));'%(mo_idx, final_g[2]))
      print('  CHECK(d2psi[%d] == Approx(%.10f));'%(mo_idx,final_l))
      print()

      print('  CHECK(all_values[%d][%d] == Approx(%.10f));'%(iat, mo_idx, final_v))
      print('  CHECK(all_grad[%d][%d][0] == Approx(%.10f));'%(iat, mo_idx, final_g[0]))
      print('  CHECK(all_grad[%d][%d][1] == Approx(%.10f));'%(iat, mo_idx, final_g[1]))
      print('  CHECK(all_grad[%d][%d][2] == Approx(%.10f));'%(iat, mo_idx, final_g[2]))
      print('  CHECK(all_lap[%d][%d] == Approx(%.10f));'%(iat, mo_idx, final_l))
      print()

def gen_wavefunction_plus_correction_for_ethanol():
  pos_list, elements = read_qmcpack.read_structure_file("ethanol.structure.xml")
  basis_sets, MO_matrix = read_qmcpack.parse_qmc_wf("ethanol.wfnoj.xml", elements)

  gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)

  phi_list, eta_list, C_no_S = split_by_angular_momentum(pos_list, elements, basis_sets, MO_matrix)

  spo_name, cusp_data = read_qmcpack.read_cusp_correction_file("ethanol_downdet.cuspInfo.xml")

  cusp = CuspCorrection(cusp_data, pos_list, gtos, phi_list)

  # Bulk of the MO's - the non-cusp corrected parts
  other_mo = gaussian_orbitals.MolecularOrbital(gtos, C_no_S)

  # set electron position so it falls within rc near O atom
  xyzgrid = [(-2.1, 0.5, 0.0)]
  #xyzgrid = [(0.0, 0.0, 0.0)]
  iat = 0
  for i,epos in enumerate(xyzgrid):
    mo_v, mo_g, mo_l = other_mo.eval_vgl(epos[0], epos[1], epos[2])
    for mo_idx in range(13):
      final_v = mo_v[mo_idx]
      final_g = mo_g[mo_idx, :]
      final_l = mo_l[mo_idx]
      for center_idx in range(9):
        cv,cg,cl = cusp.compute_cusp_correction(center_idx, mo_idx, epos)
        final_v += cv
        final_g += cg
        final_l += cl

      print('  # MO %d'%mo_idx)
      print('  CHECK(values[%d] == Approx(%.10f));'%(mo_idx,final_v))
      print('  CHECK(dpsi[%d][0] == Approx(%.10f));'%(mo_idx, final_g[0]))
      print('  CHECK(dpsi[%d][1] == Approx(%.10f));'%(mo_idx, final_g[1]))
      print('  CHECK(dpsi[%d][2] == Approx(%.10f));'%(mo_idx, final_g[2]))
      print('  CHECK(d2psi[%d] == Approx(%.10f));'%(mo_idx,final_l))
      print()

      print('  CHECK(all_values[%d][%d] == Approx(%.10f));'%(iat, mo_idx, final_v))
      print('  CHECK(all_grad[%d][%d][0] == Approx(%.10f));'%(iat, mo_idx, final_g[0]))
      print('  CHECK(all_grad[%d][%d][1] == Approx(%.10f));'%(iat, mo_idx, final_g[1]))
      print('  CHECK(all_grad[%d][%d][2] == Approx(%.10f));'%(iat, mo_idx, final_g[2]))
      print('  CHECK(all_lap[%d][%d] == Approx(%.10f));'%(iat, mo_idx, final_l))
      print()



if __name__ == '__main__':
  #simple_X_vals()
  #values_for_He()
  #gen_correction_for_hcn()
  #gen_wavefunction_plus_correction_for_hcn()
  gen_wavefunction_plus_correction_for_ethanol()
