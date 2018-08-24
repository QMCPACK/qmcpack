

# Evaluate GTO's starting from a symbolic representation
# see qmc_algorithms/Wavefunctions/GaussianOrbitals.ipynb

from sympy import *
from collections import namedtuple, defaultdict
import numpy as np


CG_basis = namedtuple('CG_basis',['orbtype','nbasis','zeta','contraction_coeff'])

class GTO:
  def __init__(self, basis=None, center=[0.0, 0.0, 0.0]):
    x,y,z = symbols('x y z')
    alpha = Symbol('alpha', positive=True, real=True)
    r = Symbol('r',real=True,nonnegative=True)
    i,j,k = symbols('i j k',integer=True)
    N = Symbol('N')
    self.N = N
    norm1 = factorial(i)*factorial(j)*factorial(k)
    norm2 = factorial(2*i)*factorial(2*j)*factorial(2*k)
    norm_sym = (2*alpha/pi)**(3/S(4)) * sqrt((8*alpha)**(i+j+k)*norm1/norm2)

    gto_sym_raw = N * x**i * y**j * z**k * exp(-alpha *r**2)

    gto_sym = gto_sym_raw.subs(N, norm_sym)

    self.alpha = alpha
    self.x = x; self.y = y; self.z = z
    self.r = r
    self.i = i; self.j = j; self.k = k
    self.gto_sym = gto_sym
    self.compute_grad()

    self.ijk = get_ijk_by_type()

    self.basis = basis
    self.center = center


  def compute_grad(self):
    r2 = self.x**2 + self.y**2 + self.z**2
    gto_xyz = self.gto_sym.subs(self.r**2, r2)
    #print gto_xyz
    self.grad = [0]*3
    self.grad[0] = diff(gto_xyz, self.x).subs(r2, self.r**2)
    self.grad[1] = diff(gto_xyz, self.y).subs(r2, self.r**2)
    self.grad[2] = diff(gto_xyz, self.z).subs(r2, self.r**2)
    lap = diff(gto_xyz, self.x, 2) + \
          diff(gto_xyz, self.y, 2) + \
          diff(gto_xyz, self.z, 2)
    # Need to expand to avoid NaN's
    self.lap = expand(lap.subs(r2, self.r**2))

  def set_basis(self, basis):
    self.basis = basis

  def make_subs_list(self, i, j, k, x, y, z, alpha):
      subs_list = {self.i:i, self.j:j, self.k:k}
      r2 = x**2 + y**2 + z**2
      r_val = sqrt(r2)
      subs_list[self.r] = r_val
      subs_list[self.r**2] = r2
      subs_list[self.alpha] = alpha
      subs_list[self.x] = x
      subs_list[self.y] = y
      subs_list[self.z] = z
      return subs_list

  def eval_single_v(self, i, j, k, x, y, z, alpha):
      xc = x - self.center[0]
      yc = y - self.center[1]
      zc = z - self.center[2]
      sl1 = self.make_subs_list(i,j,k,xc,yc,zc,alpha)
      v = self.gto_sym.subs(sl1).evalf()
      return v

  def eval_single_vgl(self, i, j, k, x, y, z, alpha):
      xc = x - self.center[0]
      yc = y - self.center[1]
      zc = z - self.center[2]
      sl1 = self.make_subs_list(i,j,k,xc,yc,zc,alpha)
      v = self.gto_sym.subs(sl1).evalf()
      g = [grad.subs(sl1).evalf() for grad in self.grad]
      lap = self.lap.subs(sl1).evalf()
      return v,g,lap

  def eval_contraction_v(self, x, y, z, basis):
      vals = []
      # Expand each basis type by the angular momentum state
      angular_list = self.ijk[basis.orbtype]
      for i,j,k,name in angular_list:
        val = 0.0
        for idx in range(basis.nbasis):
          val += basis.contraction_coeff[idx] * self.eval_single_v(i,j,k,x,y,z,basis.zeta[idx])
        vals.append(val)

      return vals

  def eval_contraction_vgl(self, x, y, z, basis):
      vals = []
      grads = []
      laps = []

      angular_list = self.ijk[basis.orbtype]
      for i,j,k,name in angular_list:
        val = 0.0
        grad = [0.0, 0.0, 0.0]
        lap = 0.0
        for idx in range(basis.nbasis):
          c = basis.contraction_coeff[idx]
          v,g,l = self.eval_single_vgl(i,j,k,x,y,z,basis.zeta[idx])
          val += c*v
          lap += c*l
          grad = [c*g[m] + grad[m] for m in range(3)]
        vals.append(val)
        grads.append(grad)
        laps.append(lap)
      return vals, grads, laps

  def eval_v(self, x, y, z):
    vs = []
    for basis in self.basis:
        v = self.eval_contraction_v(x, y, z, basis)
        vs.extend(v)
    return vs

  def eval_vgl(self, x, y, z):
    vs = []
    grads = []
    lapls = []
    for basis in self.basis:
        v,g,l = self.eval_contraction_vgl(x, y, z, basis)
        vs.extend(v)
        grads.extend(g)
        lapls.extend(l)
    return vs, grads, lapls

# generated from qmcpack src/Numerics/codegen/read_order.py
# Only part of the function included for now
def get_ijk():
  ijk = []
  # S
  ijk.append( (0,0,0,"S") )
  # P
  ijk.append( (1,0,0,"X") )
  ijk.append( (0,1,0,"Y") )
  ijk.append( (0,0,1,"Z") )
  # D
  ijk.append( (2,0,0,"XX") )
  ijk.append( (0,2,0,"YY") )
  ijk.append( (0,0,2,"ZZ") )
  ijk.append( (1,1,0,"XY") )
  ijk.append( (1,0,1,"XZ") )
  ijk.append( (0,1,1,"YZ") )

  return ijk

def get_ijk_by_type():
  ijk_list = get_ijk()

  by_type = defaultdict(list)
  for i,j,k,name in ijk_list:
    L = i + j + k
    by_type[L].append((i,j,k,name))

  return by_type

def get_ijk_inverse_index(basis_set):
  ijk_list = get_ijk_by_type()
  by_index = list()
  for basis in basis_set:
    angular_list = ijk_list[basis.orbtype]
    for angular_info in angular_list:
      by_index.append( (basis, angular_info) )
  return by_index



# Collection of atoms with different types and positions
class GTO_centers:
  def __init__(self, pos_list, elements, basis_sets):
    self.gto_list = []
    for pos, element in zip(pos_list, elements):
      gto = GTO(basis_sets[element], pos)
      self.gto_list.append(gto)

  def eval_v(self, x, y, z):
    vs = []
    for gto in self.gto_list:
      v = gto.eval_v(x, y, z)
      vs.extend(v)
    return vs

  def eval_vgl(self, x, y, z):
    vs = []
    grads = []
    laps = []
    for gto in self.gto_list:
      v,g,l = gto.eval_vgl(x,y,z)
      vs.extend(v)
      grads.extend(g)
      laps.extend(l)

    return vs, grads, laps

def get_center_and_ijk_by_index(pos_list, elements, basis_sets):
  index = []
  for pos_idx, (pos, element) in enumerate(zip(pos_list, elements)):
    index_for_one = get_ijk_inverse_index(basis_sets[element])
    for basis,angular_info in index_for_one:
      index.append((pos_idx, basis, angular_info))
  return index

class MolecularOrbital:
  def __init__(self, gto, MO_matrix):
    self.gto = gto
    self.MO_matrix = MO_matrix

  def eval_v(self, x, y, z):
    ao_vals = self.gto.eval_v(x, y, z)
    mo_vals = np.dot(self.MO_matrix, ao_vals)
    return mo_vals

  def eval_v_one_MO(self, x, y, z, mo_idx):
    ao_vals = self.gto.eval_v(x, y, z)
    mo_val = np.dot(self.MO_matrix[mo_idx, :], ao_vals)
    return mo_val

  def eval_vgl(self, x, y, z):
    ao_v, ao_g, ao_l = self.gto.eval_vgl(x, y, z)
    mo_v = np.dot(self.MO_matrix, ao_v)
    mo_g = np.dot(self.MO_matrix, ao_g)
    mo_l = np.dot(self.MO_matrix, ao_l)

    return mo_v, mo_g, mo_l

  def eval_vgl_one_MO(self, x, y, z, mo_idx):
    ao_vals, ao_g, ao_l = self.gto.eval_vgl(x, y, z)
    mo_val = np.dot(self.MO_matrix[mo_idx, :], ao_vals)
    mo_g = np.dot(self.MO_matrix[mo_idx, :], ao_g)
    mo_l = np.dot(self.MO_matrix[mo_idx, :], ao_l)
    return mo_val, mo_g, mo_l



if __name__ == '__main__':
  gto = GTO()
