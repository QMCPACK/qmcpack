
from collections import namedtuple, defaultdict
from sympy import *

# See the GaussianOrbitals notebook in the qmc_algorithms repo for more explanation,
#  especially about the normalization.

def single_gaussian():
  r = Symbol('r',nonnegative=True)
  alpha = Symbol('alpha')

  alpha_val = 3.0
  r_val = 1.2
  print 'alpha = ',alpha_val
  print 'r = ',r_val

  phi = exp(-alpha*r**2)
  print phi
  print 'f == ',phi.subs(alpha, alpha_val).subs(r,r_val)

  d_phi = diff(phi, r)
  print 'symbolic df = ',d_phi
  print 'df == ',d_phi.subs(alpha, alpha_val).subs(r,r_val)

  dd_phi = diff(phi, r, 2)
  print 'symbolic d2f = ',dd_phi
  print 'd2f == ',dd_phi.subs(alpha, alpha_val).subs(r,r_val)

  d3d_phi = diff(phi, r, 3)
  print 'symbolic d3f = ',d3d_phi
  print 'd3f == ',d3d_phi.subs(alpha, alpha_val).subs(r,r_val)


CG_basis = namedtuple('CG_basis',['orbtype','nbasis','zeta','contraction_coeff'])

# Read basis sets in Gamess format
def read_gms_basis(fname):
    element = ''
    basis_sets = dict()
    with open(fname,'r') as f:
        lines = f.readlines()
        idx = 0
        while idx < len(lines):
            line = lines[idx].strip()
            if line.startswith('!') or len(line) == 0:
                idx += 1
                continue
            if line.startswith('$DATA'):
                idx += 1
                element = lines[idx].strip()
                basis_sets[element] = defaultdict(list)
                while idx < len(lines):
                  idx += 1
                  if lines[idx].startswith('$END'):
                    break
                  orbtype, s_nval = lines[idx].strip().split()
                  nval = int(s_nval)
                  coef_list = []
                  zeta_list = []
                  for i in range(nval):
                      idx += 1
                      s_orbidx,s_zeta,s_c = lines[idx].strip().split()
                      orbidx = int(s_orbidx)
                      zeta = float(s_zeta)
                      c = float(s_c)
                      zeta_list.append(zeta)
                      coef_list.append(c)
                  assert(len(coef_list) == nval)
                  assert(len(zeta_list) == nval)
                  cg = CG_basis(orbtype, nval, zeta_list, coef_list)
                  basis_sets[element][orbtype].append(cg)
            idx +=1
    return basis_sets

def create_full_gto_norm_symbolic():
  # full normalization
  n1 = sympify('(2*alpha/pi)**(3/4)')
  n2 = sympify('(8*alpha)**(i+j+k)')

  i = Symbol('i')
  j = Symbol('j')
  k = Symbol('k')
  n3 = factorial(i)*factorial(j)*factorial(k)
  d1 = factorial(2*i)*factorial(2*j)*factorial(2*k)
  norm = n1*sqrt(n2*n3/d1)
  return norm

def create_radial_gto_norm_symbolic():
  # just the radial part of the normalization
  n1 = sympify('2*alpha**(3/4)*2**(3/4)* (1/pi)**(1/4)')
  #n2 = sympify('(8*alpha)**(i+j+k)')
  n2 = sympify('(8*alpha)**L')
  #i = Symbol('i')
  #j = Symbol('j')
  #k = Symbol('k')
  L = Symbol('L')
  #n3 = factorial(i)*factorial(j)*factorial(k)
  #d1 = factorial(2*i)*factorial(2*j)*factorial(2*k)
  n3 = factorial(L)
  d1 = factorial(2*L)
  #d2 = 2*(i+j+k) + 1

  d2 = 2*L + 1
  norm = n1*sqrt(n2*n3/d1/d2)
  return norm

def create_angular_gto_norm_symbolic():
  return create_full_gto_norm_symbolic() / create_radial_gto_norm_symbolic()

def create_gto_symbolic():
    pre = sympify('x**i * y**j * z**k * exp(-alpha *r**2)')
    return pre


def create_sym(ijk=[0,0,0]):
  gto_sym = create_gto_symbolic()
  gto = gto_sym.subs({Symbol('i'):ijk[0], Symbol('j'):ijk[1],Symbol('k'):ijk[2]})
  norm = create_radial_gto_norm_symbolic()
  l_val = sum(ijk)
  norm_s = norm.subs({Symbol('i'):ijk[0], Symbol('j'):ijk[1],Symbol('k'):ijk[2], Symbol('L'):l_val})
  print 'norm_s',norm_s
  norm = lambdify(Symbol('alpha'), norm_s)

  i = Symbol('i',integer=True)
  c = IndexedBase('c')
  alpha2 = IndexedBase('alpha')
  norm2 = IndexedBase('N')
  N_basis = Symbol('N_b', integer=True)
  cg_sym = Sum(norm2[i]*c[i]*gto.subs(Symbol('alpha'),alpha2[i]),(i,1,N_basis))
  return cg_sym,norm,norm2,c,alpha2,N_basis

def eval_sym(cg_sym, norm, norm2, N_basis, c, alpha2, h_basis):
  cg_unroll = cg_sym.subs(N_basis, h_basis.nbasis).doit()
  for i in range(h_basis.nbasis):
    cc = h_basis.contraction_coeff[i]
    cz = h_basis.zeta[i]
    cg_unroll = cg_unroll.subs(c[i+1],cc).subs(alpha2[i+1],cz).subs(norm2[i+1],norm(cz))
    print cc,cz,norm(cz),'normL',norm(1.0)
  return cg_unroll

def compute_from_sym(h_basis, ijk=[0,0,0]):

  cg_sym, norm, norm2, c, alpha2, N_basis = create_sym(ijk)
  print 'norm',norm
  cg = eval_sym(cg_sym, norm, norm2, N_basis, c, alpha2, h_basis)
  r = Symbol('r')
  x = Symbol('x')
  # setting x to 1.0 is important to compute just the radial part
  slist = {r:1.3, x:1.0}

  print cg
  print 'f = ',cg.subs(slist)

  d_cg = diff(cg, r);
  print 'df = ',d_cg.subs(slist)

  dd_cg = diff(cg, r, 2);
  print 'ddf = ',dd_cg.subs(slist)

  d3_cg = diff(cg, r, 3);
  print 'd3f = ',d3_cg.subs(slist)

# get i,j,k values.  Mirrors getABS in CartesianTensor.h

def convert_powers(s):
  i = s.count('x')
  j = s.count('y')
  k = s.count('z')
  return (i,j,k)

def get_ijk():
  ijk = []
  # S
  ijk.append('')
  # P
  ijk.append('x')
  ijk.append('y')
  ijk.append('z')
  # D
  ijk.append('xx')
  ijk.append('yy')
  ijk.append('zz')
  ijk.append('xy')
  ijk.append('xz')
  ijk.append('yz')
  # F
  ijk.append('xxx')
  ijk.append('yyy')
  ijk.append('zzz')
  ijk.append('xxy')
  ijk.append('xxz')
  ijk.append('yyx')
  ijk.append('yyz')
  ijk.append('zzx')
  ijk.append('zzy')
  ijk.append('xyz')
  # G
  ijk.append('xxxx')
  ijk.append('yyyy')
  ijk.append('zzzz')

  ijk.append('xxxy')
  ijk.append('xxxz')

  ijk.append('yyyx')
  ijk.append('yyyz')

  ijk.append('zzzx')
  ijk.append('zzzy')

  ijk.append('xxyy')
  ijk.append('xxzz')
  ijk.append('yyzz')
  ijk.append('xxyz')
  ijk.append('yyxz')
  ijk.append('zzxy')

  return [convert_powers(s) for s in ijk]



def compute_radial_values(p):
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)
  #print gto_s

  norm_s = create_angular_gto_norm_symbolic()
  #print 'norm',norm_s

  print_value = True
  print_grad = True
  print_lap = True
  print_hess = False
  print_ggg = False  # third derivative tensor

  ijk_s = get_ijk()
  for idx, (i,j,k) in enumerate(ijk_s):
    print ''
    expr = gto_s * norm_s
    l_val = i + j + k
    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k, Symbol('L'):l_val}
    rlist = {Symbol('x'):p[0], Symbol('y'):p[1], Symbol('z'):p[2]}
    val =  expr.subs(slist).subs(rlist).evalf()
    if print_value:
      print "  REQUIRE(ct.getYlm(%d) == Approx(%.12g));"%(idx, val)
    dx = diff(expr, Symbol('x'))
    dy = diff(expr, Symbol('y'))
    dz = diff(expr, Symbol('z'))
    dx_val =  dx.subs(slist).subs(rlist).evalf()
    dy_val =  dy.subs(slist).subs(rlist).evalf()
    dz_val =  dz.subs(slist).subs(rlist).evalf()
    if print_grad:
      print "  REQUIRE(ct.getGradYlm(%d)[0] == Approx(%.12g));"%(idx, dx_val)
      print "  REQUIRE(ct.getGradYlm(%d)[1] == Approx(%.12g));"%(idx, dy_val)
      print "  REQUIRE(ct.getGradYlm(%d)[2] == Approx(%.12g));"%(idx, dz_val)

    lap = diff(expr, Symbol('x'), 2) + diff(expr, Symbol('y'), 2) + diff(expr, Symbol('z'), 2)
    lap_val = lap.subs(slist).subs(rlist).evalf()
    if print_lap:
      print "  REQUIRE(ct.getLaplYlm(%d) == Approx(%.12g));"%(idx, lap_val)


    if print_hess:
      print ''
      axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
      for ii,si in enumerate(axis_syms):
        for jj,sj in enumerate(axis_syms):
          h_s = diff(diff(expr, si), sj)
          hess_val = h_s.subs(slist).subs(rlist).evalf()

          #print ii,jj,hess_val
          #if hess_val != 0:
          #  print "  REQUIRE(ct.getHessYlm(%d)(%d,%d) == Approx(%.12g));"%(idx, ii, jj, hess_val)
          print "  REQUIRE(ct.getHessYlm(%d)(%d,%d) == Approx(%.12g));"%(idx, ii, jj, hess_val)

      print ''

    if print_ggg:
      print ''
      axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
      for ii,si in enumerate(axis_syms):
        for jj,sj in enumerate(axis_syms):
          for kk,sk in enumerate(axis_syms):
            ggg_s = diff(diff(diff(expr, si), sj), sk)
            ggg_val = ggg_s.subs(slist).subs(rlist).evalf()

            print "  REQUIRE(ct.getGGGYlm(%d)[%d](%d,%d) == Approx(%.12g));"%(idx, ii, jj, kk, ggg_val)

      print ''



if __name__ == '__main__':
    # Generate data for unit tests below

    # For test_gaussian_basis

    # Data for 'Gaussian Combo'
    #basis_sets = read_gms_basis('sto3g.txt')
    #compute_from_sym(basis_sets['HYDROGEN']['S'][0])

    # Data for 'Gaussian Combo P'
    #basis_sets_cc = read_gms_basis('cc-pVDZ.txt')
    #compute_from_sym(basis_sets_cc['CARBON']['P'][0],[1,0,0])

    # Data for 'Gaussian Combo D'
    #basis_sets_cc = read_gms_basis('cc-pVDZ.txt')
    #compute_from_sym(basis_sets_cc['CARBON']['D'][0],[2,0,0])


    # Data for test_cartesian_tensor
    compute_radial_values([1.3,1.2,-0.5])

