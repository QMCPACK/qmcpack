
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

# generated from read_order.py
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
  # F
  ijk.append( (3,0,0,"XXX") )
  ijk.append( (0,3,0,"YYY") )
  ijk.append( (0,0,3,"ZZZ") )
  ijk.append( (2,1,0,"XXY") )
  ijk.append( (2,0,1,"XXZ") )
  ijk.append( (1,2,0,"YYX") )
  ijk.append( (0,2,1,"YYZ") )
  ijk.append( (1,0,2,"ZZX") )
  ijk.append( (0,1,2,"ZZY") )
  ijk.append( (1,1,1,"XYZ") )
  # G
  ijk.append( (4,0,0,"XXXX") )
  ijk.append( (0,4,0,"YYYY") )
  ijk.append( (0,0,4,"ZZZZ") )
  ijk.append( (3,1,0,"XXXY") )
  ijk.append( (3,0,1,"XXXZ") )
  ijk.append( (1,3,0,"YYYX") )
  ijk.append( (0,3,1,"YYYZ") )
  ijk.append( (1,0,3,"ZZZX") )
  ijk.append( (0,1,3,"ZZZY") )
  ijk.append( (2,2,0,"XXYY") )
  ijk.append( (2,0,2,"XXZZ") )
  ijk.append( (0,2,2,"YYZZ") )
  ijk.append( (2,1,1,"XXYZ") )
  ijk.append( (1,2,1,"YYXZ") )
  ijk.append( (1,1,2,"ZZXY") )
  # H
  ijk.append( (5,0,0,"XXXXX") )
  ijk.append( (0,5,0,"YYYYY") )
  ijk.append( (0,0,5,"ZZZZZ") )
  ijk.append( (4,1,0,"XXXXY") )
  ijk.append( (4,0,1,"XXXXZ") )
  ijk.append( (1,4,0,"YYYYX") )
  ijk.append( (0,4,1,"YYYYZ") )
  ijk.append( (1,0,4,"ZZZZX") )
  ijk.append( (0,1,4,"ZZZZY") )
  ijk.append( (3,2,0,"XXXYY") )
  ijk.append( (3,0,2,"XXXZZ") )
  ijk.append( (2,3,0,"YYYXX") )
  ijk.append( (0,3,2,"YYYZZ") )
  ijk.append( (2,0,3,"ZZZXX") )
  ijk.append( (0,2,3,"ZZZYY") )
  ijk.append( (3,1,1,"XXXYZ") )
  ijk.append( (1,3,1,"YYYXZ") )
  ijk.append( (1,1,3,"ZZZXY") )
  ijk.append( (2,2,1,"XXYYZ") )
  ijk.append( (2,1,2,"XXZZY") )
  ijk.append( (1,2,2,"YYZZX") )
  # I
  ijk.append( (6,0,0,"X6") )
  ijk.append( (0,6,0,"Y6") )
  ijk.append( (0,0,6,"Z6") )
  ijk.append( (5,1,0,"X5Y") )
  ijk.append( (5,0,1,"X5Z") )
  ijk.append( (1,5,0,"Y5X") )
  ijk.append( (0,5,1,"Y5Z") )
  ijk.append( (1,0,5,"Z5X") )
  ijk.append( (0,1,5,"Z5Y") )
  ijk.append( (4,2,0,"X4Y2") )
  ijk.append( (4,0,2,"X4Z2") )
  ijk.append( (2,4,0,"Y4X2") )
  ijk.append( (0,4,2,"Y4Z2") )
  ijk.append( (2,0,4,"Z4X2") )
  ijk.append( (0,2,4,"Z4Y2") )
  ijk.append( (4,1,1,"X4YZ") )
  ijk.append( (1,4,1,"Y4XZ") )
  ijk.append( (1,1,4,"Z4XY") )
  ijk.append( (3,3,0,"X3Y3") )
  ijk.append( (3,0,3,"X3Z3") )
  ijk.append( (0,3,3,"Y3Z3") )
  ijk.append( (3,2,1,"X3Y2Z") )
  ijk.append( (3,1,2,"X3Z2Y") )
  ijk.append( (2,3,1,"Y3X2Z") )
  ijk.append( (1,3,2,"Y3Z2X") )
  ijk.append( (2,1,3,"Z3X2Y") )
  ijk.append( (1,2,3,"Z3Y2X") )
  ijk.append( (2,2,2,"X2Y2Z2") )

  return ijk


# To create test cases for each routine
# evaluate                 print_value = True, others False
# evaluateAll              print_value,print_grad,print_lap = True, others False
# evaluateWithHessian      print_value,print_grad,print_hess = True, others False
# evaluateWithThirdDeriv   print_value,print_grad,print_hess,print_ggg = True
# evaluateThirdDerivOnly   print_ggg = True, others False

def compute_radial_values(p, print_value=False, print_grad=False, print_lap=False,
                          print_hess=False, print_ggg=False, using_soa=False):

  out = ''
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)
  #print gto_s

  norm_s = create_angular_gto_norm_symbolic()
  #print 'norm',norm_s


  ijk_s = get_ijk()
  for idx, (i,j,k,s) in enumerate(ijk_s):
    expr = gto_s * norm_s
    l_val = i + j + k
    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k, Symbol('L'):l_val}
    rlist = {Symbol('x'):p[0], Symbol('y'):p[1], Symbol('z'):p[2]}
    val =  expr.subs(slist).subs(rlist).evalf()
    if print_value:
      if using_soa:
        out += "  CHECK(XYZ[%d] == Approx(%.12g));\n"%(idx, val)
      else:
        out += "  CHECK(ct.getYlm(%d) == Approx(%.12g));\n"%(idx, val)
    dx = diff(expr, Symbol('x'))
    dy = diff(expr, Symbol('y'))
    dz = diff(expr, Symbol('z'))
    dx_val =  dx.subs(slist).subs(rlist).evalf()
    dy_val =  dy.subs(slist).subs(rlist).evalf()
    dz_val =  dz.subs(slist).subs(rlist).evalf()
    if print_grad:
      if using_soa:
        out += "  CHECK(gr0[%d] == Approx(%.12g));\n"%(idx, dx_val)
        out += "  CHECK(gr1[%d] == Approx(%.12g));\n"%(idx, dy_val)
        out += "  CHECK(gr2[%d] == Approx(%.12g));\n"%(idx, dz_val)
      else:
        out += "  CHECK(ct.getGradYlm(%d)[0] == Approx(%.12g));\n"%(idx, dx_val)
        out += "  CHECK(ct.getGradYlm(%d)[1] == Approx(%.12g));\n"%(idx, dy_val)
        out += "  CHECK(ct.getGradYlm(%d)[2] == Approx(%.12g));\n"%(idx, dz_val)

    lap = diff(expr, Symbol('x'), 2) + diff(expr, Symbol('y'), 2) + diff(expr, Symbol('z'), 2)
    lap_val = lap.subs(slist).subs(rlist).evalf()
    if print_lap:
      if using_soa:
        out += "  CHECK(lap[%d] == Approx(%.12g));\n"%(idx, lap_val)
      else:
        out += "  CHECK(ct.getLaplYlm(%d) == Approx(%.12g));\n"%(idx, lap_val)


    if print_hess:
      out += '\n'
      axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
      for ii,si in enumerate(axis_syms):
        for jj,sj in enumerate(axis_syms):
          h_s = diff(diff(expr, si), sj)
          hess_val = h_s.subs(slist).subs(rlist).evalf()

          #print ii,jj,hess_val
          #if hess_val != 0:
          #  print "  CHECK(ct.getHessYlm(%d)(%d,%d) == Approx(%.12g));"%(idx, ii, jj, hess_val)
          if using_soa:
            if ii <= jj:
              out += "  CHECK(h%d%d[%d] == Approx(%.12g));\n"%(ii, jj, idx, hess_val)
          else:
            out += "  CHECK(ct.getHessYlm(%d)(%d,%d) == Approx(%.12g));\n"%(idx, ii, jj, hess_val)

      out += '\n'

    if print_ggg:
      out += '\n'
      axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
      for ii,si in enumerate(axis_syms):
        for jj,sj in enumerate(axis_syms):
          for kk,sk in enumerate(axis_syms):
            ggg_s = diff(diff(diff(expr, si), sj), sk)
            ggg_val = ggg_s.subs(slist).subs(rlist).evalf()

            out += "  CHECK(ct.getGGGYlm(%d)[%d](%d,%d) == Approx(%.12g));\n"%(idx, ii, jj, kk, ggg_val)

      out += '\n'

  return out

# A simple template replacement engine.
# Template items to be replaced start on a line with '%'.
def run_template(fname_in, fname_out, bodies):
  out = ''
  with open(fname_in, 'r') as f:
    for line in f:
      if line.startswith('%'):
        key = line.strip()[1:]
        if key in bodies:
          line = bodies[key]
        else:
          print 'Error, template item not found, key:',key, ' line = ',line
      out += line

  with open(fname_out, 'w') as f:
    f.write(out)


def create_test_full_cartesian_tensor():
  p = [1.3,1.2,-0.5]
    #compute_radial_values([1.3,1.2,-0.5])
  bodies = dict()
  bodies['test_evaluate'] = compute_radial_values(p, print_value=True)
  bodies['test_evaluate_all'] = compute_radial_values(p, print_value=True,
                                                         print_grad=True,
                                                         print_lap=True)
  bodies['test_evaluate_with_hessian'] = compute_radial_values(p, print_value=True,
                                                                  print_grad=True,
                                                                  print_hess=True)
  bodies['test_evaluate_with_third_deriv'] = compute_radial_values(p, print_value=True,
                                                                      print_grad=True,
                                                                      print_hess=True,
                                                                      print_ggg=True)
  bodies['test_evaluate_third_deriv_only'] = compute_radial_values(p, print_ggg=True)

  fname_in = 'test_full_cartesian_tensor.cpp.in'
  fname_out = 'test_full_cartesian_tensor.cpp'

  run_template(fname_in, fname_out, bodies)

def create_test_full_soa_cartesian_tensor():
  p = [1.3,1.2,-0.5]
    #compute_radial_values([1.3,1.2,-0.5])
  bodies = dict()
  bodies['test_evaluateV'] = compute_radial_values(p, print_value=True, using_soa=True)
  bodies['test_evaluateVGL'] = compute_radial_values(p, print_value=True,
                                                         print_grad=True,
                                                         print_lap=True,
                                                         using_soa=True)
  bodies['test_evaluateVGH'] = compute_radial_values(p, print_value=True,
                                                                  print_grad=True,
                                                                  print_hess=True,
                                                                  using_soa=True)

  fname_in = 'test_full_soa_cartesian_tensor.cpp.in'
  fname_out = 'test_full_soa_cartesian_tensor.cpp'

  run_template(fname_in, fname_out, bodies)




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
    #compute_radial_values([1.3,1.2,-0.5])

    # Create full test
    create_test_full_cartesian_tensor()
    create_test_full_soa_cartesian_tensor()

