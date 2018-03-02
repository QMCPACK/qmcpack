
# Generate Bspline Jastrow values to test against

# Cut and paste the relevant part of the output into test_bspline_jastrow.cpp

from sympy import *
from collections import defaultdict
import numpy as np

def to_interval(ival):
    """Convert relational expression to an Interval"""
    min_val = None
    lower_open = False
    max_val = None
    upper_open = True
    if isinstance(ival, And):
        for rel in ival.args:
            if isinstance(rel, StrictGreaterThan):
                min_val = rel.args[1]
                #lower_open = True
            elif isinstance(rel, GreaterThan):
                min_val = rel.args[1]
                #lower_open = False
            elif isinstance(rel, StrictLessThan):
                max_val = rel.args[1]
                #upper_open = True
            elif isinstance(rel, LessThan):
                max_val = rel.args[1]
                #upper_open = False
            else:
                print('unhandled ',rel)

    if min_val == None or max_val == None:
        print('error',ival)
    return Interval(min_val, max_val, lower_open, upper_open)

# Transpose the interval and coefficients
#  Note that interval [0,1) has the polynomial coefficients found in the einspline code
#  The other intervals could be shifted, and they would also have the same polynomials
def transpose_interval_and_coefficients(sym_basis):
    cond_map = defaultdict(list)

    i1 = Interval(0,5, False, False) # interval for evaluation
    for idx, s0 in enumerate(sym_basis):
        for expr, cond in s0.args:
            if cond != True:
                i2 = to_interval(cond)
                if not i1.is_disjoint(i2):
                    cond_map[i2].append( (idx, expr) )
    return cond_map

# Create piecewise expression from the transposed intervals
# basis_map - map of interval to list of spline expressions for that interval
#         c - coefficient symbol (needs to allow indexing)
#        xs - symbol for the position variable ('x')
def recreate_piecewise(basis_map, c, xs):
    args = []
    for cond, exprs in basis_map.items():
        e = 0
        for idx, b in exprs:
            e += c[idx] * b
        args.append( (e, cond.as_relational(xs)))
    args.append((0, True))
    return Piecewise(*args)





def gen_bspline_jastrow(nknots, rcut_val, param, cusp_val):
  xs = Symbol('x')

  Delta = Symbol('Delta',positive=True)
  knots = [i*Delta for i in range(nknots)]
  print 'knots = ',knots
  all_knots = [i*Delta for i in range(-3,nknots+3)]
  rcut = (nknots-1)*Delta

  # Third-order bspline
  jastrow_sym_basis = bspline_basis_set(3, all_knots, xs)
  print("Number of basis functions = ",len(jastrow_sym_basis))
  #jastrow_sym_basis

  # Rearrange the basis and conditionals into a more useful form
  jastrow_cond_map = transpose_interval_and_coefficients(jastrow_sym_basis)
  c = IndexedBase('c',shape=(nknots+3))
  #c = MatrixSymbol('c',nknots+2,1)  # better for code-gen
  jastrow_spline = recreate_piecewise(jastrow_cond_map, c, xs)


  Delta_val = rcut_val*1.0/(nknots+1)
  #print 'Delta = ',Delta_val
  #print 'coeff size = ',nknots+4

  coeffs = np.zeros(nknots+4)
  coeffs[0] = -2*cusp_val*Delta_val + param[1]
  coeffs[1] = param[0]
  coeffs[2] = param[1]
  coeffs[3:-3] = param[2:]
  #print 'coeffs',coeffs

  deriv_jastrow_spline = diff(jastrow_spline, xs)
  deriv2_jastrow_spline = diff(jastrow_spline, xs, 2)

  vals = []
  for i in range(20):
    x = 0.6*i
    jv = jastrow_spline.subs({xs:x,Delta:Delta_val})
    jd = deriv_jastrow_spline.subs({xs:x,Delta:Delta_val})
    jdd = deriv2_jastrow_spline.subs({xs:x,Delta:Delta_val})
    #print jj
    subslist = dict()
    for i in range(12):
      subslist[c[i]] = coeffs[i]
    #print x,jv.subs(subslist),jd.subs(subslist),jdd.subs(subslist)
    vals.append((x,jv.subs(subslist),jd.subs(subslist),jdd.subs(subslist)))

 # Assumes
 # struct JValues
 # {
 #  double r;
 #  double u;
 #  double du;
 #  double ddu;
 # };
  tmpl = """
 const int N = {N};
 JValues Vals[N] = {{
   {values}
 }};
"""
  fmt_values = ',\n  '.join("{%.2f, %15.10g, %15.10g, %15.10g}"%(r,u,du,ddu) for r,u,du,ddu in vals)
  s = tmpl.format(N=len(vals), values=fmt_values)
  print s


# Generate output for these parameters

#<jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\"> \
#   <correlation rcut=\"10\" size=\"10\" speciesA=\"u\" speciesB=\"d\"> \
#      <coefficients id=\"ud\" type=\"Array\"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients> \
#    </correlation> \
#</jastrow> \

def gen_case_two_body():
  rcut_val = 10
  param = np.array([0.02904699284, -0.1004179, -0.1752703883, -0.2232576505, -0.2728029201, -0.3253286875, -0.3624525145, -0.3958223107, -0.4268582166, -0.4394531176])
  nknots = 10
  cusp_val = -0.5
  gen_bspline_jastrow(nknots, rcut_val, param, cusp_val)


#   <jastrow type=\"One-Body\" name=\"J1\" function=\"bspline\" source=\"ion0\" print=\"yes\"> \
#       <correlation elementType=\"C\" size=\"8\" cusp=\"0.0\"> \
#               <coefficients id=\"eC\" type=\"Array\"> \
#-0.2032153051 -0.1625595974 -0.143124599 -0.1216434956 -0.09919771951 -0.07111729038 \
#-0.04445345869 -0.02135082917 \
#               </coefficients> \
#            </correlation> \
#         </jastrow> \


def gen_case_one_body():
  rcut_val = 10
  param = np.array([-0.2032153051, -0.1625595974, -0.143124599, -0.1216434956, -0.09919771951, -0.07111729038, -0.04445345869, -0.02135082917])
  nknots = 8
  cusp_val = 0.0
  gen_bspline_jastrow(nknots, rcut_val, param, cusp_val)

if __name__ == '__main__':
  #gen_case_two_body()
  gen_case_one_body()

