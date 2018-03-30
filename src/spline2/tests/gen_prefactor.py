
# Generate spline prefactors to check code in MultiBspline.hpp


from sympy import *

from bspline_funcs import  transpose_interval_and_coefficients, get_base_interval


def gen_prefactor():
  xs = Symbol('x')

  Delta = Symbol('Delta', positive=True)
  nknots = 2
  all_knots = [i*Delta for i in range(-3, nknots+3)]

  # Third-order bspline
  sym_basis = bspline_basis_set(3, all_knots, xs)
  print("Number of basis functions = ",len(sym_basis))

  cond_map = transpose_interval_and_coefficients(sym_basis)

  spline_exprs = get_base_interval(cond_map)
  spline_exprs = [(idx,s.subs(Delta, 1)) for idx,s in spline_exprs]


  # For values
  xval = 0.3
  print
  print '  tx = %g'%xval
  for idx,s in spline_exprs:
    print '  REQUIRE(a[%d] == Approx(%g));'%(idx, s.subs(xs,xval))

  # For values, first and second derivatives
  print '  tx = %g'%xval
  for idx,s in spline_exprs:
    print '  REQUIRE(a[%d] == Approx(%g));'%(idx, s.subs(xs,xval))
    print '  REQUIRE(da[%d] == Approx(%g));'%(idx, diff(s,xs).subs(xs,xval))
    print '  REQUIRE(d2a[%d] == Approx(%g));'%(idx, diff(s,xs,2).subs(xs,xval))

if __name__ == '__main__':
  gen_prefactor()
