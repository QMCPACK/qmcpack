
from __future__ import print_function
from sympy import *
from eqn_manip import *

# Solve for cubic spline coefficients using a straightforward (but inefficient) derivation
#  from the defining equations.

# Generates unit tests for the spline coefficient solver (CubicSplineSolve in SplineSolvers.h)
# and the evaluation routines (OneDimCubicSpline in OneDimCubicSpline.h)
# Run this script when these unit tests need updating or expanding.
#   (The computation for the 4-knot case may take a few minutes)
# Put the output of this script into test_one_dim_cubic_spline.cpp, and then run through
# clang-format to clean it up.

# See for more explanation of the derivation, see
# https://github.com/QMCPACK/qmc_algorithms/blob/master/Wavefunctions/Cubic%20Splines%20Basic.ipynb



# Symbols useful enough to be global

# Distance to knot, for non-uniform knots
t = IndexedBase('t')

y = IndexedBase('y')
x = IndexedBase('x')


# Create the equations that define cubic splines and solve by brute force
def create_solution_for_val(nknots, naturalBC=(True, True), firstDeriv=(0.0, 0.0)):

  n = Symbol('n', integer=True)
  i = Symbol('i', integer=True)

  a,b,c,d = [IndexedBase(s) for s in 'a b c d'.split()]
  # Non-uniform knots
  si = a[i] + b[i]*t[i] + c[i]*t[i]*t[i] + d[i]*t[i]**3

  # Value at knots (t=0)
  eq1 = Eq(si.subs(t[i],0), y[i])

  # Value at knots (t=1)
  eq2 = Eq(si.subs(t[i],x[i+1]-x[i]), y[i+1])

  # Continuity of first derivatives at knots
  dsi = diff(si, t[i])
  eq3 = Eq(dsi.subs(t[i],x[i+1]-x[i]), dsi.subs(i, i+1).subs(t[i+1],0))

  # Continuity of second derivatives at knots
  d2si = diff(si, t[i], 2)
  eq4 = Eq(d2si.subs(t[i],x[i+1]-x[i]).subs(t[0],x[0]), d2si.subs(i, i+1).subs(t[i+1],0))


  sym_lhs = [a[i],b[i],c[i],d[i],a[i+1],b[i+1],c[i+1],d[i+1]]
  sym_rhs = [y[i]]
  eq3b = divide_terms(eq3, sym_lhs, sym_rhs)
  eq4b = divide_terms(eq4, sym_lhs, sym_rhs)

  if naturalBC[0]:
    # Natural BC (second deriv is zero at both ends)
    eq5 = Eq(d2si.subs(i,0).subs(t[0],0), 0)
  else:
    # Specified first derivatives at boundaries
    eq5 = Eq(dsi.subs(i,0).subs(t[0],0), firstDeriv[0])


  if naturalBC[1]:
    eq6 = Eq(d2si.subs(t[i],x[i+1]-x[i]).subs(i,n-1), 0)
  else:
    eq6 = Eq(dsi.subs(t[i],x[i+1]-x[i]).subs(i,n-1), firstDeriv[1])



  nval = nknots - 1
  neqn = 4*nval
  m = Matrix.eye(neqn, neqn)
  rhs = Matrix.eye(neqn,1)

  symlist = list()
  for idx in range(nval):
      symlist.append(a[idx])
      symlist.append(b[idx])
      symlist.append(c[idx])
      symlist.append(d[idx])
  #print(symlist)

  for idx,sym in enumerate(symlist):
      m[0,idx] = get_coeff_for(eq5.lhs.subs(i,0), sym, symlist)
      m[1,idx] = get_coeff_for(eq1.lhs.subs(i,0), sym, symlist)
      m[2,idx] = get_coeff_for(eq2.lhs.subs(i,0), sym, symlist)
  rhs[0] = eq5.rhs.subs(i,0)
  rhs[1] = eq1.rhs.subs(i,0)
  rhs[2] = eq2.rhs.subs(i,0)

  jdx = 3
  for nv in range(1,nval):
      for idx,sym in enumerate(symlist):
          m[jdx+0,idx] = get_coeff_for(eq1.lhs.subs(i,nv), sym, symlist)
          m[jdx+1,idx] = get_coeff_for(eq2.lhs.subs(i,nv), sym, symlist)
          m[jdx+2,idx] = get_coeff_for(eq3b.lhs.subs(i,nv-1), sym, symlist)
          m[jdx+3,idx] = get_coeff_for(eq4b.lhs.subs(i,nv-1), sym, symlist)
      rhs[jdx+0] = eq1.rhs.subs(i,nv)
      rhs[jdx+1] = eq2.rhs.subs(i,nv)
      rhs[jdx+2] = eq3b.rhs.subs(i,nv)
      rhs[jdx+3] = eq4b.rhs.subs(i,nv)
      jdx += 4

  for idx,sym in enumerate(symlist):
      m[jdx,idx] = get_coeff_for(eq6.lhs.subs(n,nval), sym , symlist)
  rhs[jdx] = eq6.rhs.subs(n,nval)

  #for idx in range(4*nval):
  #   print('m = ',m.row(idx))
  #print('m = ',m)

  sln=m.LUsolve(rhs)

  subslist={}
  for idx in range(nval):
      subslist[a[idx]] = sln[4*idx + 0]
      subslist[b[idx]] = sln[4*idx + 1]
      subslist[c[idx]] = sln[4*idx + 2]
      subslist[d[idx]] = sln[4*idx + 3]

  spline = dict()
  for idx in range(nval):
    spline[idx] = si.subs(i,idx).subs(subslist)
  return spline


# Create C++ brace initializer list from a python list
def convert_to_brace_list(xvals, vertical_layout=False, constructor=None):
  joiner = ","
  if vertical_layout:
    joiner = ",\n"
  s = "{"
  if constructor is None:
    s += joiner.join([str(x) for x in xvals])
  else:
    # This assumes the representation (str(x)) already encloses the value in parentheses.
    # This is true if the entries in xvals are tuples.
    s += joiner.join([constructor + str(x)  for x in xvals])
  s+= "}"
  return s;

# --------------------------------------------------------------
# Tests of the coefficients from the cubic spline solver routine
# --------------------------------------------------------------

spline_solver_template ="""
// Generated from gen_cubic_spline.py
TEST_CASE("spline_function_{test_case}", "[numerics]")
{{
  const int n = {n};
  double x[n] = {xvals};
  double y[n] = {yvals};
  double y2[n];

  CubicSplineSolve(x, y, n, {leftDeriv}, {rightDeriv}, y2);

{y2_checks}

}}
"""


# Test the coefficient values produced by the cubic spline solver
def generate_cubic_spline_coeff_test_case(n, xvals, yvals, naturalBC=(True,True),
                                          firstDeriv=(0.0, 0.0), test_idx=0):
  sln = create_solution_for_val(n, naturalBC=naturalBC, firstDeriv=firstDeriv)

  x_init_expr = convert_to_brace_list(xvals)
  y_init_expr = convert_to_brace_list(yvals)

  ysubs = {y[i]:yv for i,yv in enumerate(yvals)}
  xsubs = {x[i]:xv for i,xv in enumerate(xvals)}

  nd = 1e33 # signal for natural derivative

  checks_str = ""
  for i in range(n-1):
    si2 = diff(sln[i], t[i], 2)
    sval = si2.subs(xsubs).subs(ysubs).subs(t[i], 0.0)
    #print('sval = ',sval)
    check = "  CHECK(y2[%d] == Approx(%g));"%(i, sval)
    checks_str += check + "\n"

  # Last second derivative must be obtained from the end of the last interval
  si2 = diff(sln[n-2], t[n-2], 2)
  sval = si2.subs(xsubs).subs(ysubs).subs(t[n-2], xvals[n-1]-xvals[n-2])
  #print('last sval = ',sval)
  check = "  CHECK(y2[%d] == Approx(%g));"%(n-1, sval)
  checks_str += check + "\n"

  leftDeriv = nd if naturalBC[0] else firstDeriv[0]
  rightDeriv = nd if naturalBC[1] else firstDeriv[1]
  out = spline_solver_template.format(test_case=test_idx, n=n, xvals=x_init_expr, yvals=y_init_expr,
                                      leftDeriv=leftDeriv, rightDeriv=rightDeriv,y2_checks=checks_str)
  print(out)


# Create multiple test cases with different number of knots and different boundary conditions
def generate_cubic_spline_coeff_test_cases():
  xvals = [0.0, 1.0, 2.0]
  yvals = [1.0, 2.0, 1.5]
  generate_cubic_spline_coeff_test_case(len(xvals), xvals, yvals, test_idx=1)
  generate_cubic_spline_coeff_test_case(len(xvals), xvals, yvals, naturalBC=(True,False), firstDeriv=(0.0, 1.0), test_idx=2)
  generate_cubic_spline_coeff_test_case(len(xvals), xvals, yvals, naturalBC=(False,False), firstDeriv=(1.0, 2.0), test_idx=3)

  xvals = [0.0, 1.2, 2.4, 3.0]
  yvals = [1.0, 2.0, 1.5, 1.8]
  generate_cubic_spline_coeff_test_case(len(xvals), xvals, yvals, test_idx=4)
  generate_cubic_spline_coeff_test_case(len(xvals), xvals, yvals, naturalBC=(False,False), firstDeriv=(0.0, 1.0), test_idx=5)



# ---------------------------------------------
# Tests of the cubic spline evaluation routines
# ---------------------------------------------

spline_evaluate_data_structure = """
// Structure for holding value, first and second derivatives
struct D2U
{
  D2U(double val_, double du_, double d2u_) : val(val_), du(du_), d2u(d2u_) {}
  double val;
  double du;
  double d2u;
};
"""


spline_evaluate_template ="""
// Generated from gen_cubic_spline.py
TEST_CASE("one_dim_cubic_spline_{test_case}", "[numerics]")
{{

  const int n = {n};
  std::vector<double> yvals = {yvals};

  LinearGrid<double> grid;
  grid.set(0.0, 2.0, n);

  OneDimCubicSpline<double> cubic_spline(&grid, yvals);

  int imin = 0;
  int imax = {n}-1;

  double yp0 = {leftDeriv};
  double ypn = {rightDeriv};

  cubic_spline.spline(imin, yp0, imax, ypn);

  std::vector<double> check_xvals = {check_xvals};
  std::vector<double> check_yvals;
  std::vector<D2U> check_yvals_d2u;

  for (int i = 0; i < check_xvals.size(); i++) {{
    double r = check_xvals[i];
    double val = cubic_spline.splint(r);
    check_yvals.push_back(val);

    double du, d2u;
    double val2 = cubic_spline.splint(r, du, d2u);
    check_yvals_d2u.push_back(D2U(val2, du, d2u));

    //std::cout << i << " r = " << r << " val = " << val << " " << check_yvals[i] << std::endl;
  }}

{output_check}

}}
"""

def generate_cubic_spline_evaluation_tests():
  xvals = [0.0, 1.0, 2.0]
  yvals = [1.0, 2.0, 1.5]
  leftDeriv = 1.0
  rightDeriv = 2.0

  n = len(xvals)

  sln = create_solution_for_val(n, naturalBC=(False, False), firstDeriv=(leftDeriv, rightDeriv))

  ysubs = {y[i]:yv for i,yv in enumerate(yvals)}
  xsubs = {x[i]:xv for i,xv in enumerate(xvals)}

  checks_str = ""
  for i in range(n-1):
    si2 = diff(sln[i], t[i], 2)
    sval = si2.subs(xsubs).subs(ysubs).subs(t[i], 0.0)
    #print('2nd deriv %i %g'%(i,sval))

  # Shrink the range slightly so the last point doesn't fall exactly on the boundary
  grid_range = xvals[-1] - xvals[0] - 0.0000000001
  ncheck = 6
  delta = grid_range/(ncheck-1)
  check_xvals = []
  check_yvals = []
  check_y2vals = []
  check_y3vals = []
  for i in range(ncheck):
    rval = i*delta + xvals[0]
    check_xvals.append(rval)
    # Find the correct spline interval
    tval = 0.0
    for idx in range(n-1):
      tval = rval - xvals[idx]
      if rval >= xvals[idx] and rval < xvals[idx+1]:
        break

    # Compute spline value and derivatives
    val = sln[idx].subs(xsubs).subs(ysubs).subs(t[idx], tval)
    dval = diff(sln[idx],t[idx]).subs(xsubs).subs(ysubs).subs(t[idx], tval)
    d2val = diff(sln[idx],t[idx],2).subs(xsubs).subs(ysubs).subs(t[idx], tval)

    check_yvals.append(val)
    check_y2vals.append((val,dval,d2val))
    #print(rval,val)

  y_expr = convert_to_brace_list(yvals)
  x_check_expr = convert_to_brace_list(check_xvals, vertical_layout=True)
  #y_check_expr = convert_to_brace_list(check_yvals, vertical_layout=True)
  #y2_check_expr = convert_to_brace_list(check_y2vals, vertical_layout=True, constructor='D2U')

  output_check = ''
  for i in range(ncheck):
    output_check += "  CHECK(check_yvals[%d] == Approx(%12.8g));\n"%(i,check_yvals[i])

  output_check += '\n'

  for i in range(ncheck):
    output_check += "  CHECK(check_yvals_d2u[%d].val == Approx(%12.8g));\n"%(i,check_y2vals[i][0])
    output_check += "  CHECK(check_yvals_d2u[%d].du == Approx(%12.8g));\n"%(i,check_y2vals[i][1])
    output_check += "  CHECK(check_yvals_d2u[%d].d2u == Approx(%12.8g));\n"%(i,check_y2vals[i][2])
  output_check += '\n'

  test_idx = 1
  out = spline_evaluate_template.format(test_case=test_idx, n=n, yvals=y_expr,
                                        check_xvals=x_check_expr,
                                      leftDeriv=leftDeriv, rightDeriv=rightDeriv, output_check=output_check)

  print(spline_evaluate_data_structure)
  print("\n")
  print(out)

if __name__ == '__main__':
  # Put these in test_one_dim_cubic_spline.cpp, and then run through clang-format
  generate_cubic_spline_coeff_test_cases()
  generate_cubic_spline_evaluation_tests()
