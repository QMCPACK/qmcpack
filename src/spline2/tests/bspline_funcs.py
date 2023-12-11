

from sympy import *
from collections import defaultdict

# Collection of routines for generating symbolic spline functions
# Some of these routines should be common with QMCWavefunctions/tests/gen_bspline_jastrow.py

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


# Get the values corresponding to the interval starting at 0
# basis_map - map of interval to list of spline expressions for that interval
def get_base_interval(basis_map):
    for cond, exprs in basis_map.items():
        if cond.start == 0:
          return exprs
    return None


# Create a set of spline functions
def create_spline(nknots, xs, c, Delta):
  Delta = Symbol('Delta', positive=True)
  all_knots = [i*Delta for i in range(-3, nknots+3)]

  # Third-order bspline
  sym_basis = bspline_basis_set(3, all_knots, xs)
  cond_map = transpose_interval_and_coefficients(sym_basis)
  spline = recreate_piecewise(cond_map, c, xs)

  return spline

