"""Collection of routines for generating symbolic spline functions

Some of these routines should be common with QMCWavefunctions/tests/gen_bspline_jastrow.py
"""

from collections import defaultdict

from sympy import (
    And,
    Expr,
    GreaterThan,
    IndexedBase,
    Interval,
    LessThan,
    Piecewise,
    StrictGreaterThan,
    StrictLessThan,
    Symbol,
    bspline_basis_set,
)


def to_interval(ival: And) -> Interval | None:
    """Convert relational expression to an Interval."""
    if not isinstance(ival, And):
        print("error", ival)

    min_val    = None
    lower_open = False
    max_val    = None
    upper_open = True
    for rel in ival.args:
        match rel:
            case StrictGreaterThan():
                min_val = rel.args[1]
                #lower_open = True
            case GreaterThan():
                min_val = rel.args[1]
                #lower_open = False
            case StrictLessThan():
                max_val = rel.args[1]
                #upper_open = True
            case LessThan():
                max_val = rel.args[1]
                #upper_open = False
            case _:
                print("unhandled ", rel)

    if min_val is None or max_val is None:
        print("error", ival)
    return Interval(min_val, max_val, lower_open, upper_open)


def transpose_interval_and_coefficients(
    sym_basis: list[Piecewise]
) -> defaultdict[Interval, tuple[int, Expr] | list]:
    """Transpose the interval and coefficients.

    Parameters
    ----------
    sym_basis : list of Piecewise
        3rd-order bspline basis.
        (output of ``sympy.bspline_basis_set(d, knots, x)`` with ``d=3``)

    Notes
    -----
    The interval [0,1) has the polynomial coefficients found in the
    einspline code. The other intervals could be shifted, and they would
    also have the same polynomials.
    """
    cond_map = defaultdict(list)

    i1 = Interval(0, 5, False, False) # Interval for evaluation
    for idx, s0 in enumerate(sym_basis):
        for expr, cond in s0.args:
            if cond == True:
                continue

            i2 = to_interval(cond)
            if not i1.is_disjoint(i2):
                cond_map[i2].append((idx, expr))
    return cond_map


def recreate_piecewise(basis_map: dict, c: IndexedBase, xs: Symbol) -> Piecewise:
    """Create piecewise expression from the transposed intervals.

    Parameters
    ----------
    basis_map : dict
        Map of interval to list of spline expressions for that interval.
    c : IndexedBase
        Coefficient symbol (needs to allow indexing).
    xs : Symbol
        Symbol for the position variable ('x').

    See Also
    --------
    transpose_interval_and_coefficients
    """

    args = []
    for cond, exprs in basis_map.items():
        e = 0
        for idx, b in exprs:
            e += c[idx] * b
        args.append((e, cond.as_relational(xs)))
    args.append((0, True))
    return Piecewise(*args)


def get_base_interval(
    basis_map: dict[Interval, tuple[int, Expr]]
) -> tuple[int, Expr] | None:
    """Get the values corresponding to the interval starting at 0.

    Parameters
    ----------
    basis_map : dict
        Map of interval to list of spline expressions for that interval.
    """
    for cond, exprs in basis_map.items():
        if cond.start == 0:
          return exprs
    return None


def create_spline(
    nknots: int,
    xs: Symbol,
    c: IndexedBase,
    Delta: Symbol | None = None,
) -> Piecewise:
    """Create a set of spline functions."""
    if Delta is None:
        Delta = Symbol('Delta', positive=True)
    all_knots = [i*Delta for i in range(-3, nknots+3)]

    # Third-order bspline
    sym_basis = bspline_basis_set(3, all_knots, xs)
    cond_map = transpose_interval_and_coefficients(sym_basis)
    spline = recreate_piecewise(cond_map, c, xs)

    return spline
