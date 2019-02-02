from __future__ import print_function

from sympy import Eq

# Perform equation manipulations as one might use when working by hand to put
# equations into a particular form - e.g. moving all the terms of one type to
# one side.
# There is also a function to extract coefficients for a particular symbol,
# used for creating the matrix of coefficients from a set of equations.

# Used by scripts that symbolically derive equations (cubic splines, etc)

# Move symbols in sym_list from left hand side of equation to right hand side
def move_terms(eqn, sym_list):
    new_lhs = eqn.lhs
    new_rhs = eqn.rhs
    for sym in sym_list:
        c = eqn.lhs.coeff(sym)
        new_lhs = new_lhs - c*sym
        new_rhs = new_rhs - c*sym
    return Eq(new_lhs, new_rhs)

# Move symbols in sym_list from right hand side of equation to left hand side
def move_terms_left(eqn, sym_list):
    new_lhs = eqn.lhs
    new_rhs = eqn.rhs
    for sym in sym_list:
        c = eqn.rhs.coeff(sym)
        new_lhs = new_lhs - c*sym
        new_rhs = new_rhs - c*sym
    return Eq(new_lhs, new_rhs)

# Move all there terms in the symbol lists to the respective side of the equation
def divide_terms(eqn, sym_list_left, sym_list_right):
    #print 'start',eqn
    eqn1 = move_terms(eqn, sym_list_right)
    #print 'middle ',eqn1
    eqn2 = move_terms_left(eqn1, sym_list_left)
    return eqn2

# Multiply equation by term
def mult_eqn(eqn, e):
    return Eq(eqn.lhs*e, eqn.rhs*e)

# Extract coefficient from an expression for the symbol 'sym'.
# Works by setting all values in symlist to zero, except the target in 'sym'.
def get_coeff_for(expr, sym, symlist):
    # expression
    # symbol to get coefficient for
    # symlist - total list of symbols
    subslist = {}
    subslist[sym] = 1
    for s in symlist:
        if s != sym:
            subslist[s] = 0
    coeff = expr.subs(subslist)
    return coeff
