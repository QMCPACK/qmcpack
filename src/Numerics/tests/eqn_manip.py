from __future__ import print_function

from sympy import Eq

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

def divide_terms(eqn, sym_list_left, sym_list_right):
    #print 'start',eqn
    eqn1 = move_terms(eqn, sym_list_right)
    #print 'middle ',eqn1
    eqn2 = move_terms_left(eqn1, sym_list_left)
    return eqn2

# Multiply equation by term
def mult_eqn(eqn, e):
    return Eq(eqn.lhs*e, eqn.rhs*e)

# for all values other than the target, set to zero
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
