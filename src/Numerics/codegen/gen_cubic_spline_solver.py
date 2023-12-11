from __future__ import print_function
import sys
from codegen_extras import *
from sympy import *
from sympy.codegen.ast import For,CodeBlock, Comment
from sympy.codegen.cnodes import void

sys.path.append("../tests")
from eqn_manip import *

# Generate cubic spline solver routine
# To use:
# - Run this script and save the output
# - Clean up some of the non-code text at the beginning
# - Run through clang-format
# - Put into SplineSolvers.h


# See this Jupyter notebook for more details:
# https://github.com/QMCPACK/qmc_algorithms/blob/master/Wavefunctions/CubicSplineSolver.ipynb


#
# Tridiagonal solver
# From Wikipedia : https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
#
# It would be good to derive these using Gaussian elimination, but for now they
# will be treated as given.
#
n = Symbol('n', integer=True)
i = Symbol('i', integer=True)
x = IndexedBase('x',shape=(n,))
dp = IndexedBase("d'",shape=(n,))
cp = IndexedBase("c'",shape=(n,))
a = IndexedBase("a",shape=(n,))
b = IndexedBase("b",shape=(n,))
c = IndexedBase("c",shape=(n,))
d = IndexedBase("d",shape=(n,))

# Use the C++ range 0,n-1
start = 0
end = n-1

# forward sweep
teq1 = Eq(cp[start], c[start]/b[start])
teq2 = Eq(dp[start], d[start]/b[start])
teq3 = Eq(dp[i],(d[i] - dp[i-1]*a[i])/ (b[i] - cp[i-1]*a[i]))
teq4 = Eq(cp[i],c[i]/(b[i] - cp[i-1]*a[i]))

# backward sweep
teq5 = Eq(x[end],dp[end])
teq6 = Eq(x[i],dp[i] - cp[i]*x[i+1])


#
# Cubic spline equation derivation
#


# Distance from the previous knot, for the case of uniform knot spacing
t = Symbol('t')

# Number of knots
n = Symbol('n', integer=True)

# Function values to intepolated at the knots
y = IndexedBase('y',shape=(n,))

# Coefficients of the spline function
a,b,c,d = [IndexedBase(s, shape=(n,)) for s in 'a b c d'.split()]

# Knot locations
x = IndexedBase('x',shape=(n,))

# Spacing between knots
L = IndexedBase('L',shape=(n,))   # L[i] = x[i+1] - x[i]

# Cubic spline equation
si = a[i] + b[i]*t + c[i]*t*t + d[i]*t**3
print(si)

# Value at knots (t=0)
sp1 = Eq(si.subs(t,0), y[i])

# Value at next knot
sp2 = Eq(si.subs(t,L[i]), y[i+1])

# Express the second derivative at the beginning of the interval in terms of E
E = IndexedBase('E',shape=(n,))
sp3 = Eq(E[i], diff(si,t,2).subs(t,0))

# Express the second derivative at the end of the interval in terms of E
sp4 = Eq(E[i+1], diff(si,t,2).subs(t,L[i]))

# Solve for spline coefficients in terms of E's
sln = solve([sp1,sp2,sp3,sp4], [a[i],b[i],c[i],d[i]])

# also for i+1
sln1 = {k.subs(i,i+1):v.subs(i,i+1) for k,v in sln.items()}

# Continuity of first derivatives at knots
# This will define the tridiagonal system to be solved
sp5 = Eq(diff(si,t).subs(t,L[i]), diff(si,t).subs(i, i+1).subs(t,0))

sp6 = sp5.subs(sln).subs(sln1)
sp7 = expand(sp6)
sp8 = divide_terms(sp7, [E[i],E[i+1],E[i+2]], [y[i],y[i+1],y[i+2]])
sp9 = mult_eqn(sp8,6)
print(sp9)

# The index 'i' used in the cubic spline equations is not the same 'i' used
# in the tridigonal solver.   Here we need to make them match.
# The first foundry condition will the equation at index at 0.
# Adjust the indexing on this equation so i=1 is the index of the first continuity interval match
sp9 = sp9.subs(i,i-1)

# Extract the three coefficients in each row for the general case
symlist = [E[i-1],E[i],E[i+1],E[i+2]]
coeff1 = get_coeff_for(sp9.lhs, E[i-1], symlist)
coeff2 = get_coeff_for(sp9.lhs, E[i], symlist)
coeff3 = get_coeff_for(sp9.lhs, E[i+1], symlist)


# Now get the coefficients for the boundary conditions (first row and last row)

# Natural BC
bc_natural_start = Eq(E[i].subs(i,0),0)
bc_natural_end = Eq(E[i].subs(i,end),0)

# The coefficients and RHS for this BC are pretty simple. but we will follow
# a deterministic path for derivation anyway.
bc_natural_start_coeff1 = get_coeff_for(bc_natural_start.lhs, E[start],[E[start]])
bc_natural_start_coeff2 = get_coeff_for(bc_natural_start.lhs, E[start+1],[E[start],E[start+1]])
bc_natural_end_coeff1 = get_coeff_for(bc_natural_end.lhs, E[end-1],[E[end]])
bc_natural_end_coeff2 = get_coeff_for(bc_natural_end.lhs, E[end],[E[end]])

# BC - first derivative specified at the beginning of the range
yp0 = Symbol('yp0')
eqbc1=Eq(diff(si,t).subs(t,0).subs(sln).subs(i,0), yp0)
eqbc1b = divide_terms(expand(eqbc1),[E[0],E[1]],[y[0],y[1],yp0])
eqbc1c = mult_eqn(eqbc1b, 6)
bc_firstd_start_coeff1 = get_coeff_for(eqbc1c.lhs, E[0], [E[0],E[1]])
bc_firstd_start_coeff2 = get_coeff_for(eqbc1c.lhs, E[1], [E[0],E[1]])


# For the general algorithm, the input parameters for the boundary conditions are
#  - first derivative, if value is less than cutoff
#  - second derivative is zero, if vlaue is greater than cutoff

bc_cutoff = 0.99e30

tbc_start_coeff1 = Piecewise((bc_firstd_start_coeff1, yp0 < bc_cutoff),(bc_natural_start_coeff1,True))
tbc_start_coeff2 = Piecewise((bc_firstd_start_coeff2, yp0 < bc_cutoff),(bc_natural_start_coeff2,True))

sym_bc_start_coeff1 = Symbol('bc_start1')
sym_bc_start_coeff2 = Symbol('bc_start2')
bc_eqs = [Eq(sym_bc_start_coeff1, tbc_start_coeff1)]
bc_eqs.append(Eq(sym_bc_start_coeff2, tbc_start_coeff2))

# BC - first derivative specified at the end of the range
ypn = Symbol('ypn')
eqbc2=Eq(diff(si,t).subs(t,L[end-1]).subs(sln).subs(i,end-1),ypn)
eqbc2b = divide_terms(expand(eqbc2),[E[end-1],E[end]],[y[end-1],y[end],ypn])
eqbc2c = mult_eqn(eqbc2b, 6)
bc_firstd_end_coeff1 = get_coeff_for(eqbc2c.lhs, E[end-1],[E[end-1],E[end]])
bc_firstd_end_coeff2 = get_coeff_for(eqbc2c.lhs, E[end],[E[end-1],E[end]])

# Create the conditional expression for the end BC
tbc_end_coeff1 = Piecewise((bc_firstd_end_coeff1, ypn < bc_cutoff),(bc_natural_end_coeff1, True))
sym_bc_end_coeff1 = Symbol('bc_end1')
bc_eqs.append(Eq(sym_bc_end_coeff1, tbc_end_coeff1))
tbc_end_coeff2 = Piecewise((bc_firstd_end_coeff2, ypn < bc_cutoff),(bc_natural_end_coeff2, True))
sym_bc_end_coeff2 = Symbol('bc_end2')
bc_eqs.append(Eq(sym_bc_end_coeff2, tbc_end_coeff2))

# conditional expressions for RHS for boundary conditions
rhs_start = Piecewise((eqbc1c.rhs,yp0 < bc_cutoff),(bc_natural_start.rhs,True))
rhs_end = Piecewise((eqbc2c.rhs, ypn < bc_cutoff), (bc_natural_end.rhs, True))

sym_rhs_start = Symbol('rhs_start')
sym_rhs_end = Symbol('rhs_end')
bc_eqs.append(Eq(sym_rhs_start, rhs_start))
bc_eqs.append(Eq(sym_rhs_end, rhs_end))


#
# Substitutions for the tridiagonal solver
#

subslist = {
    a[start] : 0,
    a[i] : coeff1,
    a[end] : sym_bc_end_coeff1,

    b[start] : sym_bc_start_coeff1,
    b[i] : coeff2,
    b[end] : sym_bc_end_coeff2,

    c[start] : sym_bc_start_coeff2,
    c[i] : coeff3,
    c[end] : 0,

    d[start] : sym_rhs_start,
    d[i] : sp9.rhs,
    d[end] : sym_rhs_end,
}

# Replace knot spacing with differences between knot locations
subsL = {
  L[i] : x[i+1] - x[i],
  L[i+1] : x[i+2] - x[i+1],
  L[i-1] : x[i] - x[i-1],
  L[start] : x[start+1]-x[start],
  L[start+1] : x[start+2]-x[start+1],
  L[end-1] : x[end] - x[end-1],
}


# Substitute into the tridiagonal solver
teq2b = teq2.subs(subslist).subs(subsL)
teq3b = simplify(teq3.subs(subslist).subs(subsL))
teq4b = teq4.subs(subslist).subs(subsL)
teq5b = Eq(teq5.lhs,teq5.rhs.subs(dp[end],teq3.rhs).subs(i,end).subs(subslist))


# Extract sub-expressions
subexpr, final_expr = cse([simplify(teq3b),simplify(teq4b)],symbols=numbered_symbols('z'))

# Substitute knot spacing into the boundary conditions
bc_eqs2 = [eq.subs(subsL) for eq in bc_eqs]


# Use temporary storage for cp, and reuse output vector for dp
tmp = IndexedBase('u',shape=(n,))
y2 = IndexedBase('y2',shape=(n,))
storage_subs = {cp:y2, dp:tmp}
#storage_subs = {}
teq1c = teq1.subs(subslist).subs(storage_subs)
teq2c = teq2b.subs(subslist).subs(storage_subs)
teq3c = final_expr[0].subs(storage_subs)
teq4c = final_expr[1].subs(storage_subs)
teq5c = teq5b.subs(storage_subs).subs(x,y2)
teq6c = teq6.subs(storage_subs).subs(x,y2)


#
# Code Generation
#

# Output will be a template function - this is the type
templateT = Type('T')


# forward sweep
fr = ARange(start+1,end,1)

body = []
for e in subexpr:
    body.append(Variable(e[0],type=templateT).as_Declaration(value=e[1].subs(storage_subs)))

body.append(convert_eq_to_assignment(teq3c))
body.append(convert_eq_to_assignment(teq4c))
loop1 = For(i,fr,body)


# backward sweep
br = ARangeClosedEnd(end-1,start,-1)
loop2 = For(i,br,[convert_eq_to_assignment(teq6c)])


tmp_init = VariableWithInit("n",tmp,type=Type("std::vector<T>")).as_Declaration()

bc_tmps = []
for e in bc_eqs2:
    bc_tmps.append(Variable(e.lhs, type=templateT).as_Declaration(value=e.rhs))

body = [tmp_init]
body.extend(bc_tmps)
body.extend([convert_eq_to_assignment(teq1c),
             convert_eq_to_assignment(teq2c),
             loop1,
              convert_eq_to_assignment(teq5c),
              loop2])

algo = CodeBlock(*body)

# Set up to create a template function
tx = Pointer(x,type=templateT)
ty = Pointer(y,type=templateT)
ty2 = Pointer(y2,type=templateT)
yp0_var = Variable('yp0',type=templateT)
ypn_var = Variable('ypn',type=templateT)

tf = TemplateFunctionDefinition(void, "CubicSplineSolve",[tx,ty,n,yp0_var,ypn_var,ty2],[templateT],algo)


ACP = ACodePrinter()
gen_comment = Comment("Generated by gen_cubic_spline_solver.py")
cb = CodeBlock(gen_comment, tf)
s = ACP.doprint(cb)
print(s)


