
from sympy import *

# Verify SymTrace in MultiBspline.hpp

# The first part of this script verifies the expression for the trace of the product of a
# symmetric matrix (H) and a general matrix (G).

# The second part creates body of the SymTrace test in test_multi_spline.cpp.

G = MatrixSymbol('G', 3, 3)
H = MatrixSymbol('H', 3, 3)

h1 = Matrix(H)
print 'H = ',h1
# Symmetrize H
h1 = h1.subs(H[1,0], H[0,1])
h1 = h1.subs(H[2,1], H[1,2])
h1 = h1.subs(H[2,0], H[0,2])
print 'Symmetrized H = ',h1
print

e =  Trace(h1*Matrix(G)).doit()
print 'Trace = ',e

h00 = Symbol('h00')
h01 = Symbol('h01')
h02 = Symbol('h02')
h11 = Symbol('h11')
h12 = Symbol('h12')
h22 = Symbol('h22')
print
e = e.subs(H[0,0], h00)
e = e.subs(H[0,1], h01)
e = e.subs(H[0,2], h02)
e = e.subs(H[1,1], h11)
e = e.subs(H[1,2], h12)
e = e.subs(H[2,2], h22)
print 'Trace = ',e
print
e2 =  e.collect([h00,h01,h02,h12,h11,h22])
print 'Trace =',e2

g = IndexedBase('g')
e2 = e2.subs(G[0,0], g[0])
e2 = e2.subs(G[0,1]+G[1,0], g[1])
e2 = e2.subs(G[0,2]+G[2,0], g[2])
e2 = e2.subs(G[1,1], g[1])
e2 = e2.subs(G[1,2]+G[2,1], g[4])
e2 = e2.subs(G[2,2], g[5])
print
print 'Replace with symmetrized G'
print 'Trace =',e2

vH = Matrix([[1.0, 2.0, 3.0],
             [2.0, 4.4, 1.1],
             [3.0, 1.1, 0.9]])
vG = Matrix([[0.1, 0.2, 0.3],
             [1.4, 2.3, 8.0],
             [0.9, 1.4, 2.3]])

tr = Trace(vH*vG).doit()
print 'Trace = ',tr


print
print 'Concrete values for unit test'
print '// Code from here to the end of the function generate by gen_trace.py'
print
print '  double h00 = %g;'%vH[0,0]
print '  double h01 = %g;'%vH[0,1]
print '  double h02 = %g;'%vH[0,2]
print '  double h11 = %g;'%vH[1,1]
print '  double h12 = %g;'%vH[1,2]
print '  double h22 = %g;'%vH[2,2]
print
print '  double gg[6] = {%g, %g, %g, %g, %g, %g};'%(vG[0,0], vG[0,1] + vG[1,0], vG[0,2] + vG[2,0], vG[1,1], vG[1,2] + vG[2,1], vG[2,2])
print
print '  double tr = SymTrace(h00, h01, h02, h11, h12, h22, gg);'
print '  CHECK(tr == Approx(%g));'%tr
print
