
# Compute spline functions symbolically.  Also solve for the coefficients and evaluate spline
# functions at a few positions.

# Generates code for test_splines() in test_multi_spline.cpp

# This file can get slow when computing the 3D coefficients - run under Pypy to speed it up

from sympy import *

from bspline_funcs import create_spline


# Construct and solve the matrix equation for the coefficients
# This produces a minimal length coefficient array.  For uniformity of access,
#  the actual coefficient array has the first three values replicated onto the end.
def gen_coeffs_1D():
  xs = Symbol('x')

  Delta = Symbol('Delta', positive=True)
  nknots = 5

  c = IndexedBase('c',shape=(nknots+3))

  spline = create_spline(nknots, xs, c, Delta)

  m = nknots
  # Construct the matrix of coefficients and positions
  mat = Matrix.eye(m)
  for i in range(m):
    for j in range(m):
      subslist = {}
      for k in range(nknots+3):
        if k%m == i:
          subslist[c[k]] = 1
        else:
          subslist[c[k]] = 0

      subslist[xs] = j*Delta
      print 'subslist',subslist

      e = spline.subs(subslist)
      print j,i,e
      mat[j,i] = e

  #print mat
  for i in range(m):
    print mat.row(i)

  tpi = 2*S.Pi
  delta = 1.0/nknots
  rhs = []
  for i in range(m):
    x = delta*(i%(nknots))
    val = sin(tpi*x).evalf()
    print i%nknots,x,val
    rhs.append(val)

  b = Matrix(rhs)
  print 'Function values = ',b
  res = mat.LUsolve(b)
  print 'Coefficients = ',res


# Construct and solve the matrix equation for coefficients

def gen_coeffs_3D():
  xs = Symbol('x')
  ys = Symbol('y')
  zs = Symbol('z')

  Delta = Symbol('Delta', positive=True)
  nknots = 5

  cx = IndexedBase('cx',shape=(nknots+3))
  cy = IndexedBase('cy',shape=(nknots+3))
  cz = IndexedBase('cz',shape=(nknots+3))

  spline_x = create_spline(nknots, xs, cx, Delta)
  spline_y = create_spline(nknots, ys, cy, Delta)
  spline_z = create_spline(nknots, zs, cz, Delta)

  # Construct the matrix of coefficients and positions
  m = nknots
  mat = Matrix.eye(m**3)
  for ix in range(m):
    for iy in range(m):
      for iz in range(m):
        for jx in range(m):
          for jy in range(m):
            for jz in range(m):
              i = ix*m*m + iy*m + iz
              j = jx*m*m + jy*m + jz
              subslist_x = {}
              for k in range(nknots+3):
                if k%m == ix:
                  subslist_x[cx[k]] = 1
                else:
                  subslist_x[cx[k]] = 0

              subslist_x[xs] = jx*Delta
              #print 'subslist',subslist

              subslist_y = {}
              for k in range(nknots+3):
                if k%m == iy:
                  subslist_y[cy[k]] = 1
                else:
                  subslist_y[cy[k]] = 0

              subslist_y[ys] = jy*Delta

              subslist_z = {}
              for k in range(nknots+3):
                if k%m == iz:
                  subslist_z[cz[k]] = 1
                else:
                  subslist_z[cz[k]] = 0

              subslist_z[zs] = jz*Delta

              #print j,i,'X',spline_x
              #print j,i,'Y',spline_y
              #print j,i,'Z',spline_z
              e = spline_x.subs(subslist_x) * spline_y.subs(subslist_y) * spline_z.subs(subslist_z)
              #print j,i,e
              mat[j,i] = e
        print j,i,e

  #print mat
  for i in range(m):
    print mat.row(i)

  tpi = 2*S.Pi
  delta = 1.0/nknots
  rhs = []
  for ix in range(m):
    for iy in range(m):
      for iz in range(m):
        x = delta*(ix%(nknots))
        y = delta*(iy%(nknots))
        z = delta*(iz%(nknots))
        val = sin(tpi*x).evalf() + sin(3*tpi*y).evalf() + sin(4*tpi*z).evalf()
        i = ix*m*m + iy*m + iz
        rhs.append(val)

  b = Matrix(rhs)
  print 'b = ',b
  res = mat.LUsolve(b)
  mp = nknots + 3
  subslist = dict()
  #print 'coefficients =',res
  for ix in range(m):
    for iy in range(m):
      for iz in range(m):
        idx = ix*m*m + iy*m + iz
        j = ix*mp*mp + iy*mp + iz
        subslist[cx[ix]*cy[iy]*cz[iz]] = res[idx]

  with open('coefs.dat','w') as f:
    f.write('# Coefficients for %d x %d x %d spline\n'%(m,m,m))
    f.write('# idx  ix iy iz  coefficient\n')
    for ix in range(m):
      for iy in range(m):
        for iz in range(m):
          idx = ix*m*m + iy*m + iz
          f.write('%d %d %d %d %15.10g\n'%(idx, ix, iy, iz, res[idx]))

def read_coefs():
  coef_list = list()
  coef_map = dict()
  with open('coefs.dat','r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('#') or len(line) == 0:
        continue

      lc = line.split()
      idx = int(lc[0])
      ix = int(lc[1])
      iy = int(lc[2])
      iz = int(lc[3])
      coef_val = float(lc[4])

      coef_list.append(coef_val)
      assert(len(coef_list) == idx+1)

      coef_map[(ix,iy,iz)] = coef_val

  return coef_list, coef_map


# Evaluate the spline and derivatives at various points
# and generate code for test_splines() in test_multi_spline.cpp

def evaluate_spline_3D():
  xs = Symbol('x')
  ys = Symbol('y')
  zs = Symbol('z')

  Delta = Symbol('Delta', positive=True)
  nknots = 5

  delta = 1.0/nknots

  cx = IndexedBase('cx',shape=(nknots+3))
  cy = IndexedBase('cy',shape=(nknots+3))
  cz = IndexedBase('cz',shape=(nknots+3))

  spline_x = create_spline(nknots, xs, cx, Delta)
  spline_y = create_spline(nknots, ys, cy, Delta)
  spline_z = create_spline(nknots, zs, cz, Delta)

  coef_list, coef_map = read_coefs()

  m = nknots

  # Create a substitution list for the coefficients
  subslist = dict()
  for ix in range(m):
    for iy in range(m):
      for iz in range(m):
        subslist[cx[ix]*cy[iy]*cz[iz]] = coef_map[(ix,iy,iz)]

  e = (spline_x * spline_y * spline_z).subs(Delta, delta)


  print
  print '//  Code from here to the end of the function is generated by gen_bspline_values.py'
  print
  pos = {xs:0, ys:0, zs:0}
  print '  TinyVector<T,3> pos = {%g, %g, %g};'%(pos[xs], pos[ys], pos[zs])
  print
  print '  // symbolic value at pos = ',e.subs(pos)
  print
  print '  aligned_vector<T> v(npad);'
  print '  bs.evaluate(pos, v);'
  # the 'expand' after substituting pos is necessary for the triples of coefficients (cx*cy*cz) to
  # be formed so the substitution of coefficients works
  val = e.subs(pos).expand().subs(subslist)
  print '  CHECK(v[0] == Approx(%15.10g));'%(val)
  print

  print '  VectorSoaContainer<T,3> dv(npad);'
  print '  VectorSoaContainer<T,6> hess(npad); // 6 - number of unique hessian components'
  print '  bs.evaluate_vgh(pos, v, dv, hess);'


  ex = diff(e, xs).subs(pos).expand().subs(subslist)
  ey = diff(e, ys).subs(pos).expand().subs(subslist)
  ez = diff(e, zs).subs(pos).expand().subs(subslist)
  print '  // Gradient'
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(0,ex)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(1,ey)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(2,ez)

  print
  print '  // Hessian'
  print '  for (int i = 0; i < 6; i++) {'
  print '    CHECK(hess[0][i] == Approx(0.0));'
  print '  }'


  print
  pos1 = {xs:0.1, ys:0.2, zs:0.3}
  print '  pos = {%g, %g, %g};'%(pos1[xs], pos1[ys], pos1[zs])
  print '  bs.evaluate(pos, v);'
  print

  val = e.subs(pos1).expand().subs(subslist)
  print '  // Value'
  print '  CHECK(v[0] == Approx(%15.10g));'%(val)
  print
  print '  bs.evaluate_vgh(pos, v, dv, hess);'
  print '  // Value'
  print '  CHECK(v[0] == Approx(%15.10g));'%(val)
  ex = diff(e, xs).subs(pos1).expand().subs(subslist)
  ey = diff(e, ys).subs(pos1).expand().subs(subslist)
  ez = diff(e, zs).subs(pos1).expand().subs(subslist)
  print '  // Gradient'
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(0,ex)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(1,ey)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(2,ez)

  print '  // Hessian'
  for idx,(d1,d2) in enumerate([(xs,xs), (xs,ys), (xs,zs), (ys, ys), (ys, zs), (zs, zs)]):
      hess = diff(diff(e,d1), d2).subs(pos1).expand().subs(subslist)
      print '  CHECK(hess[0][%d] == Approx(%15.10g));'%(idx, hess)
      #print d1,d2,hess

  print
  print
  print '  VectorSoaContainer<T,3> lap(npad);'
  print '  bs.evaluate_vgl(pos, v, dv, lap);'

  lap = diff(e,xs,2) + diff(e,ys,2) + diff(e,zs,2)
  lap_val = lap.subs(pos1).expand().subs(subslist)
  print '  // Value'
  print '  CHECK(v[0] == Approx(%15.10g));'%(val)
  print '  // Gradient'
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(0,ex)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(1,ey)
  print '  CHECK(dv[0][%d] == Approx(%15.10g));'%(2,ez)
  print '  // Laplacian'
  print '  CHECK(lap[0][0] == Approx(%15.10g));'%(lap_val)


if __name__ == '__main__':
  #gen_coeffs_1D()

  #  Run gen_coeffs_3D to create the coefs.dat file.
  #    This function takes a few minutes.  Run under Pypy to speed it up.
  gen_coeffs_3D()
  #  The evaluate routine runs faster by reading the coefficients from the file.
  evaluate_spline_3D()
