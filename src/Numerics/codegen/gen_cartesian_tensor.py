
# Generate angular tensors using symbolic expressions for the GTO's and derivatives

# There are two steps to the code generation, and only the second is automated.
#  1. Use read_order.py to generate:
#      A. A python version of get_ijk, which is pasted into this file.
#      B. A C++ version of get_ABC, which is pasted into CartesianTensor.h.in
#      C. Repeat B for a version of get_ABC to be pasted into SoaCartesianTensor.h.in
#  2. Run create_cartesian_tensor.py
#      A. Copy CartesianTensor.h up on directory (to src/Numerics)
#      B. Copy SoaCartesianTensor.h up to src/QMCWaveFunctions/LCAO)
#      C. Apply clang-format on modified source files.


from collections import namedtuple, defaultdict
from sympy import *

# See the GaussianOrbitals notebook in the qmc_algorithms repo for more explanation,
#  especially about the normalization.

def create_gto_symbolic():
    pre = sympify('x**i * y**j * z**k * exp(-alpha *r**2)')
    return pre

# this function generated from 'read_order.py', should match GAMESS order
def get_ijk():
  ijk = []
  # S
  ijk.append( (0,0,0,"S") )
  # P
  ijk.append( (1,0,0,"X") )
  ijk.append( (0,1,0,"Y") )
  ijk.append( (0,0,1,"Z") )
  # D
  ijk.append( (2,0,0,"XX") )
  ijk.append( (0,2,0,"YY") )
  ijk.append( (0,0,2,"ZZ") )
  ijk.append( (1,1,0,"XY") )
  ijk.append( (1,0,1,"XZ") )
  ijk.append( (0,1,1,"YZ") )
  # F
  ijk.append( (3,0,0,"XXX") )
  ijk.append( (0,3,0,"YYY") )
  ijk.append( (0,0,3,"ZZZ") )
  ijk.append( (2,1,0,"XXY") )
  ijk.append( (2,0,1,"XXZ") )
  ijk.append( (1,2,0,"YYX") )
  ijk.append( (0,2,1,"YYZ") )
  ijk.append( (1,0,2,"ZZX") )
  ijk.append( (0,1,2,"ZZY") )
  ijk.append( (1,1,1,"XYZ") )
  # G
  ijk.append( (4,0,0,"XXXX") )
  ijk.append( (0,4,0,"YYYY") )
  ijk.append( (0,0,4,"ZZZZ") )
  ijk.append( (3,1,0,"XXXY") )
  ijk.append( (3,0,1,"XXXZ") )
  ijk.append( (1,3,0,"YYYX") )
  ijk.append( (0,3,1,"YYYZ") )
  ijk.append( (1,0,3,"ZZZX") )
  ijk.append( (0,1,3,"ZZZY") )
  ijk.append( (2,2,0,"XXYY") )
  ijk.append( (2,0,2,"XXZZ") )
  ijk.append( (0,2,2,"YYZZ") )
  ijk.append( (2,1,1,"XXYZ") )
  ijk.append( (1,2,1,"YYXZ") )
  ijk.append( (1,1,2,"ZZXY") )
  # H
  ijk.append( (5,0,0,"XXXXX") )
  ijk.append( (0,5,0,"YYYYY") )
  ijk.append( (0,0,5,"ZZZZZ") )
  ijk.append( (4,1,0,"XXXXY") )
  ijk.append( (4,0,1,"XXXXZ") )
  ijk.append( (1,4,0,"YYYYX") )
  ijk.append( (0,4,1,"YYYYZ") )
  ijk.append( (1,0,4,"ZZZZX") )
  ijk.append( (0,1,4,"ZZZZY") )
  ijk.append( (3,2,0,"XXXYY") )
  ijk.append( (3,0,2,"XXXZZ") )
  ijk.append( (2,3,0,"YYYXX") )
  ijk.append( (0,3,2,"YYYZZ") )
  ijk.append( (2,0,3,"ZZZXX") )
  ijk.append( (0,2,3,"ZZZYY") )
  ijk.append( (3,1,1,"XXXYZ") )
  ijk.append( (1,3,1,"YYYXZ") )
  ijk.append( (1,1,3,"ZZZXY") )
  ijk.append( (2,2,1,"XXYYZ") )
  ijk.append( (2,1,2,"XXZZY") )
  ijk.append( (1,2,2,"YYZZX") )
  # I
  ijk.append( (6,0,0,"X6") )
  ijk.append( (0,6,0,"Y6") )
  ijk.append( (0,0,6,"Z6") )
  ijk.append( (5,1,0,"X5Y") )
  ijk.append( (5,0,1,"X5Z") )
  ijk.append( (1,5,0,"Y5X") )
  ijk.append( (0,5,1,"Y5Z") )
  ijk.append( (1,0,5,"Z5X") )
  ijk.append( (0,1,5,"Z5Y") )
  ijk.append( (4,2,0,"X4Y2") )
  ijk.append( (4,0,2,"X4Z2") )
  ijk.append( (2,4,0,"Y4X2") )
  ijk.append( (0,4,2,"Y4Z2") )
  ijk.append( (2,0,4,"Z4X2") )
  ijk.append( (0,2,4,"Z4Y2") )
  ijk.append( (4,1,1,"X4YZ") )
  ijk.append( (1,4,1,"Y4XZ") )
  ijk.append( (1,1,4,"Z4XY") )
  ijk.append( (3,3,0,"X3Y3") )
  ijk.append( (3,0,3,"X3Z3") )
  ijk.append( (0,3,3,"Y3Z3") )
  ijk.append( (3,2,1,"X3Y2Z") )
  ijk.append( (3,1,2,"X3Z2Y") )
  ijk.append( (2,3,1,"Y3X2Z") )
  ijk.append( (1,3,2,"Y3Z2X") )
  ijk.append( (2,1,3,"Z3X2Y") )
  ijk.append( (1,2,3,"Z3Y2X") )
  ijk.append( (2,2,2,"X2Y2Z2") )

  return ijk

#  Input is a list of i,j,k,s ( descriptive string )
#  Output is a list that adds the maximum L (sum of i,j,k) up to that point.
def gen_lmax(ijk_list):
  list_with_lmax = []
  lmax = -1
  for i,j,k,s in ijk_list:
    current_l = i+j+k
    if current_l > lmax:
      lmax = current_l
    list_with_lmax.append( (i,j,k,s,lmax) )
  return list_with_lmax


# Replace powers with pre-computed values.
# One reason for this replacement is the current code prints powers as '**' and so is not even valid C++.
#  Using the C printer and printing as 'pow' could lead to poor performance.
# This does a manual common subexpression elimination.  Could imagine using the cse module to do more
# automatically.  The tricky part might be only computing pieces needed for lmax.
def replace_common_subexpressions(expr, slist=None):
  x,y,z = symbols('x y z')
  x2,y2,z2 = symbols('x2 y2 z2')
  x3,y3,z3 = symbols('x3 y3 z3')
  x4,y4,z4 = symbols('x4 y4 z4')
  x5,y5,z5 = symbols('x5 y5 z5')

  rlist0 = {x**6:x*x5, y**6:y*y5, z**6:z*z5}
  rlist1 = {x**5:x*x4, y**5:y*y4, z**5:z*z4}
  rlist2 = {x**4:x4, y**4:y4, z**4: z4}
  rlist3 = {x**3:x3, y**3:y3, z**3:z3}
  rlist4 = {x**2:x2, y**2:y2, z**2:z2}

  if slist:
    expr = expr.subs(slist)

  return expr.subs(rlist0).subs(rlist1).subs(rlist2).subs(rlist3).subs(rlist4)

def gen_evaluate():
  out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluate(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  value_type x4=x3*x, y4=y3*y, z4=z3*z;
  value_type x5=x4*x, y5=y4*y, z5=z4*z;
  switch(Lmax)
  {
%s
  }
  for (int i=0; i<XYZ.size(); i++)
    XYZ[i]*= NormFactor[i];
}
"""
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}
    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s; // %s\n"%(idx,val,s)

  return out_str%body_str

def gen_evaluate_all():
  out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateAll(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  value_type x4=x3*x, y4=y3*y, z4=z3*z;
  value_type x5=x4*x, y5=y4*y, z5=z4*z;
  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
    gradXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
    laplXYZ[i]=0.0;

  switch(Lmax)
  {
%s
  }
  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    laplXYZ[i]*= NormFactor[i];

}
"""
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    # Compute derivatives symbolically
    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)
    lap = diff(gto_s, x, 2) + diff(gto_s, y, 2) + diff(gto_s, z, 2)

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}
    val = replace_common_subexpressions(gto_s, slist)
    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    lap_val = replace_common_subexpressions(lap, slist)


    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    if gx != 0:
      body_str += "    gradXYZ[%d][0] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gradXYZ[%d][1] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gradXYZ[%d][2] = %s;\n"%(idx,gz)

    if lap_val != 0:
      body_str += "    laplXYZ[%d] = %s;\n"%(idx,lap_val)

  return out_str%body_str

def gen_evaluate_with_hessian():
  out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateWithHessian(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  value_type x4=x3*x, y4=y3*y, z4=z3*z;
  value_type x5=x4*x, y5=y4*y, z5=z4*z;
  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
    gradXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
    hessXYZ[i]=0.0;

  switch(Lmax)
  {
%s
  }
  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    hessXYZ[i]*= NormFactor[i];
}
"""

  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}

    # Compute derivatives symbolically
    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)

    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    if gx != 0:
      body_str += "    gradXYZ[%d][0] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gradXYZ[%d][1] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gradXYZ[%d][2] = %s;\n"%(idx,gz)

    axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        # Compute Hessian elements symbolically
        h_s = diff(diff(gto_s, si), sj)
        hess_val = replace_common_subexpressions(h_s, slist)
        if hess_val != 0:
          body_str += "    hessXYZ[%d](%d,%d) = %s;\n"%(idx,ii,jj,hess_val)

  return out_str%body_str

def gen_evaluate_with_third_deriv():
  out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateWithThirdDeriv(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  value_type x4=x3*x, y4=y3*y, z4=z3*z;
  value_type x5=x4*x, y5=y4*y, z5=z4*z;

  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
    gradXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
    hessXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
  {
    gggXYZ[i][0]=0.0;
    gggXYZ[i][1]=0.0;
    gggXYZ[i][2]=0.0;
  }

  switch(Lmax)
  {
%s
  }

  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    hessXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
  {
    gggXYZ[i][0] *= NormFactor[i];
    gggXYZ[i][1] *= NormFactor[i];
    gggXYZ[i][2] *= NormFactor[i];
  }
}

"""

  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}

    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)

    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    if gx != 0:
      body_str += "    gradXYZ[%d][0] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gradXYZ[%d][1] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gradXYZ[%d][2] = %s;\n"%(idx,gz)

    axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        h_s = diff(diff(gto_s, si), sj)
        hess_val = replace_common_subexpressions(h_s, slist)
        if hess_val != 0:
          body_str += "    hessXYZ[%d](%d,%d) = %s;\n"%(idx,ii,jj,hess_val)

    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        for kk,sk in enumerate(axis_syms):
          ggg_s = diff(diff(diff(gto_s, si), sj), sk)
          ggg_val = replace_common_subexpressions(ggg_s, slist)

          if ggg_val != 0:
            body_str += "    gggXYZ[%d][%d](%d,%d) = %s;\n"%(idx,ii,jj,kk,ggg_val)

  return out_str%body_str


def gen_evaluate_third_deriv_only():
  out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateThirdDerivOnly(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
  {
    gggXYZ[i][0]=0.0;
    gggXYZ[i][1]=0.0;
    gggXYZ[i][2]=0.0;
  }

  switch(Lmax)
  {
%s
  }

  for (int i=0; i<ntot; i++)
  {
    gggXYZ[i][0] *= NormFactor[i];
    gggXYZ[i][1] *= NormFactor[i];
    gggXYZ[i][2] *= NormFactor[i];
  }
}

"""

  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  case_has_content = False
  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      body_str += "  case %d:\n"%lmax
      if not case_has_content and curr_lmax != -1:
        body_str += "        ; // empty statement\n"
      curr_lmax = lmax
      case_has_content = False

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}

    axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        for kk,sk in enumerate(axis_syms):
          ggg_s = diff(diff(diff(gto_s, si), sj), sk)
          ggg_val = replace_common_subexpressions(ggg_s, slist)

          if ggg_val != 0:
            body_str += "    gggXYZ[%d][%d](%d,%d) = %s;\n"%(idx,ii,jj,kk,ggg_val)
            case_has_content = True

  return out_str%body_str

def gen_soa_evaluate_bare():
  out_str = """
template<class T>
void SoaCartesianTensor<T>::evaluate_bare(T x, T y, T z, T* restrict XYZ) const
{
  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;
  const T x4=x3*x, y4=y3*y, z4=z3*z;
  const T x5=x4*x, y5=y4*y, z5=z4*z;
  switch(Lmax)
  {
%s
  }
}
"""
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}
    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s; // %s\n"%(idx,val,s)

  return out_str%body_str


def gen_soa_evaluate_vgl():
  out_str = """
template<class T>
void SoaCartesianTensor<T>::evaluateVGL(T x, T y, T z)
{

  constexpr T czero(0);
  cXYZ=czero;

  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;
  const T x4=x3*x, y4=y3*y, z4=z3*z;
  const T x5=x4*x, y5=y4*y, z5=z4*z;
  T* restrict XYZ=cXYZ.data(0);
  T* restrict gr0=cXYZ.data(1);
  T* restrict gr1=cXYZ.data(2);
  T* restrict gr2=cXYZ.data(3);
  T* restrict lap=cXYZ.data(4);

  switch(Lmax)
  {
%s
  }

  const size_t ntot=NormFactor.size();
  for (size_t i=0; i<ntot; i++)
  {
    XYZ[i]*= NormFactor[i];
    gr0[i]*= NormFactor[i];
    gr1[i]*= NormFactor[i];
    gr2[i]*= NormFactor[i];
    lap[i]*= NormFactor[i];
  }
}
"""
  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  # put index and values in a list so it can be reversed
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    # Compute derivatives symbolically
    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)
    lap = diff(gto_s, x, 2) + diff(gto_s, y, 2) + diff(gto_s, z, 2)

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}
    val = replace_common_subexpressions(gto_s, slist)
    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    lap_val = replace_common_subexpressions(lap, slist)


    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    if gx != 0:
      body_str += "    gr0[%d] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gr1[%d] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gr2[%d] = %s;\n"%(idx,gz)

    if lap_val != 0:
      body_str += "    lap[%d] = %s;\n"%(idx,lap_val)

  return out_str%body_str

def gen_soa_evaluate_vgh():
  out_str = """
template<class T>
void SoaCartesianTensor<T>::evaluateVGH(T x, T y, T z)
{
  constexpr T czero(0);
  cXYZ=czero;

  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;
  const T x4=x3*x, y4=y3*y, z4=z3*z;
  const T x5=x4*x, y5=y4*y, z5=z4*z;

  T* restrict XYZ=cXYZ.data(0);
  T* restrict gr0=cXYZ.data(1);
  T* restrict gr1=cXYZ.data(2);
  T* restrict gr2=cXYZ.data(3);
  T* restrict h00=cXYZ.data(4);
  T* restrict h01=cXYZ.data(5);
  T* restrict h02=cXYZ.data(6);
  T* restrict h11=cXYZ.data(7);
  T* restrict h12=cXYZ.data(8);
  T* restrict h22=cXYZ.data(9);


  switch(Lmax)
  {
%s
  }

  const size_t ntot=cXYZ.size();
  for(size_t i=0; i<ntot; ++i)
  {
    XYZ[i]*= NormFactor[i];
    gr0[i]*= NormFactor[i];
    gr1[i]*= NormFactor[i];
    gr2[i]*= NormFactor[i];
    h00[i]*= NormFactor[i];
    h01[i]*= NormFactor[i];
    h02[i]*= NormFactor[i];
    h11[i]*= NormFactor[i];
    h12[i]*= NormFactor[i];
    h22[i]*= NormFactor[i];
  }

}
"""

  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}

    # Compute derivatives symbolically
    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)

    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    if gx != 0:
      body_str += "    gr0[%d] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gr1[%d] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gr2[%d] = %s;\n"%(idx,gz)

    axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        if ii <= jj:
          # Compute Hessian elements symbolically
          h_s = diff(diff(gto_s, si), sj)
          hess_val = replace_common_subexpressions(h_s, slist)
          if hess_val != 0 :
            body_str += "    h%d%d[%d] = %s;\n"%(ii,jj,idx,hess_val)

  return out_str%body_str

def gen_soa_evaluate_vghgh():
  out_str = """
template<class T>
void SoaCartesianTensor<T>::evaluateVGHGH(T x, T y, T z)
{
  constexpr T czero(0);
  cXYZ=czero;

  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;
  const T x4=x3*x, y4=y3*y, z4=z3*z;
  const T x5=x4*x, y5=y4*y, z5=z4*z;

  T* restrict XYZ   = cXYZ.data(0);
  T* restrict gr0   = cXYZ.data(1);
  T* restrict gr1   = cXYZ.data(2);
  T* restrict gr2   = cXYZ.data(3);
  T* restrict h00   = cXYZ.data(4);
  T* restrict h01   = cXYZ.data(5);
  T* restrict h02   = cXYZ.data(6);
  T* restrict h11   = cXYZ.data(7);
  T* restrict h12   = cXYZ.data(8);
  T* restrict h22   = cXYZ.data(9);
  T* restrict gh000 = cXYZ.data(10);
  T* restrict gh001 = cXYZ.data(11);
  T* restrict gh002 = cXYZ.data(12);
  T* restrict gh011 = cXYZ.data(13);
  T* restrict gh012 = cXYZ.data(14);
  T* restrict gh022 = cXYZ.data(15);
  T* restrict gh111 = cXYZ.data(16);
  T* restrict gh112 = cXYZ.data(17);
  T* restrict gh122 = cXYZ.data(18);
  T* restrict gh222 = cXYZ.data(19);

  switch(Lmax)
  {
%s
  }

  const size_t ntot=cXYZ.size();
  for(size_t i=0; i<ntot; ++i)
  {
    XYZ[i]   *= NormFactor[i];
    gr0[i]   *= NormFactor[i];
    gr1[i]   *= NormFactor[i];
    gr2[i]   *= NormFactor[i];
    h00[i]   *= NormFactor[i];
    h01[i]   *= NormFactor[i];
    h02[i]   *= NormFactor[i];
    h11[i]   *= NormFactor[i];
    h12[i]   *= NormFactor[i];
    h22[i]   *= NormFactor[i];
    gh000[i] *= NormFactor[i];
    gh001[i] *= NormFactor[i];
    gh002[i] *= NormFactor[i];
    gh011[i] *= NormFactor[i];
    gh012[i] *= NormFactor[i];
    gh022[i] *= NormFactor[i];
    gh111[i] *= NormFactor[i];
    gh112[i] *= NormFactor[i];
    gh122[i] *= NormFactor[i];
    gh222[i] *= NormFactor[i];
  }

}
"""

  gto_s = create_gto_symbolic()
  # just the 'angular' part
  gto_s = gto_s.subs(Symbol('alpha'),0)

  ijk_with_lmax = gen_lmax(get_ijk())
  body_str = ''
  curr_lmax = -1
  ijk_l = [(idx,c) for idx,c in enumerate(ijk_with_lmax)]

  x,y,z = symbols('x y z')

  for idx, (i,j,k,s,lmax) in reversed(ijk_l):
    if lmax != curr_lmax:
      curr_lmax = lmax
      body_str += "  case %d:\n"%curr_lmax

    slist = {Symbol('i'):i, Symbol('j'):j, Symbol('k'):k}

    # Compute derivatives symbolically
    dx = diff(gto_s, x)
    dy = diff(gto_s, y)
    dz = diff(gto_s, z)

    val = replace_common_subexpressions(gto_s, slist)
    body_str += "    XYZ[%d] = %s;     // %s\n"%(idx,val,s)

    gx = replace_common_subexpressions(dx, slist)
    gy = replace_common_subexpressions(dy, slist)
    gz = replace_common_subexpressions(dz, slist)
    if gx != 0:
      body_str += "    gr0[%d] = %s;\n"%(idx,gx)
    if gy != 0:
      body_str += "    gr1[%d] = %s;\n"%(idx,gy)
    if gz != 0:
      body_str += "    gr2[%d] = %s;\n"%(idx,gz)

    axis_syms = [Symbol('x'), Symbol('y'), Symbol('z')]
    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
        if ii <= jj:
          # Compute Hessian elements symbolically
          h_s = diff(diff(gto_s, si), sj)
          hess_val = replace_common_subexpressions(h_s, slist)
          if hess_val != 0 :
            body_str += "    h%d%d[%d] = %s;\n"%(ii,jj,idx,hess_val)

    for ii,si in enumerate(axis_syms):
      for jj,sj in enumerate(axis_syms):
         for kk,sk in enumerate(axis_syms):
           if ii <= jj and jj <= kk:
           # Compute Grad Hessian elements symbolically
             ghess_s = diff(diff(diff(gto_s, si), sj),sk)
             ghess_val = replace_common_subexpressions(ghess_s, slist)
             if ghess_val != 0 :
               body_str += "    gh%d%d%d[%d] = %s;\n"%(ii,jj,kk,idx,ghess_val)

  return out_str%body_str

# A simple template replacement engine.
# Template items to be replaced start on a line with '%'.
def run_template(fname_in, fname_out, bodies):
  out = ''
  with open(fname_in, 'r') as f:
    for line in f:
      if line.startswith('%'):
        key = line.strip()[1:]
        if key in bodies:
          line = bodies[key]
        else:
          print 'Error, template item not found, key:',key, ' line = ',line
      out += line

  with open(fname_out, 'w') as f:
    f.write(out)

dire_codegen_text = """
/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/Numerics/codegen/%(script_name)s and %(template_file_name)s

 Edit %(template_file_name)s, rerun %(script_name)s, and copy the generated file here.
*/
"""

def create_cartesian_tensor_h():
  bodies = dict()
  bodies['evaluate'] = gen_evaluate()
  bodies['evaluate_all'] = gen_evaluate_all()
  bodies['evaluate_with_hessian'] = gen_evaluate_with_hessian()
  bodies['evaluate_with_third_deriv'] = gen_evaluate_with_third_deriv()
  bodies['evaluate_third_deriv_only'] = gen_evaluate_third_deriv_only()
  fname_out = 'CartesianTensor.h'
  fname_in= 'CartesianTensor.h.in'

  bodies['dire_codegen_warning'] = dire_codegen_text%({'script_name':'gen_cartesian_tensor.py', 'template_file_name':fname_in})

  run_template(fname_in, fname_out, bodies)

def create_soa_cartesian_tensor_h():
  bodies = dict()
  bodies['evaluate_bare'] = gen_soa_evaluate_bare()
  bodies['evaluate_vgl'] = gen_soa_evaluate_vgl()
  bodies['evaluate_vgh'] = gen_soa_evaluate_vgh()
  bodies['evaluate_vghgh'] = gen_soa_evaluate_vghgh()
  fname_in= 'SoaCartesianTensor.h.in'
  fname_out = 'SoaCartesianTensor.h'

  bodies['dire_codegen_warning'] = dire_codegen_text%({'script_name':'gen_cartesian_tensor.py', 'template_file_name':fname_in})

  run_template(fname_in, fname_out, bodies)


if __name__ == '__main__':
    #print gen_evaluate()
    #print gen_evaluate_all()
    #print gen_evaluate_with_hessian()
    #print gen_evaluate_with_third_deriv()
    #print gen_evaluate_third_deriv_only()

    # Create CartesianTensor.h from CartesianTensor.h.in
    create_cartesian_tensor_h()

    # Create SoaCartesianTensor.h from SoaCartesianTensor.h.in
    create_soa_cartesian_tensor_h()

