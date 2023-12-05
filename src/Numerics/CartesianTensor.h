//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/Numerics/codegen/gen_cartesian_tensor.py and CartesianTensor.h.in

 Edit CartesianTensor.h.in, rerun gen_cartesian_tensor.py, and copy the generated file here.
*/


#ifndef QMCPLUSPLUS_CARTESIAN_TENSOR_H
#define QMCPLUSPLUS_CARTESIAN_TENSOR_H

#include "OhmmsPETE/Tensor.h"
#include "Configuration.h"

/**
 *
 * The template parameters
 * - T, the value_type, e.g. double
 * - Point_t, a vector type to provide xyz coordinate.
 * Point_t must have the operator[] defined, e.g., TinyVector\<double,3\>.
 *
 * Array ordered as [S,X,Y,Z,XX,YY,ZZ,XY,XZ,YZ,...]
 *    (following Gamess order)
 */
template<class T,
         class Point_t,
         class Tensor_t = qmcplusplus::Tensor<T, 3>,
         class GGG_t    = qmcplusplus::TinyVector<Tensor_t, 3>>
class CartesianTensor
{
public:
  using value_type = T;
  using pos_type   = Point_t;
  using hess_type  = Tensor_t;
  using ggg_type   = GGG_t;
  using This_t     = CartesianTensor<T, Point_t, Tensor_t>;

  /** constructor
   * @param l_max maximum angular momentum
   *
   * Evaluate all the constants and prefactors.
  */
  explicit CartesianTensor(const int l_max);

  ///makes a table of \f$ N(a,b,c) x^a y^b z^c \f$ and their gradients up to Lmax.
  void evaluate(const Point_t& p);

  ///makes a table of \f$ N(a,b,c) x^a y^b z^c \f$ and their gradients up to Lmax.
  void evaluateAll(const Point_t& p);

  void evaluateTest(const Point_t& p);

  ///makes a table of \f$ N(a,b,c) x^a y^b z^c \f$ and their gradients and hessians up to Lmax.
  void evaluateWithHessian(const Point_t& p);

  ///makes a table of \f$ N(a,b,c) x^a y^b z^c \f$ and their gradients and hessians and third derivatives up to Lmax.
  void evaluateWithThirdDeriv(const Point_t& p);

  ///makes a table of Third derivatives of \f$ N(a,b,c) x^a y^b z^c \f$
  void evaluateThirdDerivOnly(const Point_t& p);

  inline value_type getYlm(int lm) const { return XYZ[lm]; }

  inline Point_t getGradYlm(int lm) const { return gradXYZ[lm]; }

  inline value_type getLaplYlm(int lm) const { return laplXYZ[lm]; }

  inline Tensor_t getHessYlm(int lm) const { return hessXYZ[lm]; }

  inline GGG_t getGGGYlm(int lm) const { return gggXYZ[lm]; }

  inline int size() const { return XYZ.size(); }

  inline int lmax() const { return Lmax; }

  inline void getABC(int n, int& a, int& b, int& c);

  int DFactorial(int num) { return (num < 2) ? 1 : num * DFactorial(num - 2); }

  ///maximum angular momentum for the center
  int Lmax;

  ///values  XYZ\f$=x^a y^b z^c \f$
  std::vector<value_type> XYZ;
  /// Normalization factors
  std::vector<value_type> NormFactor;

  std::vector<Point_t> gradXYZ;
  // laplacian
  std::vector<value_type> laplXYZ;

  std::vector<hess_type> hessXYZ;

  std::vector<ggg_type> gggXYZ;
};
template<class T, class Point_t, class Tensor_t, class GGG_t>
CartesianTensor<T, Point_t, Tensor_t, GGG_t>::CartesianTensor(const int l_max) : Lmax(l_max)
{
  if (Lmax < 0 || Lmax > 6)
  {
    std::cerr << "CartesianTensor can't handle Lmax > 6 or Lmax < 0.\n";
    APP_ABORT("");
  }
  int ntot = 0;
  for (int i = 0; i <= Lmax; i++)
    ntot += (i + 1) * (i + 2) / 2;
  XYZ.resize(ntot);
  gradXYZ.resize(ntot);
  laplXYZ.resize(ntot);
  hessXYZ.resize(ntot);
  gggXYZ.resize(ntot);
  NormFactor.resize(ntot, 1);
  int p = 0;
  int a = 0, b = 0, c = 0;
  const double pi = 4.0 * atan(1.0);
  for (int l = 0; l <= Lmax; l++)
  {
    int n = (l + 1) * (l + 2) / 2;
    for (int k = 0; k < n; k++)
    {
      getABC(p, a, b, c);
      // factor of (alpha^(l+3/2))^(1/2) goes into the radial function
      // mmorales: HACK HACK HACK, to avoid modifyng the radial functions,
      //           I add a term to the normalization to cancel the term
      //           coming from the Spherical Harmonics
      //           NormL = pow(2,L+1)*sqrt(2.0/static_cast<real_type>(DFactorial(2*l+1)))*pow(2.0/pi,0.25)
      double L = static_cast<double>(l);
      double NormL =
          std::pow(2, L + 1) * sqrt(2.0 / static_cast<double>(DFactorial(2 * l + 1))) * std::pow(2.0 / pi, 0.25);
      NormFactor[p++] = std::pow(2.0 / pi, 0.75) * std::pow(4.0, 0.5 * (a + b + c)) *
          std::sqrt(1.0 /
                    static_cast<double>((DFactorial(2 * a - 1) * DFactorial(2 * b - 1) * DFactorial(2 * c - 1)))) /
          NormL;
    }
  }
}


template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::evaluate(const Point_t& p)
{
  value_type x = p[0], y = p[1], z = p[2];
  value_type x2 = x * x, y2 = y * y, z2 = z * z;
  value_type x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  value_type x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  value_type x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;
  switch (Lmax)
  {
  case 6:
    XYZ[83] = x2 * y2 * z2; // X2Y2Z2
    XYZ[82] = x * y2 * z3;  // Z3Y2X
    XYZ[81] = x2 * y * z3;  // Z3X2Y
    XYZ[80] = x * y3 * z2;  // Y3Z2X
    XYZ[79] = x2 * y3 * z;  // Y3X2Z
    XYZ[78] = x3 * y * z2;  // X3Z2Y
    XYZ[77] = x3 * y2 * z;  // X3Y2Z
    XYZ[76] = y3 * z3;      // Y3Z3
    XYZ[75] = x3 * z3;      // X3Z3
    XYZ[74] = x3 * y3;      // X3Y3
    XYZ[73] = x * y * z4;   // Z4XY
    XYZ[72] = x * y4 * z;   // Y4XZ
    XYZ[71] = x4 * y * z;   // X4YZ
    XYZ[70] = y2 * z4;      // Z4Y2
    XYZ[69] = x2 * z4;      // Z4X2
    XYZ[68] = y4 * z2;      // Y4Z2
    XYZ[67] = x2 * y4;      // Y4X2
    XYZ[66] = x4 * z2;      // X4Z2
    XYZ[65] = x4 * y2;      // X4Y2
    XYZ[64] = y * z * z4;   // Z5Y
    XYZ[63] = x * z * z4;   // Z5X
    XYZ[62] = y * y4 * z;   // Y5Z
    XYZ[61] = x * y * y4;   // Y5X
    XYZ[60] = x * x4 * z;   // X5Z
    XYZ[59] = x * x4 * y;   // X5Y
    XYZ[58] = z * z5;       // Z6
    XYZ[57] = y * y5;       // Y6
    XYZ[56] = x * x5;       // X6
  case 5:
    XYZ[55] = x * y2 * z2; // YYZZX
    XYZ[54] = x2 * y * z2; // XXZZY
    XYZ[53] = x2 * y2 * z; // XXYYZ
    XYZ[52] = x * y * z3;  // ZZZXY
    XYZ[51] = x * y3 * z;  // YYYXZ
    XYZ[50] = x3 * y * z;  // XXXYZ
    XYZ[49] = y2 * z3;     // ZZZYY
    XYZ[48] = x2 * z3;     // ZZZXX
    XYZ[47] = y3 * z2;     // YYYZZ
    XYZ[46] = x2 * y3;     // YYYXX
    XYZ[45] = x3 * z2;     // XXXZZ
    XYZ[44] = x3 * y2;     // XXXYY
    XYZ[43] = y * z4;      // ZZZZY
    XYZ[42] = x * z4;      // ZZZZX
    XYZ[41] = y4 * z;      // YYYYZ
    XYZ[40] = x * y4;      // YYYYX
    XYZ[39] = x4 * z;      // XXXXZ
    XYZ[38] = x4 * y;      // XXXXY
    XYZ[37] = z * z4;      // ZZZZZ
    XYZ[36] = y * y4;      // YYYYY
    XYZ[35] = x * x4;      // XXXXX
  case 4:
    XYZ[34] = x * y * z2; // ZZXY
    XYZ[33] = x * y2 * z; // YYXZ
    XYZ[32] = x2 * y * z; // XXYZ
    XYZ[31] = y2 * z2;    // YYZZ
    XYZ[30] = x2 * z2;    // XXZZ
    XYZ[29] = x2 * y2;    // XXYY
    XYZ[28] = y * z3;     // ZZZY
    XYZ[27] = x * z3;     // ZZZX
    XYZ[26] = y3 * z;     // YYYZ
    XYZ[25] = x * y3;     // YYYX
    XYZ[24] = x3 * z;     // XXXZ
    XYZ[23] = x3 * y;     // XXXY
    XYZ[22] = z4;         // ZZZZ
    XYZ[21] = y4;         // YYYY
    XYZ[20] = x4;         // XXXX
  case 3:
    XYZ[19] = x * y * z; // XYZ
    XYZ[18] = y * z2;    // ZZY
    XYZ[17] = x * z2;    // ZZX
    XYZ[16] = y2 * z;    // YYZ
    XYZ[15] = x * y2;    // YYX
    XYZ[14] = x2 * z;    // XXZ
    XYZ[13] = x2 * y;    // XXY
    XYZ[12] = z3;        // ZZZ
    XYZ[11] = y3;        // YYY
    XYZ[10] = x3;        // XXX
  case 2:
    XYZ[9] = y * z; // YZ
    XYZ[8] = x * z; // XZ
    XYZ[7] = x * y; // XY
    XYZ[6] = z2;    // ZZ
    XYZ[5] = y2;    // YY
    XYZ[4] = x2;    // XX
  case 1:
    XYZ[3] = z; // Z
    XYZ[2] = y; // Y
    XYZ[1] = x; // X
  case 0:
    XYZ[0] = 1; // S
  }
  for (int i = 0; i < XYZ.size(); i++)
    XYZ[i] *= NormFactor[i];
}


template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::evaluateAll(const Point_t& p)
{
  value_type x = p[0], y = p[1], z = p[2];
  value_type x2 = x * x, y2 = y * y, z2 = z * z;
  value_type x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  value_type x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  value_type x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;
  int ntot = XYZ.size();
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] = 0.0;
  for (int i = 0; i < ntot; i++)
    laplXYZ[i] = 0.0;

  switch (Lmax)
  {
  case 6:
    XYZ[83]        = x2 * y2 * z2; // X2Y2Z2
    gradXYZ[83][0] = 2 * x * y2 * z2;
    gradXYZ[83][1] = 2 * x2 * y * z2;
    gradXYZ[83][2] = 2 * x2 * y2 * z;
    laplXYZ[83]    = 2 * x2 * y2 + 2 * x2 * z2 + 2 * y2 * z2;
    XYZ[82]        = x * y2 * z3; // Z3Y2X
    gradXYZ[82][0] = y2 * z3;
    gradXYZ[82][1] = 2 * x * y * z3;
    gradXYZ[82][2] = 3 * x * y2 * z2;
    laplXYZ[82]    = 6 * x * y2 * z + 2 * x * z3;
    XYZ[81]        = x2 * y * z3; // Z3X2Y
    gradXYZ[81][0] = 2 * x * y * z3;
    gradXYZ[81][1] = x2 * z3;
    gradXYZ[81][2] = 3 * x2 * y * z2;
    laplXYZ[81]    = 6 * x2 * y * z + 2 * y * z3;
    XYZ[80]        = x * y3 * z2; // Y3Z2X
    gradXYZ[80][0] = y3 * z2;
    gradXYZ[80][1] = 3 * x * y2 * z2;
    gradXYZ[80][2] = 2 * x * y3 * z;
    laplXYZ[80]    = 6 * x * y * z2 + 2 * x * y3;
    XYZ[79]        = x2 * y3 * z; // Y3X2Z
    gradXYZ[79][0] = 2 * x * y3 * z;
    gradXYZ[79][1] = 3 * x2 * y2 * z;
    gradXYZ[79][2] = x2 * y3;
    laplXYZ[79]    = 6 * x2 * y * z + 2 * y3 * z;
    XYZ[78]        = x3 * y * z2; // X3Z2Y
    gradXYZ[78][0] = 3 * x2 * y * z2;
    gradXYZ[78][1] = x3 * z2;
    gradXYZ[78][2] = 2 * x3 * y * z;
    laplXYZ[78]    = 6 * x * y * z2 + 2 * x3 * y;
    XYZ[77]        = x3 * y2 * z; // X3Y2Z
    gradXYZ[77][0] = 3 * x2 * y2 * z;
    gradXYZ[77][1] = 2 * x3 * y * z;
    gradXYZ[77][2] = x3 * y2;
    laplXYZ[77]    = 6 * x * y2 * z + 2 * x3 * z;
    XYZ[76]        = y3 * z3; // Y3Z3
    gradXYZ[76][1] = 3 * y2 * z3;
    gradXYZ[76][2] = 3 * y3 * z2;
    laplXYZ[76]    = 6 * y * z3 + 6 * y3 * z;
    XYZ[75]        = x3 * z3; // X3Z3
    gradXYZ[75][0] = 3 * x2 * z3;
    gradXYZ[75][2] = 3 * x3 * z2;
    laplXYZ[75]    = 6 * x * z3 + 6 * x3 * z;
    XYZ[74]        = x3 * y3; // X3Y3
    gradXYZ[74][0] = 3 * x2 * y3;
    gradXYZ[74][1] = 3 * x3 * y2;
    laplXYZ[74]    = 6 * x * y3 + 6 * x3 * y;
    XYZ[73]        = x * y * z4; // Z4XY
    gradXYZ[73][0] = y * z4;
    gradXYZ[73][1] = x * z4;
    gradXYZ[73][2] = 4 * x * y * z3;
    laplXYZ[73]    = 12 * x * y * z2;
    XYZ[72]        = x * y4 * z; // Y4XZ
    gradXYZ[72][0] = y4 * z;
    gradXYZ[72][1] = 4 * x * y3 * z;
    gradXYZ[72][2] = x * y4;
    laplXYZ[72]    = 12 * x * y2 * z;
    XYZ[71]        = x4 * y * z; // X4YZ
    gradXYZ[71][0] = 4 * x3 * y * z;
    gradXYZ[71][1] = x4 * z;
    gradXYZ[71][2] = x4 * y;
    laplXYZ[71]    = 12 * x2 * y * z;
    XYZ[70]        = y2 * z4; // Z4Y2
    gradXYZ[70][1] = 2 * y * z4;
    gradXYZ[70][2] = 4 * y2 * z3;
    laplXYZ[70]    = 12 * y2 * z2 + 2 * z4;
    XYZ[69]        = x2 * z4; // Z4X2
    gradXYZ[69][0] = 2 * x * z4;
    gradXYZ[69][2] = 4 * x2 * z3;
    laplXYZ[69]    = 12 * x2 * z2 + 2 * z4;
    XYZ[68]        = y4 * z2; // Y4Z2
    gradXYZ[68][1] = 4 * y3 * z2;
    gradXYZ[68][2] = 2 * y4 * z;
    laplXYZ[68]    = 12 * y2 * z2 + 2 * y4;
    XYZ[67]        = x2 * y4; // Y4X2
    gradXYZ[67][0] = 2 * x * y4;
    gradXYZ[67][1] = 4 * x2 * y3;
    laplXYZ[67]    = 12 * x2 * y2 + 2 * y4;
    XYZ[66]        = x4 * z2; // X4Z2
    gradXYZ[66][0] = 4 * x3 * z2;
    gradXYZ[66][2] = 2 * x4 * z;
    laplXYZ[66]    = 12 * x2 * z2 + 2 * x4;
    XYZ[65]        = x4 * y2; // X4Y2
    gradXYZ[65][0] = 4 * x3 * y2;
    gradXYZ[65][1] = 2 * x4 * y;
    laplXYZ[65]    = 12 * x2 * y2 + 2 * x4;
    XYZ[64]        = y * z * z4; // Z5Y
    gradXYZ[64][1] = z * z4;
    gradXYZ[64][2] = 5 * y * z4;
    laplXYZ[64]    = 20 * y * z3;
    XYZ[63]        = x * z * z4; // Z5X
    gradXYZ[63][0] = z * z4;
    gradXYZ[63][2] = 5 * x * z4;
    laplXYZ[63]    = 20 * x * z3;
    XYZ[62]        = y * y4 * z; // Y5Z
    gradXYZ[62][1] = 5 * y4 * z;
    gradXYZ[62][2] = y * y4;
    laplXYZ[62]    = 20 * y3 * z;
    XYZ[61]        = x * y * y4; // Y5X
    gradXYZ[61][0] = y * y4;
    gradXYZ[61][1] = 5 * x * y4;
    laplXYZ[61]    = 20 * x * y3;
    XYZ[60]        = x * x4 * z; // X5Z
    gradXYZ[60][0] = 5 * x4 * z;
    gradXYZ[60][2] = x * x4;
    laplXYZ[60]    = 20 * x3 * z;
    XYZ[59]        = x * x4 * y; // X5Y
    gradXYZ[59][0] = 5 * x4 * y;
    gradXYZ[59][1] = x * x4;
    laplXYZ[59]    = 20 * x3 * y;
    XYZ[58]        = z * z5; // Z6
    gradXYZ[58][2] = 6 * z * z4;
    laplXYZ[58]    = 30 * z4;
    XYZ[57]        = y * y5; // Y6
    gradXYZ[57][1] = 6 * y * y4;
    laplXYZ[57]    = 30 * y4;
    XYZ[56]        = x * x5; // X6
    gradXYZ[56][0] = 6 * x * x4;
    laplXYZ[56]    = 30 * x4;
  case 5:
    XYZ[55]        = x * y2 * z2; // YYZZX
    gradXYZ[55][0] = y2 * z2;
    gradXYZ[55][1] = 2 * x * y * z2;
    gradXYZ[55][2] = 2 * x * y2 * z;
    laplXYZ[55]    = 2 * x * y2 + 2 * x * z2;
    XYZ[54]        = x2 * y * z2; // XXZZY
    gradXYZ[54][0] = 2 * x * y * z2;
    gradXYZ[54][1] = x2 * z2;
    gradXYZ[54][2] = 2 * x2 * y * z;
    laplXYZ[54]    = 2 * x2 * y + 2 * y * z2;
    XYZ[53]        = x2 * y2 * z; // XXYYZ
    gradXYZ[53][0] = 2 * x * y2 * z;
    gradXYZ[53][1] = 2 * x2 * y * z;
    gradXYZ[53][2] = x2 * y2;
    laplXYZ[53]    = 2 * x2 * z + 2 * y2 * z;
    XYZ[52]        = x * y * z3; // ZZZXY
    gradXYZ[52][0] = y * z3;
    gradXYZ[52][1] = x * z3;
    gradXYZ[52][2] = 3 * x * y * z2;
    laplXYZ[52]    = 6 * x * y * z;
    XYZ[51]        = x * y3 * z; // YYYXZ
    gradXYZ[51][0] = y3 * z;
    gradXYZ[51][1] = 3 * x * y2 * z;
    gradXYZ[51][2] = x * y3;
    laplXYZ[51]    = 6 * x * y * z;
    XYZ[50]        = x3 * y * z; // XXXYZ
    gradXYZ[50][0] = 3 * x2 * y * z;
    gradXYZ[50][1] = x3 * z;
    gradXYZ[50][2] = x3 * y;
    laplXYZ[50]    = 6 * x * y * z;
    XYZ[49]        = y2 * z3; // ZZZYY
    gradXYZ[49][1] = 2 * y * z3;
    gradXYZ[49][2] = 3 * y2 * z2;
    laplXYZ[49]    = 6 * y2 * z + 2 * z3;
    XYZ[48]        = x2 * z3; // ZZZXX
    gradXYZ[48][0] = 2 * x * z3;
    gradXYZ[48][2] = 3 * x2 * z2;
    laplXYZ[48]    = 6 * x2 * z + 2 * z3;
    XYZ[47]        = y3 * z2; // YYYZZ
    gradXYZ[47][1] = 3 * y2 * z2;
    gradXYZ[47][2] = 2 * y3 * z;
    laplXYZ[47]    = 6 * y * z2 + 2 * y3;
    XYZ[46]        = x2 * y3; // YYYXX
    gradXYZ[46][0] = 2 * x * y3;
    gradXYZ[46][1] = 3 * x2 * y2;
    laplXYZ[46]    = 6 * x2 * y + 2 * y3;
    XYZ[45]        = x3 * z2; // XXXZZ
    gradXYZ[45][0] = 3 * x2 * z2;
    gradXYZ[45][2] = 2 * x3 * z;
    laplXYZ[45]    = 6 * x * z2 + 2 * x3;
    XYZ[44]        = x3 * y2; // XXXYY
    gradXYZ[44][0] = 3 * x2 * y2;
    gradXYZ[44][1] = 2 * x3 * y;
    laplXYZ[44]    = 6 * x * y2 + 2 * x3;
    XYZ[43]        = y * z4; // ZZZZY
    gradXYZ[43][1] = z4;
    gradXYZ[43][2] = 4 * y * z3;
    laplXYZ[43]    = 12 * y * z2;
    XYZ[42]        = x * z4; // ZZZZX
    gradXYZ[42][0] = z4;
    gradXYZ[42][2] = 4 * x * z3;
    laplXYZ[42]    = 12 * x * z2;
    XYZ[41]        = y4 * z; // YYYYZ
    gradXYZ[41][1] = 4 * y3 * z;
    gradXYZ[41][2] = y4;
    laplXYZ[41]    = 12 * y2 * z;
    XYZ[40]        = x * y4; // YYYYX
    gradXYZ[40][0] = y4;
    gradXYZ[40][1] = 4 * x * y3;
    laplXYZ[40]    = 12 * x * y2;
    XYZ[39]        = x4 * z; // XXXXZ
    gradXYZ[39][0] = 4 * x3 * z;
    gradXYZ[39][2] = x4;
    laplXYZ[39]    = 12 * x2 * z;
    XYZ[38]        = x4 * y; // XXXXY
    gradXYZ[38][0] = 4 * x3 * y;
    gradXYZ[38][1] = x4;
    laplXYZ[38]    = 12 * x2 * y;
    XYZ[37]        = z * z4; // ZZZZZ
    gradXYZ[37][2] = 5 * z4;
    laplXYZ[37]    = 20 * z3;
    XYZ[36]        = y * y4; // YYYYY
    gradXYZ[36][1] = 5 * y4;
    laplXYZ[36]    = 20 * y3;
    XYZ[35]        = x * x4; // XXXXX
    gradXYZ[35][0] = 5 * x4;
    laplXYZ[35]    = 20 * x3;
  case 4:
    XYZ[34]        = x * y * z2; // ZZXY
    gradXYZ[34][0] = y * z2;
    gradXYZ[34][1] = x * z2;
    gradXYZ[34][2] = 2 * x * y * z;
    laplXYZ[34]    = 2 * x * y;
    XYZ[33]        = x * y2 * z; // YYXZ
    gradXYZ[33][0] = y2 * z;
    gradXYZ[33][1] = 2 * x * y * z;
    gradXYZ[33][2] = x * y2;
    laplXYZ[33]    = 2 * x * z;
    XYZ[32]        = x2 * y * z; // XXYZ
    gradXYZ[32][0] = 2 * x * y * z;
    gradXYZ[32][1] = x2 * z;
    gradXYZ[32][2] = x2 * y;
    laplXYZ[32]    = 2 * y * z;
    XYZ[31]        = y2 * z2; // YYZZ
    gradXYZ[31][1] = 2 * y * z2;
    gradXYZ[31][2] = 2 * y2 * z;
    laplXYZ[31]    = 2 * y2 + 2 * z2;
    XYZ[30]        = x2 * z2; // XXZZ
    gradXYZ[30][0] = 2 * x * z2;
    gradXYZ[30][2] = 2 * x2 * z;
    laplXYZ[30]    = 2 * x2 + 2 * z2;
    XYZ[29]        = x2 * y2; // XXYY
    gradXYZ[29][0] = 2 * x * y2;
    gradXYZ[29][1] = 2 * x2 * y;
    laplXYZ[29]    = 2 * x2 + 2 * y2;
    XYZ[28]        = y * z3; // ZZZY
    gradXYZ[28][1] = z3;
    gradXYZ[28][2] = 3 * y * z2;
    laplXYZ[28]    = 6 * y * z;
    XYZ[27]        = x * z3; // ZZZX
    gradXYZ[27][0] = z3;
    gradXYZ[27][2] = 3 * x * z2;
    laplXYZ[27]    = 6 * x * z;
    XYZ[26]        = y3 * z; // YYYZ
    gradXYZ[26][1] = 3 * y2 * z;
    gradXYZ[26][2] = y3;
    laplXYZ[26]    = 6 * y * z;
    XYZ[25]        = x * y3; // YYYX
    gradXYZ[25][0] = y3;
    gradXYZ[25][1] = 3 * x * y2;
    laplXYZ[25]    = 6 * x * y;
    XYZ[24]        = x3 * z; // XXXZ
    gradXYZ[24][0] = 3 * x2 * z;
    gradXYZ[24][2] = x3;
    laplXYZ[24]    = 6 * x * z;
    XYZ[23]        = x3 * y; // XXXY
    gradXYZ[23][0] = 3 * x2 * y;
    gradXYZ[23][1] = x3;
    laplXYZ[23]    = 6 * x * y;
    XYZ[22]        = z4; // ZZZZ
    gradXYZ[22][2] = 4 * z3;
    laplXYZ[22]    = 12 * z2;
    XYZ[21]        = y4; // YYYY
    gradXYZ[21][1] = 4 * y3;
    laplXYZ[21]    = 12 * y2;
    XYZ[20]        = x4; // XXXX
    gradXYZ[20][0] = 4 * x3;
    laplXYZ[20]    = 12 * x2;
  case 3:
    XYZ[19]        = x * y * z; // XYZ
    gradXYZ[19][0] = y * z;
    gradXYZ[19][1] = x * z;
    gradXYZ[19][2] = x * y;
    XYZ[18]        = y * z2; // ZZY
    gradXYZ[18][1] = z2;
    gradXYZ[18][2] = 2 * y * z;
    laplXYZ[18]    = 2 * y;
    XYZ[17]        = x * z2; // ZZX
    gradXYZ[17][0] = z2;
    gradXYZ[17][2] = 2 * x * z;
    laplXYZ[17]    = 2 * x;
    XYZ[16]        = y2 * z; // YYZ
    gradXYZ[16][1] = 2 * y * z;
    gradXYZ[16][2] = y2;
    laplXYZ[16]    = 2 * z;
    XYZ[15]        = x * y2; // YYX
    gradXYZ[15][0] = y2;
    gradXYZ[15][1] = 2 * x * y;
    laplXYZ[15]    = 2 * x;
    XYZ[14]        = x2 * z; // XXZ
    gradXYZ[14][0] = 2 * x * z;
    gradXYZ[14][2] = x2;
    laplXYZ[14]    = 2 * z;
    XYZ[13]        = x2 * y; // XXY
    gradXYZ[13][0] = 2 * x * y;
    gradXYZ[13][1] = x2;
    laplXYZ[13]    = 2 * y;
    XYZ[12]        = z3; // ZZZ
    gradXYZ[12][2] = 3 * z2;
    laplXYZ[12]    = 6 * z;
    XYZ[11]        = y3; // YYY
    gradXYZ[11][1] = 3 * y2;
    laplXYZ[11]    = 6 * y;
    XYZ[10]        = x3; // XXX
    gradXYZ[10][0] = 3 * x2;
    laplXYZ[10]    = 6 * x;
  case 2:
    XYZ[9]        = y * z; // YZ
    gradXYZ[9][1] = z;
    gradXYZ[9][2] = y;
    XYZ[8]        = x * z; // XZ
    gradXYZ[8][0] = z;
    gradXYZ[8][2] = x;
    XYZ[7]        = x * y; // XY
    gradXYZ[7][0] = y;
    gradXYZ[7][1] = x;
    XYZ[6]        = z2; // ZZ
    gradXYZ[6][2] = 2 * z;
    laplXYZ[6]    = 2;
    XYZ[5]        = y2; // YY
    gradXYZ[5][1] = 2 * y;
    laplXYZ[5]    = 2;
    XYZ[4]        = x2; // XX
    gradXYZ[4][0] = 2 * x;
    laplXYZ[4]    = 2;
  case 1:
    XYZ[3]        = z; // Z
    gradXYZ[3][2] = 1;
    XYZ[2]        = y; // Y
    gradXYZ[2][1] = 1;
    XYZ[1]        = x; // X
    gradXYZ[1][0] = 1;
  case 0:
    XYZ[0] = 1; // S
  }
  for (int i = 0; i < ntot; i++)
    XYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    laplXYZ[i] *= NormFactor[i];
}

#if 0
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
  case 6:
    XYZ[83] = x2*y2*z2;     // X2Y2Z2
    gradXYZ[83][0] = 2*x*y2*z2;
    gradXYZ[83][1] = 2*x2*y*z2;
    gradXYZ[83][2] = 2*x2*y2*z;

  }
  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    laplXYZ[i]*= NormFactor[i];
}
#endif


template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::evaluateWithHessian(const Point_t& p)
{
  value_type x = p[0], y = p[1], z = p[2];
  value_type x2 = x * x, y2 = y * y, z2 = z * z;
  value_type x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  value_type x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  value_type x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;
  int ntot = XYZ.size();
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] = 0.0;
  for (int i = 0; i < ntot; i++)
    hessXYZ[i] = 0.0;

  switch (Lmax)
  {
  case 6:
    XYZ[83]           = x2 * y2 * z2; // X2Y2Z2
    gradXYZ[83][0]    = 2 * x * y2 * z2;
    gradXYZ[83][1]    = 2 * x2 * y * z2;
    gradXYZ[83][2]    = 2 * x2 * y2 * z;
    hessXYZ[83](0, 0) = 2 * y2 * z2;
    hessXYZ[83](0, 1) = 4 * x * y * z2;
    hessXYZ[83](0, 2) = 4 * x * y2 * z;
    hessXYZ[83](1, 0) = 4 * x * y * z2;
    hessXYZ[83](1, 1) = 2 * x2 * z2;
    hessXYZ[83](1, 2) = 4 * x2 * y * z;
    hessXYZ[83](2, 0) = 4 * x * y2 * z;
    hessXYZ[83](2, 1) = 4 * x2 * y * z;
    hessXYZ[83](2, 2) = 2 * x2 * y2;
    XYZ[82]           = x * y2 * z3; // Z3Y2X
    gradXYZ[82][0]    = y2 * z3;
    gradXYZ[82][1]    = 2 * x * y * z3;
    gradXYZ[82][2]    = 3 * x * y2 * z2;
    hessXYZ[82](0, 1) = 2 * y * z3;
    hessXYZ[82](0, 2) = 3 * y2 * z2;
    hessXYZ[82](1, 0) = 2 * y * z3;
    hessXYZ[82](1, 1) = 2 * x * z3;
    hessXYZ[82](1, 2) = 6 * x * y * z2;
    hessXYZ[82](2, 0) = 3 * y2 * z2;
    hessXYZ[82](2, 1) = 6 * x * y * z2;
    hessXYZ[82](2, 2) = 6 * x * y2 * z;
    XYZ[81]           = x2 * y * z3; // Z3X2Y
    gradXYZ[81][0]    = 2 * x * y * z3;
    gradXYZ[81][1]    = x2 * z3;
    gradXYZ[81][2]    = 3 * x2 * y * z2;
    hessXYZ[81](0, 0) = 2 * y * z3;
    hessXYZ[81](0, 1) = 2 * x * z3;
    hessXYZ[81](0, 2) = 6 * x * y * z2;
    hessXYZ[81](1, 0) = 2 * x * z3;
    hessXYZ[81](1, 2) = 3 * x2 * z2;
    hessXYZ[81](2, 0) = 6 * x * y * z2;
    hessXYZ[81](2, 1) = 3 * x2 * z2;
    hessXYZ[81](2, 2) = 6 * x2 * y * z;
    XYZ[80]           = x * y3 * z2; // Y3Z2X
    gradXYZ[80][0]    = y3 * z2;
    gradXYZ[80][1]    = 3 * x * y2 * z2;
    gradXYZ[80][2]    = 2 * x * y3 * z;
    hessXYZ[80](0, 1) = 3 * y2 * z2;
    hessXYZ[80](0, 2) = 2 * y3 * z;
    hessXYZ[80](1, 0) = 3 * y2 * z2;
    hessXYZ[80](1, 1) = 6 * x * y * z2;
    hessXYZ[80](1, 2) = 6 * x * y2 * z;
    hessXYZ[80](2, 0) = 2 * y3 * z;
    hessXYZ[80](2, 1) = 6 * x * y2 * z;
    hessXYZ[80](2, 2) = 2 * x * y3;
    XYZ[79]           = x2 * y3 * z; // Y3X2Z
    gradXYZ[79][0]    = 2 * x * y3 * z;
    gradXYZ[79][1]    = 3 * x2 * y2 * z;
    gradXYZ[79][2]    = x2 * y3;
    hessXYZ[79](0, 0) = 2 * y3 * z;
    hessXYZ[79](0, 1) = 6 * x * y2 * z;
    hessXYZ[79](0, 2) = 2 * x * y3;
    hessXYZ[79](1, 0) = 6 * x * y2 * z;
    hessXYZ[79](1, 1) = 6 * x2 * y * z;
    hessXYZ[79](1, 2) = 3 * x2 * y2;
    hessXYZ[79](2, 0) = 2 * x * y3;
    hessXYZ[79](2, 1) = 3 * x2 * y2;
    XYZ[78]           = x3 * y * z2; // X3Z2Y
    gradXYZ[78][0]    = 3 * x2 * y * z2;
    gradXYZ[78][1]    = x3 * z2;
    gradXYZ[78][2]    = 2 * x3 * y * z;
    hessXYZ[78](0, 0) = 6 * x * y * z2;
    hessXYZ[78](0, 1) = 3 * x2 * z2;
    hessXYZ[78](0, 2) = 6 * x2 * y * z;
    hessXYZ[78](1, 0) = 3 * x2 * z2;
    hessXYZ[78](1, 2) = 2 * x3 * z;
    hessXYZ[78](2, 0) = 6 * x2 * y * z;
    hessXYZ[78](2, 1) = 2 * x3 * z;
    hessXYZ[78](2, 2) = 2 * x3 * y;
    XYZ[77]           = x3 * y2 * z; // X3Y2Z
    gradXYZ[77][0]    = 3 * x2 * y2 * z;
    gradXYZ[77][1]    = 2 * x3 * y * z;
    gradXYZ[77][2]    = x3 * y2;
    hessXYZ[77](0, 0) = 6 * x * y2 * z;
    hessXYZ[77](0, 1) = 6 * x2 * y * z;
    hessXYZ[77](0, 2) = 3 * x2 * y2;
    hessXYZ[77](1, 0) = 6 * x2 * y * z;
    hessXYZ[77](1, 1) = 2 * x3 * z;
    hessXYZ[77](1, 2) = 2 * x3 * y;
    hessXYZ[77](2, 0) = 3 * x2 * y2;
    hessXYZ[77](2, 1) = 2 * x3 * y;
    XYZ[76]           = y3 * z3; // Y3Z3
    gradXYZ[76][1]    = 3 * y2 * z3;
    gradXYZ[76][2]    = 3 * y3 * z2;
    hessXYZ[76](1, 1) = 6 * y * z3;
    hessXYZ[76](1, 2) = 9 * y2 * z2;
    hessXYZ[76](2, 1) = 9 * y2 * z2;
    hessXYZ[76](2, 2) = 6 * y3 * z;
    XYZ[75]           = x3 * z3; // X3Z3
    gradXYZ[75][0]    = 3 * x2 * z3;
    gradXYZ[75][2]    = 3 * x3 * z2;
    hessXYZ[75](0, 0) = 6 * x * z3;
    hessXYZ[75](0, 2) = 9 * x2 * z2;
    hessXYZ[75](2, 0) = 9 * x2 * z2;
    hessXYZ[75](2, 2) = 6 * x3 * z;
    XYZ[74]           = x3 * y3; // X3Y3
    gradXYZ[74][0]    = 3 * x2 * y3;
    gradXYZ[74][1]    = 3 * x3 * y2;
    hessXYZ[74](0, 0) = 6 * x * y3;
    hessXYZ[74](0, 1) = 9 * x2 * y2;
    hessXYZ[74](1, 0) = 9 * x2 * y2;
    hessXYZ[74](1, 1) = 6 * x3 * y;
    XYZ[73]           = x * y * z4; // Z4XY
    gradXYZ[73][0]    = y * z4;
    gradXYZ[73][1]    = x * z4;
    gradXYZ[73][2]    = 4 * x * y * z3;
    hessXYZ[73](0, 1) = z4;
    hessXYZ[73](0, 2) = 4 * y * z3;
    hessXYZ[73](1, 0) = z4;
    hessXYZ[73](1, 2) = 4 * x * z3;
    hessXYZ[73](2, 0) = 4 * y * z3;
    hessXYZ[73](2, 1) = 4 * x * z3;
    hessXYZ[73](2, 2) = 12 * x * y * z2;
    XYZ[72]           = x * y4 * z; // Y4XZ
    gradXYZ[72][0]    = y4 * z;
    gradXYZ[72][1]    = 4 * x * y3 * z;
    gradXYZ[72][2]    = x * y4;
    hessXYZ[72](0, 1) = 4 * y3 * z;
    hessXYZ[72](0, 2) = y4;
    hessXYZ[72](1, 0) = 4 * y3 * z;
    hessXYZ[72](1, 1) = 12 * x * y2 * z;
    hessXYZ[72](1, 2) = 4 * x * y3;
    hessXYZ[72](2, 0) = y4;
    hessXYZ[72](2, 1) = 4 * x * y3;
    XYZ[71]           = x4 * y * z; // X4YZ
    gradXYZ[71][0]    = 4 * x3 * y * z;
    gradXYZ[71][1]    = x4 * z;
    gradXYZ[71][2]    = x4 * y;
    hessXYZ[71](0, 0) = 12 * x2 * y * z;
    hessXYZ[71](0, 1) = 4 * x3 * z;
    hessXYZ[71](0, 2) = 4 * x3 * y;
    hessXYZ[71](1, 0) = 4 * x3 * z;
    hessXYZ[71](1, 2) = x4;
    hessXYZ[71](2, 0) = 4 * x3 * y;
    hessXYZ[71](2, 1) = x4;
    XYZ[70]           = y2 * z4; // Z4Y2
    gradXYZ[70][1]    = 2 * y * z4;
    gradXYZ[70][2]    = 4 * y2 * z3;
    hessXYZ[70](1, 1) = 2 * z4;
    hessXYZ[70](1, 2) = 8 * y * z3;
    hessXYZ[70](2, 1) = 8 * y * z3;
    hessXYZ[70](2, 2) = 12 * y2 * z2;
    XYZ[69]           = x2 * z4; // Z4X2
    gradXYZ[69][0]    = 2 * x * z4;
    gradXYZ[69][2]    = 4 * x2 * z3;
    hessXYZ[69](0, 0) = 2 * z4;
    hessXYZ[69](0, 2) = 8 * x * z3;
    hessXYZ[69](2, 0) = 8 * x * z3;
    hessXYZ[69](2, 2) = 12 * x2 * z2;
    XYZ[68]           = y4 * z2; // Y4Z2
    gradXYZ[68][1]    = 4 * y3 * z2;
    gradXYZ[68][2]    = 2 * y4 * z;
    hessXYZ[68](1, 1) = 12 * y2 * z2;
    hessXYZ[68](1, 2) = 8 * y3 * z;
    hessXYZ[68](2, 1) = 8 * y3 * z;
    hessXYZ[68](2, 2) = 2 * y4;
    XYZ[67]           = x2 * y4; // Y4X2
    gradXYZ[67][0]    = 2 * x * y4;
    gradXYZ[67][1]    = 4 * x2 * y3;
    hessXYZ[67](0, 0) = 2 * y4;
    hessXYZ[67](0, 1) = 8 * x * y3;
    hessXYZ[67](1, 0) = 8 * x * y3;
    hessXYZ[67](1, 1) = 12 * x2 * y2;
    XYZ[66]           = x4 * z2; // X4Z2
    gradXYZ[66][0]    = 4 * x3 * z2;
    gradXYZ[66][2]    = 2 * x4 * z;
    hessXYZ[66](0, 0) = 12 * x2 * z2;
    hessXYZ[66](0, 2) = 8 * x3 * z;
    hessXYZ[66](2, 0) = 8 * x3 * z;
    hessXYZ[66](2, 2) = 2 * x4;
    XYZ[65]           = x4 * y2; // X4Y2
    gradXYZ[65][0]    = 4 * x3 * y2;
    gradXYZ[65][1]    = 2 * x4 * y;
    hessXYZ[65](0, 0) = 12 * x2 * y2;
    hessXYZ[65](0, 1) = 8 * x3 * y;
    hessXYZ[65](1, 0) = 8 * x3 * y;
    hessXYZ[65](1, 1) = 2 * x4;
    XYZ[64]           = y * z * z4; // Z5Y
    gradXYZ[64][1]    = z * z4;
    gradXYZ[64][2]    = 5 * y * z4;
    hessXYZ[64](1, 2) = 5 * z4;
    hessXYZ[64](2, 1) = 5 * z4;
    hessXYZ[64](2, 2) = 20 * y * z3;
    XYZ[63]           = x * z * z4; // Z5X
    gradXYZ[63][0]    = z * z4;
    gradXYZ[63][2]    = 5 * x * z4;
    hessXYZ[63](0, 2) = 5 * z4;
    hessXYZ[63](2, 0) = 5 * z4;
    hessXYZ[63](2, 2) = 20 * x * z3;
    XYZ[62]           = y * y4 * z; // Y5Z
    gradXYZ[62][1]    = 5 * y4 * z;
    gradXYZ[62][2]    = y * y4;
    hessXYZ[62](1, 1) = 20 * y3 * z;
    hessXYZ[62](1, 2) = 5 * y4;
    hessXYZ[62](2, 1) = 5 * y4;
    XYZ[61]           = x * y * y4; // Y5X
    gradXYZ[61][0]    = y * y4;
    gradXYZ[61][1]    = 5 * x * y4;
    hessXYZ[61](0, 1) = 5 * y4;
    hessXYZ[61](1, 0) = 5 * y4;
    hessXYZ[61](1, 1) = 20 * x * y3;
    XYZ[60]           = x * x4 * z; // X5Z
    gradXYZ[60][0]    = 5 * x4 * z;
    gradXYZ[60][2]    = x * x4;
    hessXYZ[60](0, 0) = 20 * x3 * z;
    hessXYZ[60](0, 2) = 5 * x4;
    hessXYZ[60](2, 0) = 5 * x4;
    XYZ[59]           = x * x4 * y; // X5Y
    gradXYZ[59][0]    = 5 * x4 * y;
    gradXYZ[59][1]    = x * x4;
    hessXYZ[59](0, 0) = 20 * x3 * y;
    hessXYZ[59](0, 1) = 5 * x4;
    hessXYZ[59](1, 0) = 5 * x4;
    XYZ[58]           = z * z5; // Z6
    gradXYZ[58][2]    = 6 * z * z4;
    hessXYZ[58](2, 2) = 30 * z4;
    XYZ[57]           = y * y5; // Y6
    gradXYZ[57][1]    = 6 * y * y4;
    hessXYZ[57](1, 1) = 30 * y4;
    XYZ[56]           = x * x5; // X6
    gradXYZ[56][0]    = 6 * x * x4;
    hessXYZ[56](0, 0) = 30 * x4;
  case 5:
    XYZ[55]           = x * y2 * z2; // YYZZX
    gradXYZ[55][0]    = y2 * z2;
    gradXYZ[55][1]    = 2 * x * y * z2;
    gradXYZ[55][2]    = 2 * x * y2 * z;
    hessXYZ[55](0, 1) = 2 * y * z2;
    hessXYZ[55](0, 2) = 2 * y2 * z;
    hessXYZ[55](1, 0) = 2 * y * z2;
    hessXYZ[55](1, 1) = 2 * x * z2;
    hessXYZ[55](1, 2) = 4 * x * y * z;
    hessXYZ[55](2, 0) = 2 * y2 * z;
    hessXYZ[55](2, 1) = 4 * x * y * z;
    hessXYZ[55](2, 2) = 2 * x * y2;
    XYZ[54]           = x2 * y * z2; // XXZZY
    gradXYZ[54][0]    = 2 * x * y * z2;
    gradXYZ[54][1]    = x2 * z2;
    gradXYZ[54][2]    = 2 * x2 * y * z;
    hessXYZ[54](0, 0) = 2 * y * z2;
    hessXYZ[54](0, 1) = 2 * x * z2;
    hessXYZ[54](0, 2) = 4 * x * y * z;
    hessXYZ[54](1, 0) = 2 * x * z2;
    hessXYZ[54](1, 2) = 2 * x2 * z;
    hessXYZ[54](2, 0) = 4 * x * y * z;
    hessXYZ[54](2, 1) = 2 * x2 * z;
    hessXYZ[54](2, 2) = 2 * x2 * y;
    XYZ[53]           = x2 * y2 * z; // XXYYZ
    gradXYZ[53][0]    = 2 * x * y2 * z;
    gradXYZ[53][1]    = 2 * x2 * y * z;
    gradXYZ[53][2]    = x2 * y2;
    hessXYZ[53](0, 0) = 2 * y2 * z;
    hessXYZ[53](0, 1) = 4 * x * y * z;
    hessXYZ[53](0, 2) = 2 * x * y2;
    hessXYZ[53](1, 0) = 4 * x * y * z;
    hessXYZ[53](1, 1) = 2 * x2 * z;
    hessXYZ[53](1, 2) = 2 * x2 * y;
    hessXYZ[53](2, 0) = 2 * x * y2;
    hessXYZ[53](2, 1) = 2 * x2 * y;
    XYZ[52]           = x * y * z3; // ZZZXY
    gradXYZ[52][0]    = y * z3;
    gradXYZ[52][1]    = x * z3;
    gradXYZ[52][2]    = 3 * x * y * z2;
    hessXYZ[52](0, 1) = z3;
    hessXYZ[52](0, 2) = 3 * y * z2;
    hessXYZ[52](1, 0) = z3;
    hessXYZ[52](1, 2) = 3 * x * z2;
    hessXYZ[52](2, 0) = 3 * y * z2;
    hessXYZ[52](2, 1) = 3 * x * z2;
    hessXYZ[52](2, 2) = 6 * x * y * z;
    XYZ[51]           = x * y3 * z; // YYYXZ
    gradXYZ[51][0]    = y3 * z;
    gradXYZ[51][1]    = 3 * x * y2 * z;
    gradXYZ[51][2]    = x * y3;
    hessXYZ[51](0, 1) = 3 * y2 * z;
    hessXYZ[51](0, 2) = y3;
    hessXYZ[51](1, 0) = 3 * y2 * z;
    hessXYZ[51](1, 1) = 6 * x * y * z;
    hessXYZ[51](1, 2) = 3 * x * y2;
    hessXYZ[51](2, 0) = y3;
    hessXYZ[51](2, 1) = 3 * x * y2;
    XYZ[50]           = x3 * y * z; // XXXYZ
    gradXYZ[50][0]    = 3 * x2 * y * z;
    gradXYZ[50][1]    = x3 * z;
    gradXYZ[50][2]    = x3 * y;
    hessXYZ[50](0, 0) = 6 * x * y * z;
    hessXYZ[50](0, 1) = 3 * x2 * z;
    hessXYZ[50](0, 2) = 3 * x2 * y;
    hessXYZ[50](1, 0) = 3 * x2 * z;
    hessXYZ[50](1, 2) = x3;
    hessXYZ[50](2, 0) = 3 * x2 * y;
    hessXYZ[50](2, 1) = x3;
    XYZ[49]           = y2 * z3; // ZZZYY
    gradXYZ[49][1]    = 2 * y * z3;
    gradXYZ[49][2]    = 3 * y2 * z2;
    hessXYZ[49](1, 1) = 2 * z3;
    hessXYZ[49](1, 2) = 6 * y * z2;
    hessXYZ[49](2, 1) = 6 * y * z2;
    hessXYZ[49](2, 2) = 6 * y2 * z;
    XYZ[48]           = x2 * z3; // ZZZXX
    gradXYZ[48][0]    = 2 * x * z3;
    gradXYZ[48][2]    = 3 * x2 * z2;
    hessXYZ[48](0, 0) = 2 * z3;
    hessXYZ[48](0, 2) = 6 * x * z2;
    hessXYZ[48](2, 0) = 6 * x * z2;
    hessXYZ[48](2, 2) = 6 * x2 * z;
    XYZ[47]           = y3 * z2; // YYYZZ
    gradXYZ[47][1]    = 3 * y2 * z2;
    gradXYZ[47][2]    = 2 * y3 * z;
    hessXYZ[47](1, 1) = 6 * y * z2;
    hessXYZ[47](1, 2) = 6 * y2 * z;
    hessXYZ[47](2, 1) = 6 * y2 * z;
    hessXYZ[47](2, 2) = 2 * y3;
    XYZ[46]           = x2 * y3; // YYYXX
    gradXYZ[46][0]    = 2 * x * y3;
    gradXYZ[46][1]    = 3 * x2 * y2;
    hessXYZ[46](0, 0) = 2 * y3;
    hessXYZ[46](0, 1) = 6 * x * y2;
    hessXYZ[46](1, 0) = 6 * x * y2;
    hessXYZ[46](1, 1) = 6 * x2 * y;
    XYZ[45]           = x3 * z2; // XXXZZ
    gradXYZ[45][0]    = 3 * x2 * z2;
    gradXYZ[45][2]    = 2 * x3 * z;
    hessXYZ[45](0, 0) = 6 * x * z2;
    hessXYZ[45](0, 2) = 6 * x2 * z;
    hessXYZ[45](2, 0) = 6 * x2 * z;
    hessXYZ[45](2, 2) = 2 * x3;
    XYZ[44]           = x3 * y2; // XXXYY
    gradXYZ[44][0]    = 3 * x2 * y2;
    gradXYZ[44][1]    = 2 * x3 * y;
    hessXYZ[44](0, 0) = 6 * x * y2;
    hessXYZ[44](0, 1) = 6 * x2 * y;
    hessXYZ[44](1, 0) = 6 * x2 * y;
    hessXYZ[44](1, 1) = 2 * x3;
    XYZ[43]           = y * z4; // ZZZZY
    gradXYZ[43][1]    = z4;
    gradXYZ[43][2]    = 4 * y * z3;
    hessXYZ[43](1, 2) = 4 * z3;
    hessXYZ[43](2, 1) = 4 * z3;
    hessXYZ[43](2, 2) = 12 * y * z2;
    XYZ[42]           = x * z4; // ZZZZX
    gradXYZ[42][0]    = z4;
    gradXYZ[42][2]    = 4 * x * z3;
    hessXYZ[42](0, 2) = 4 * z3;
    hessXYZ[42](2, 0) = 4 * z3;
    hessXYZ[42](2, 2) = 12 * x * z2;
    XYZ[41]           = y4 * z; // YYYYZ
    gradXYZ[41][1]    = 4 * y3 * z;
    gradXYZ[41][2]    = y4;
    hessXYZ[41](1, 1) = 12 * y2 * z;
    hessXYZ[41](1, 2) = 4 * y3;
    hessXYZ[41](2, 1) = 4 * y3;
    XYZ[40]           = x * y4; // YYYYX
    gradXYZ[40][0]    = y4;
    gradXYZ[40][1]    = 4 * x * y3;
    hessXYZ[40](0, 1) = 4 * y3;
    hessXYZ[40](1, 0) = 4 * y3;
    hessXYZ[40](1, 1) = 12 * x * y2;
    XYZ[39]           = x4 * z; // XXXXZ
    gradXYZ[39][0]    = 4 * x3 * z;
    gradXYZ[39][2]    = x4;
    hessXYZ[39](0, 0) = 12 * x2 * z;
    hessXYZ[39](0, 2) = 4 * x3;
    hessXYZ[39](2, 0) = 4 * x3;
    XYZ[38]           = x4 * y; // XXXXY
    gradXYZ[38][0]    = 4 * x3 * y;
    gradXYZ[38][1]    = x4;
    hessXYZ[38](0, 0) = 12 * x2 * y;
    hessXYZ[38](0, 1) = 4 * x3;
    hessXYZ[38](1, 0) = 4 * x3;
    XYZ[37]           = z * z4; // ZZZZZ
    gradXYZ[37][2]    = 5 * z4;
    hessXYZ[37](2, 2) = 20 * z3;
    XYZ[36]           = y * y4; // YYYYY
    gradXYZ[36][1]    = 5 * y4;
    hessXYZ[36](1, 1) = 20 * y3;
    XYZ[35]           = x * x4; // XXXXX
    gradXYZ[35][0]    = 5 * x4;
    hessXYZ[35](0, 0) = 20 * x3;
  case 4:
    XYZ[34]           = x * y * z2; // ZZXY
    gradXYZ[34][0]    = y * z2;
    gradXYZ[34][1]    = x * z2;
    gradXYZ[34][2]    = 2 * x * y * z;
    hessXYZ[34](0, 1) = z2;
    hessXYZ[34](0, 2) = 2 * y * z;
    hessXYZ[34](1, 0) = z2;
    hessXYZ[34](1, 2) = 2 * x * z;
    hessXYZ[34](2, 0) = 2 * y * z;
    hessXYZ[34](2, 1) = 2 * x * z;
    hessXYZ[34](2, 2) = 2 * x * y;
    XYZ[33]           = x * y2 * z; // YYXZ
    gradXYZ[33][0]    = y2 * z;
    gradXYZ[33][1]    = 2 * x * y * z;
    gradXYZ[33][2]    = x * y2;
    hessXYZ[33](0, 1) = 2 * y * z;
    hessXYZ[33](0, 2) = y2;
    hessXYZ[33](1, 0) = 2 * y * z;
    hessXYZ[33](1, 1) = 2 * x * z;
    hessXYZ[33](1, 2) = 2 * x * y;
    hessXYZ[33](2, 0) = y2;
    hessXYZ[33](2, 1) = 2 * x * y;
    XYZ[32]           = x2 * y * z; // XXYZ
    gradXYZ[32][0]    = 2 * x * y * z;
    gradXYZ[32][1]    = x2 * z;
    gradXYZ[32][2]    = x2 * y;
    hessXYZ[32](0, 0) = 2 * y * z;
    hessXYZ[32](0, 1) = 2 * x * z;
    hessXYZ[32](0, 2) = 2 * x * y;
    hessXYZ[32](1, 0) = 2 * x * z;
    hessXYZ[32](1, 2) = x2;
    hessXYZ[32](2, 0) = 2 * x * y;
    hessXYZ[32](2, 1) = x2;
    XYZ[31]           = y2 * z2; // YYZZ
    gradXYZ[31][1]    = 2 * y * z2;
    gradXYZ[31][2]    = 2 * y2 * z;
    hessXYZ[31](1, 1) = 2 * z2;
    hessXYZ[31](1, 2) = 4 * y * z;
    hessXYZ[31](2, 1) = 4 * y * z;
    hessXYZ[31](2, 2) = 2 * y2;
    XYZ[30]           = x2 * z2; // XXZZ
    gradXYZ[30][0]    = 2 * x * z2;
    gradXYZ[30][2]    = 2 * x2 * z;
    hessXYZ[30](0, 0) = 2 * z2;
    hessXYZ[30](0, 2) = 4 * x * z;
    hessXYZ[30](2, 0) = 4 * x * z;
    hessXYZ[30](2, 2) = 2 * x2;
    XYZ[29]           = x2 * y2; // XXYY
    gradXYZ[29][0]    = 2 * x * y2;
    gradXYZ[29][1]    = 2 * x2 * y;
    hessXYZ[29](0, 0) = 2 * y2;
    hessXYZ[29](0, 1) = 4 * x * y;
    hessXYZ[29](1, 0) = 4 * x * y;
    hessXYZ[29](1, 1) = 2 * x2;
    XYZ[28]           = y * z3; // ZZZY
    gradXYZ[28][1]    = z3;
    gradXYZ[28][2]    = 3 * y * z2;
    hessXYZ[28](1, 2) = 3 * z2;
    hessXYZ[28](2, 1) = 3 * z2;
    hessXYZ[28](2, 2) = 6 * y * z;
    XYZ[27]           = x * z3; // ZZZX
    gradXYZ[27][0]    = z3;
    gradXYZ[27][2]    = 3 * x * z2;
    hessXYZ[27](0, 2) = 3 * z2;
    hessXYZ[27](2, 0) = 3 * z2;
    hessXYZ[27](2, 2) = 6 * x * z;
    XYZ[26]           = y3 * z; // YYYZ
    gradXYZ[26][1]    = 3 * y2 * z;
    gradXYZ[26][2]    = y3;
    hessXYZ[26](1, 1) = 6 * y * z;
    hessXYZ[26](1, 2) = 3 * y2;
    hessXYZ[26](2, 1) = 3 * y2;
    XYZ[25]           = x * y3; // YYYX
    gradXYZ[25][0]    = y3;
    gradXYZ[25][1]    = 3 * x * y2;
    hessXYZ[25](0, 1) = 3 * y2;
    hessXYZ[25](1, 0) = 3 * y2;
    hessXYZ[25](1, 1) = 6 * x * y;
    XYZ[24]           = x3 * z; // XXXZ
    gradXYZ[24][0]    = 3 * x2 * z;
    gradXYZ[24][2]    = x3;
    hessXYZ[24](0, 0) = 6 * x * z;
    hessXYZ[24](0, 2) = 3 * x2;
    hessXYZ[24](2, 0) = 3 * x2;
    XYZ[23]           = x3 * y; // XXXY
    gradXYZ[23][0]    = 3 * x2 * y;
    gradXYZ[23][1]    = x3;
    hessXYZ[23](0, 0) = 6 * x * y;
    hessXYZ[23](0, 1) = 3 * x2;
    hessXYZ[23](1, 0) = 3 * x2;
    XYZ[22]           = z4; // ZZZZ
    gradXYZ[22][2]    = 4 * z3;
    hessXYZ[22](2, 2) = 12 * z2;
    XYZ[21]           = y4; // YYYY
    gradXYZ[21][1]    = 4 * y3;
    hessXYZ[21](1, 1) = 12 * y2;
    XYZ[20]           = x4; // XXXX
    gradXYZ[20][0]    = 4 * x3;
    hessXYZ[20](0, 0) = 12 * x2;
  case 3:
    XYZ[19]           = x * y * z; // XYZ
    gradXYZ[19][0]    = y * z;
    gradXYZ[19][1]    = x * z;
    gradXYZ[19][2]    = x * y;
    hessXYZ[19](0, 1) = z;
    hessXYZ[19](0, 2) = y;
    hessXYZ[19](1, 0) = z;
    hessXYZ[19](1, 2) = x;
    hessXYZ[19](2, 0) = y;
    hessXYZ[19](2, 1) = x;
    XYZ[18]           = y * z2; // ZZY
    gradXYZ[18][1]    = z2;
    gradXYZ[18][2]    = 2 * y * z;
    hessXYZ[18](1, 2) = 2 * z;
    hessXYZ[18](2, 1) = 2 * z;
    hessXYZ[18](2, 2) = 2 * y;
    XYZ[17]           = x * z2; // ZZX
    gradXYZ[17][0]    = z2;
    gradXYZ[17][2]    = 2 * x * z;
    hessXYZ[17](0, 2) = 2 * z;
    hessXYZ[17](2, 0) = 2 * z;
    hessXYZ[17](2, 2) = 2 * x;
    XYZ[16]           = y2 * z; // YYZ
    gradXYZ[16][1]    = 2 * y * z;
    gradXYZ[16][2]    = y2;
    hessXYZ[16](1, 1) = 2 * z;
    hessXYZ[16](1, 2) = 2 * y;
    hessXYZ[16](2, 1) = 2 * y;
    XYZ[15]           = x * y2; // YYX
    gradXYZ[15][0]    = y2;
    gradXYZ[15][1]    = 2 * x * y;
    hessXYZ[15](0, 1) = 2 * y;
    hessXYZ[15](1, 0) = 2 * y;
    hessXYZ[15](1, 1) = 2 * x;
    XYZ[14]           = x2 * z; // XXZ
    gradXYZ[14][0]    = 2 * x * z;
    gradXYZ[14][2]    = x2;
    hessXYZ[14](0, 0) = 2 * z;
    hessXYZ[14](0, 2) = 2 * x;
    hessXYZ[14](2, 0) = 2 * x;
    XYZ[13]           = x2 * y; // XXY
    gradXYZ[13][0]    = 2 * x * y;
    gradXYZ[13][1]    = x2;
    hessXYZ[13](0, 0) = 2 * y;
    hessXYZ[13](0, 1) = 2 * x;
    hessXYZ[13](1, 0) = 2 * x;
    XYZ[12]           = z3; // ZZZ
    gradXYZ[12][2]    = 3 * z2;
    hessXYZ[12](2, 2) = 6 * z;
    XYZ[11]           = y3; // YYY
    gradXYZ[11][1]    = 3 * y2;
    hessXYZ[11](1, 1) = 6 * y;
    XYZ[10]           = x3; // XXX
    gradXYZ[10][0]    = 3 * x2;
    hessXYZ[10](0, 0) = 6 * x;
  case 2:
    XYZ[9]           = y * z; // YZ
    gradXYZ[9][1]    = z;
    gradXYZ[9][2]    = y;
    hessXYZ[9](1, 2) = 1;
    hessXYZ[9](2, 1) = 1;
    XYZ[8]           = x * z; // XZ
    gradXYZ[8][0]    = z;
    gradXYZ[8][2]    = x;
    hessXYZ[8](0, 2) = 1;
    hessXYZ[8](2, 0) = 1;
    XYZ[7]           = x * y; // XY
    gradXYZ[7][0]    = y;
    gradXYZ[7][1]    = x;
    hessXYZ[7](0, 1) = 1;
    hessXYZ[7](1, 0) = 1;
    XYZ[6]           = z2; // ZZ
    gradXYZ[6][2]    = 2 * z;
    hessXYZ[6](2, 2) = 2;
    XYZ[5]           = y2; // YY
    gradXYZ[5][1]    = 2 * y;
    hessXYZ[5](1, 1) = 2;
    XYZ[4]           = x2; // XX
    gradXYZ[4][0]    = 2 * x;
    hessXYZ[4](0, 0) = 2;
  case 1:
    XYZ[3]        = z; // Z
    gradXYZ[3][2] = 1;
    XYZ[2]        = y; // Y
    gradXYZ[2][1] = 1;
    XYZ[1]        = x; // X
    gradXYZ[1][0] = 1;
  case 0:
    XYZ[0] = 1; // S
  }
  for (int i = 0; i < ntot; i++)
    XYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    hessXYZ[i] *= NormFactor[i];
}


template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::evaluateWithThirdDeriv(const Point_t& p)
{
  value_type x = p[0], y = p[1], z = p[2];
  value_type x2 = x * x, y2 = y * y, z2 = z * z;
  value_type x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  value_type x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  value_type x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;

  int ntot = XYZ.size();
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] = 0.0;
  for (int i = 0; i < ntot; i++)
    hessXYZ[i] = 0.0;
  for (int i = 0; i < ntot; i++)
  {
    gggXYZ[i][0] = 0.0;
    gggXYZ[i][1] = 0.0;
    gggXYZ[i][2] = 0.0;
  }

  switch (Lmax)
  {
  case 6:
    XYZ[83]             = x2 * y2 * z2; // X2Y2Z2
    gradXYZ[83][0]      = 2 * x * y2 * z2;
    gradXYZ[83][1]      = 2 * x2 * y * z2;
    gradXYZ[83][2]      = 2 * x2 * y2 * z;
    hessXYZ[83](0, 0)   = 2 * y2 * z2;
    hessXYZ[83](0, 1)   = 4 * x * y * z2;
    hessXYZ[83](0, 2)   = 4 * x * y2 * z;
    hessXYZ[83](1, 0)   = 4 * x * y * z2;
    hessXYZ[83](1, 1)   = 2 * x2 * z2;
    hessXYZ[83](1, 2)   = 4 * x2 * y * z;
    hessXYZ[83](2, 0)   = 4 * x * y2 * z;
    hessXYZ[83](2, 1)   = 4 * x2 * y * z;
    hessXYZ[83](2, 2)   = 2 * x2 * y2;
    gggXYZ[83][0](0, 1) = 4 * y * z2;
    gggXYZ[83][0](0, 2) = 4 * y2 * z;
    gggXYZ[83][0](1, 0) = 4 * y * z2;
    gggXYZ[83][0](1, 1) = 4 * x * z2;
    gggXYZ[83][0](1, 2) = 8 * x * y * z;
    gggXYZ[83][0](2, 0) = 4 * y2 * z;
    gggXYZ[83][0](2, 1) = 8 * x * y * z;
    gggXYZ[83][0](2, 2) = 4 * x * y2;
    gggXYZ[83][1](0, 0) = 4 * y * z2;
    gggXYZ[83][1](0, 1) = 4 * x * z2;
    gggXYZ[83][1](0, 2) = 8 * x * y * z;
    gggXYZ[83][1](1, 0) = 4 * x * z2;
    gggXYZ[83][1](1, 2) = 4 * x2 * z;
    gggXYZ[83][1](2, 0) = 8 * x * y * z;
    gggXYZ[83][1](2, 1) = 4 * x2 * z;
    gggXYZ[83][1](2, 2) = 4 * x2 * y;
    gggXYZ[83][2](0, 0) = 4 * y2 * z;
    gggXYZ[83][2](0, 1) = 8 * x * y * z;
    gggXYZ[83][2](0, 2) = 4 * x * y2;
    gggXYZ[83][2](1, 0) = 8 * x * y * z;
    gggXYZ[83][2](1, 1) = 4 * x2 * z;
    gggXYZ[83][2](1, 2) = 4 * x2 * y;
    gggXYZ[83][2](2, 0) = 4 * x * y2;
    gggXYZ[83][2](2, 1) = 4 * x2 * y;
    XYZ[82]             = x * y2 * z3; // Z3Y2X
    gradXYZ[82][0]      = y2 * z3;
    gradXYZ[82][1]      = 2 * x * y * z3;
    gradXYZ[82][2]      = 3 * x * y2 * z2;
    hessXYZ[82](0, 1)   = 2 * y * z3;
    hessXYZ[82](0, 2)   = 3 * y2 * z2;
    hessXYZ[82](1, 0)   = 2 * y * z3;
    hessXYZ[82](1, 1)   = 2 * x * z3;
    hessXYZ[82](1, 2)   = 6 * x * y * z2;
    hessXYZ[82](2, 0)   = 3 * y2 * z2;
    hessXYZ[82](2, 1)   = 6 * x * y * z2;
    hessXYZ[82](2, 2)   = 6 * x * y2 * z;
    gggXYZ[82][0](1, 1) = 2 * z3;
    gggXYZ[82][0](1, 2) = 6 * y * z2;
    gggXYZ[82][0](2, 1) = 6 * y * z2;
    gggXYZ[82][0](2, 2) = 6 * y2 * z;
    gggXYZ[82][1](0, 1) = 2 * z3;
    gggXYZ[82][1](0, 2) = 6 * y * z2;
    gggXYZ[82][1](1, 0) = 2 * z3;
    gggXYZ[82][1](1, 2) = 6 * x * z2;
    gggXYZ[82][1](2, 0) = 6 * y * z2;
    gggXYZ[82][1](2, 1) = 6 * x * z2;
    gggXYZ[82][1](2, 2) = 12 * x * y * z;
    gggXYZ[82][2](0, 1) = 6 * y * z2;
    gggXYZ[82][2](0, 2) = 6 * y2 * z;
    gggXYZ[82][2](1, 0) = 6 * y * z2;
    gggXYZ[82][2](1, 1) = 6 * x * z2;
    gggXYZ[82][2](1, 2) = 12 * x * y * z;
    gggXYZ[82][2](2, 0) = 6 * y2 * z;
    gggXYZ[82][2](2, 1) = 12 * x * y * z;
    gggXYZ[82][2](2, 2) = 6 * x * y2;
    XYZ[81]             = x2 * y * z3; // Z3X2Y
    gradXYZ[81][0]      = 2 * x * y * z3;
    gradXYZ[81][1]      = x2 * z3;
    gradXYZ[81][2]      = 3 * x2 * y * z2;
    hessXYZ[81](0, 0)   = 2 * y * z3;
    hessXYZ[81](0, 1)   = 2 * x * z3;
    hessXYZ[81](0, 2)   = 6 * x * y * z2;
    hessXYZ[81](1, 0)   = 2 * x * z3;
    hessXYZ[81](1, 2)   = 3 * x2 * z2;
    hessXYZ[81](2, 0)   = 6 * x * y * z2;
    hessXYZ[81](2, 1)   = 3 * x2 * z2;
    hessXYZ[81](2, 2)   = 6 * x2 * y * z;
    gggXYZ[81][0](0, 1) = 2 * z3;
    gggXYZ[81][0](0, 2) = 6 * y * z2;
    gggXYZ[81][0](1, 0) = 2 * z3;
    gggXYZ[81][0](1, 2) = 6 * x * z2;
    gggXYZ[81][0](2, 0) = 6 * y * z2;
    gggXYZ[81][0](2, 1) = 6 * x * z2;
    gggXYZ[81][0](2, 2) = 12 * x * y * z;
    gggXYZ[81][1](0, 0) = 2 * z3;
    gggXYZ[81][1](0, 2) = 6 * x * z2;
    gggXYZ[81][1](2, 0) = 6 * x * z2;
    gggXYZ[81][1](2, 2) = 6 * x2 * z;
    gggXYZ[81][2](0, 0) = 6 * y * z2;
    gggXYZ[81][2](0, 1) = 6 * x * z2;
    gggXYZ[81][2](0, 2) = 12 * x * y * z;
    gggXYZ[81][2](1, 0) = 6 * x * z2;
    gggXYZ[81][2](1, 2) = 6 * x2 * z;
    gggXYZ[81][2](2, 0) = 12 * x * y * z;
    gggXYZ[81][2](2, 1) = 6 * x2 * z;
    gggXYZ[81][2](2, 2) = 6 * x2 * y;
    XYZ[80]             = x * y3 * z2; // Y3Z2X
    gradXYZ[80][0]      = y3 * z2;
    gradXYZ[80][1]      = 3 * x * y2 * z2;
    gradXYZ[80][2]      = 2 * x * y3 * z;
    hessXYZ[80](0, 1)   = 3 * y2 * z2;
    hessXYZ[80](0, 2)   = 2 * y3 * z;
    hessXYZ[80](1, 0)   = 3 * y2 * z2;
    hessXYZ[80](1, 1)   = 6 * x * y * z2;
    hessXYZ[80](1, 2)   = 6 * x * y2 * z;
    hessXYZ[80](2, 0)   = 2 * y3 * z;
    hessXYZ[80](2, 1)   = 6 * x * y2 * z;
    hessXYZ[80](2, 2)   = 2 * x * y3;
    gggXYZ[80][0](1, 1) = 6 * y * z2;
    gggXYZ[80][0](1, 2) = 6 * y2 * z;
    gggXYZ[80][0](2, 1) = 6 * y2 * z;
    gggXYZ[80][0](2, 2) = 2 * y3;
    gggXYZ[80][1](0, 1) = 6 * y * z2;
    gggXYZ[80][1](0, 2) = 6 * y2 * z;
    gggXYZ[80][1](1, 0) = 6 * y * z2;
    gggXYZ[80][1](1, 1) = 6 * x * z2;
    gggXYZ[80][1](1, 2) = 12 * x * y * z;
    gggXYZ[80][1](2, 0) = 6 * y2 * z;
    gggXYZ[80][1](2, 1) = 12 * x * y * z;
    gggXYZ[80][1](2, 2) = 6 * x * y2;
    gggXYZ[80][2](0, 1) = 6 * y2 * z;
    gggXYZ[80][2](0, 2) = 2 * y3;
    gggXYZ[80][2](1, 0) = 6 * y2 * z;
    gggXYZ[80][2](1, 1) = 12 * x * y * z;
    gggXYZ[80][2](1, 2) = 6 * x * y2;
    gggXYZ[80][2](2, 0) = 2 * y3;
    gggXYZ[80][2](2, 1) = 6 * x * y2;
    XYZ[79]             = x2 * y3 * z; // Y3X2Z
    gradXYZ[79][0]      = 2 * x * y3 * z;
    gradXYZ[79][1]      = 3 * x2 * y2 * z;
    gradXYZ[79][2]      = x2 * y3;
    hessXYZ[79](0, 0)   = 2 * y3 * z;
    hessXYZ[79](0, 1)   = 6 * x * y2 * z;
    hessXYZ[79](0, 2)   = 2 * x * y3;
    hessXYZ[79](1, 0)   = 6 * x * y2 * z;
    hessXYZ[79](1, 1)   = 6 * x2 * y * z;
    hessXYZ[79](1, 2)   = 3 * x2 * y2;
    hessXYZ[79](2, 0)   = 2 * x * y3;
    hessXYZ[79](2, 1)   = 3 * x2 * y2;
    gggXYZ[79][0](0, 1) = 6 * y2 * z;
    gggXYZ[79][0](0, 2) = 2 * y3;
    gggXYZ[79][0](1, 0) = 6 * y2 * z;
    gggXYZ[79][0](1, 1) = 12 * x * y * z;
    gggXYZ[79][0](1, 2) = 6 * x * y2;
    gggXYZ[79][0](2, 0) = 2 * y3;
    gggXYZ[79][0](2, 1) = 6 * x * y2;
    gggXYZ[79][1](0, 0) = 6 * y2 * z;
    gggXYZ[79][1](0, 1) = 12 * x * y * z;
    gggXYZ[79][1](0, 2) = 6 * x * y2;
    gggXYZ[79][1](1, 0) = 12 * x * y * z;
    gggXYZ[79][1](1, 1) = 6 * x2 * z;
    gggXYZ[79][1](1, 2) = 6 * x2 * y;
    gggXYZ[79][1](2, 0) = 6 * x * y2;
    gggXYZ[79][1](2, 1) = 6 * x2 * y;
    gggXYZ[79][2](0, 0) = 2 * y3;
    gggXYZ[79][2](0, 1) = 6 * x * y2;
    gggXYZ[79][2](1, 0) = 6 * x * y2;
    gggXYZ[79][2](1, 1) = 6 * x2 * y;
    XYZ[78]             = x3 * y * z2; // X3Z2Y
    gradXYZ[78][0]      = 3 * x2 * y * z2;
    gradXYZ[78][1]      = x3 * z2;
    gradXYZ[78][2]      = 2 * x3 * y * z;
    hessXYZ[78](0, 0)   = 6 * x * y * z2;
    hessXYZ[78](0, 1)   = 3 * x2 * z2;
    hessXYZ[78](0, 2)   = 6 * x2 * y * z;
    hessXYZ[78](1, 0)   = 3 * x2 * z2;
    hessXYZ[78](1, 2)   = 2 * x3 * z;
    hessXYZ[78](2, 0)   = 6 * x2 * y * z;
    hessXYZ[78](2, 1)   = 2 * x3 * z;
    hessXYZ[78](2, 2)   = 2 * x3 * y;
    gggXYZ[78][0](0, 0) = 6 * y * z2;
    gggXYZ[78][0](0, 1) = 6 * x * z2;
    gggXYZ[78][0](0, 2) = 12 * x * y * z;
    gggXYZ[78][0](1, 0) = 6 * x * z2;
    gggXYZ[78][0](1, 2) = 6 * x2 * z;
    gggXYZ[78][0](2, 0) = 12 * x * y * z;
    gggXYZ[78][0](2, 1) = 6 * x2 * z;
    gggXYZ[78][0](2, 2) = 6 * x2 * y;
    gggXYZ[78][1](0, 0) = 6 * x * z2;
    gggXYZ[78][1](0, 2) = 6 * x2 * z;
    gggXYZ[78][1](2, 0) = 6 * x2 * z;
    gggXYZ[78][1](2, 2) = 2 * x3;
    gggXYZ[78][2](0, 0) = 12 * x * y * z;
    gggXYZ[78][2](0, 1) = 6 * x2 * z;
    gggXYZ[78][2](0, 2) = 6 * x2 * y;
    gggXYZ[78][2](1, 0) = 6 * x2 * z;
    gggXYZ[78][2](1, 2) = 2 * x3;
    gggXYZ[78][2](2, 0) = 6 * x2 * y;
    gggXYZ[78][2](2, 1) = 2 * x3;
    XYZ[77]             = x3 * y2 * z; // X3Y2Z
    gradXYZ[77][0]      = 3 * x2 * y2 * z;
    gradXYZ[77][1]      = 2 * x3 * y * z;
    gradXYZ[77][2]      = x3 * y2;
    hessXYZ[77](0, 0)   = 6 * x * y2 * z;
    hessXYZ[77](0, 1)   = 6 * x2 * y * z;
    hessXYZ[77](0, 2)   = 3 * x2 * y2;
    hessXYZ[77](1, 0)   = 6 * x2 * y * z;
    hessXYZ[77](1, 1)   = 2 * x3 * z;
    hessXYZ[77](1, 2)   = 2 * x3 * y;
    hessXYZ[77](2, 0)   = 3 * x2 * y2;
    hessXYZ[77](2, 1)   = 2 * x3 * y;
    gggXYZ[77][0](0, 0) = 6 * y2 * z;
    gggXYZ[77][0](0, 1) = 12 * x * y * z;
    gggXYZ[77][0](0, 2) = 6 * x * y2;
    gggXYZ[77][0](1, 0) = 12 * x * y * z;
    gggXYZ[77][0](1, 1) = 6 * x2 * z;
    gggXYZ[77][0](1, 2) = 6 * x2 * y;
    gggXYZ[77][0](2, 0) = 6 * x * y2;
    gggXYZ[77][0](2, 1) = 6 * x2 * y;
    gggXYZ[77][1](0, 0) = 12 * x * y * z;
    gggXYZ[77][1](0, 1) = 6 * x2 * z;
    gggXYZ[77][1](0, 2) = 6 * x2 * y;
    gggXYZ[77][1](1, 0) = 6 * x2 * z;
    gggXYZ[77][1](1, 2) = 2 * x3;
    gggXYZ[77][1](2, 0) = 6 * x2 * y;
    gggXYZ[77][1](2, 1) = 2 * x3;
    gggXYZ[77][2](0, 0) = 6 * x * y2;
    gggXYZ[77][2](0, 1) = 6 * x2 * y;
    gggXYZ[77][2](1, 0) = 6 * x2 * y;
    gggXYZ[77][2](1, 1) = 2 * x3;
    XYZ[76]             = y3 * z3; // Y3Z3
    gradXYZ[76][1]      = 3 * y2 * z3;
    gradXYZ[76][2]      = 3 * y3 * z2;
    hessXYZ[76](1, 1)   = 6 * y * z3;
    hessXYZ[76](1, 2)   = 9 * y2 * z2;
    hessXYZ[76](2, 1)   = 9 * y2 * z2;
    hessXYZ[76](2, 2)   = 6 * y3 * z;
    gggXYZ[76][1](1, 1) = 6 * z3;
    gggXYZ[76][1](1, 2) = 18 * y * z2;
    gggXYZ[76][1](2, 1) = 18 * y * z2;
    gggXYZ[76][1](2, 2) = 18 * y2 * z;
    gggXYZ[76][2](1, 1) = 18 * y * z2;
    gggXYZ[76][2](1, 2) = 18 * y2 * z;
    gggXYZ[76][2](2, 1) = 18 * y2 * z;
    gggXYZ[76][2](2, 2) = 6 * y3;
    XYZ[75]             = x3 * z3; // X3Z3
    gradXYZ[75][0]      = 3 * x2 * z3;
    gradXYZ[75][2]      = 3 * x3 * z2;
    hessXYZ[75](0, 0)   = 6 * x * z3;
    hessXYZ[75](0, 2)   = 9 * x2 * z2;
    hessXYZ[75](2, 0)   = 9 * x2 * z2;
    hessXYZ[75](2, 2)   = 6 * x3 * z;
    gggXYZ[75][0](0, 0) = 6 * z3;
    gggXYZ[75][0](0, 2) = 18 * x * z2;
    gggXYZ[75][0](2, 0) = 18 * x * z2;
    gggXYZ[75][0](2, 2) = 18 * x2 * z;
    gggXYZ[75][2](0, 0) = 18 * x * z2;
    gggXYZ[75][2](0, 2) = 18 * x2 * z;
    gggXYZ[75][2](2, 0) = 18 * x2 * z;
    gggXYZ[75][2](2, 2) = 6 * x3;
    XYZ[74]             = x3 * y3; // X3Y3
    gradXYZ[74][0]      = 3 * x2 * y3;
    gradXYZ[74][1]      = 3 * x3 * y2;
    hessXYZ[74](0, 0)   = 6 * x * y3;
    hessXYZ[74](0, 1)   = 9 * x2 * y2;
    hessXYZ[74](1, 0)   = 9 * x2 * y2;
    hessXYZ[74](1, 1)   = 6 * x3 * y;
    gggXYZ[74][0](0, 0) = 6 * y3;
    gggXYZ[74][0](0, 1) = 18 * x * y2;
    gggXYZ[74][0](1, 0) = 18 * x * y2;
    gggXYZ[74][0](1, 1) = 18 * x2 * y;
    gggXYZ[74][1](0, 0) = 18 * x * y2;
    gggXYZ[74][1](0, 1) = 18 * x2 * y;
    gggXYZ[74][1](1, 0) = 18 * x2 * y;
    gggXYZ[74][1](1, 1) = 6 * x3;
    XYZ[73]             = x * y * z4; // Z4XY
    gradXYZ[73][0]      = y * z4;
    gradXYZ[73][1]      = x * z4;
    gradXYZ[73][2]      = 4 * x * y * z3;
    hessXYZ[73](0, 1)   = z4;
    hessXYZ[73](0, 2)   = 4 * y * z3;
    hessXYZ[73](1, 0)   = z4;
    hessXYZ[73](1, 2)   = 4 * x * z3;
    hessXYZ[73](2, 0)   = 4 * y * z3;
    hessXYZ[73](2, 1)   = 4 * x * z3;
    hessXYZ[73](2, 2)   = 12 * x * y * z2;
    gggXYZ[73][0](1, 2) = 4 * z3;
    gggXYZ[73][0](2, 1) = 4 * z3;
    gggXYZ[73][0](2, 2) = 12 * y * z2;
    gggXYZ[73][1](0, 2) = 4 * z3;
    gggXYZ[73][1](2, 0) = 4 * z3;
    gggXYZ[73][1](2, 2) = 12 * x * z2;
    gggXYZ[73][2](0, 1) = 4 * z3;
    gggXYZ[73][2](0, 2) = 12 * y * z2;
    gggXYZ[73][2](1, 0) = 4 * z3;
    gggXYZ[73][2](1, 2) = 12 * x * z2;
    gggXYZ[73][2](2, 0) = 12 * y * z2;
    gggXYZ[73][2](2, 1) = 12 * x * z2;
    gggXYZ[73][2](2, 2) = 24 * x * y * z;
    XYZ[72]             = x * y4 * z; // Y4XZ
    gradXYZ[72][0]      = y4 * z;
    gradXYZ[72][1]      = 4 * x * y3 * z;
    gradXYZ[72][2]      = x * y4;
    hessXYZ[72](0, 1)   = 4 * y3 * z;
    hessXYZ[72](0, 2)   = y4;
    hessXYZ[72](1, 0)   = 4 * y3 * z;
    hessXYZ[72](1, 1)   = 12 * x * y2 * z;
    hessXYZ[72](1, 2)   = 4 * x * y3;
    hessXYZ[72](2, 0)   = y4;
    hessXYZ[72](2, 1)   = 4 * x * y3;
    gggXYZ[72][0](1, 1) = 12 * y2 * z;
    gggXYZ[72][0](1, 2) = 4 * y3;
    gggXYZ[72][0](2, 1) = 4 * y3;
    gggXYZ[72][1](0, 1) = 12 * y2 * z;
    gggXYZ[72][1](0, 2) = 4 * y3;
    gggXYZ[72][1](1, 0) = 12 * y2 * z;
    gggXYZ[72][1](1, 1) = 24 * x * y * z;
    gggXYZ[72][1](1, 2) = 12 * x * y2;
    gggXYZ[72][1](2, 0) = 4 * y3;
    gggXYZ[72][1](2, 1) = 12 * x * y2;
    gggXYZ[72][2](0, 1) = 4 * y3;
    gggXYZ[72][2](1, 0) = 4 * y3;
    gggXYZ[72][2](1, 1) = 12 * x * y2;
    XYZ[71]             = x4 * y * z; // X4YZ
    gradXYZ[71][0]      = 4 * x3 * y * z;
    gradXYZ[71][1]      = x4 * z;
    gradXYZ[71][2]      = x4 * y;
    hessXYZ[71](0, 0)   = 12 * x2 * y * z;
    hessXYZ[71](0, 1)   = 4 * x3 * z;
    hessXYZ[71](0, 2)   = 4 * x3 * y;
    hessXYZ[71](1, 0)   = 4 * x3 * z;
    hessXYZ[71](1, 2)   = x4;
    hessXYZ[71](2, 0)   = 4 * x3 * y;
    hessXYZ[71](2, 1)   = x4;
    gggXYZ[71][0](0, 0) = 24 * x * y * z;
    gggXYZ[71][0](0, 1) = 12 * x2 * z;
    gggXYZ[71][0](0, 2) = 12 * x2 * y;
    gggXYZ[71][0](1, 0) = 12 * x2 * z;
    gggXYZ[71][0](1, 2) = 4 * x3;
    gggXYZ[71][0](2, 0) = 12 * x2 * y;
    gggXYZ[71][0](2, 1) = 4 * x3;
    gggXYZ[71][1](0, 0) = 12 * x2 * z;
    gggXYZ[71][1](0, 2) = 4 * x3;
    gggXYZ[71][1](2, 0) = 4 * x3;
    gggXYZ[71][2](0, 0) = 12 * x2 * y;
    gggXYZ[71][2](0, 1) = 4 * x3;
    gggXYZ[71][2](1, 0) = 4 * x3;
    XYZ[70]             = y2 * z4; // Z4Y2
    gradXYZ[70][1]      = 2 * y * z4;
    gradXYZ[70][2]      = 4 * y2 * z3;
    hessXYZ[70](1, 1)   = 2 * z4;
    hessXYZ[70](1, 2)   = 8 * y * z3;
    hessXYZ[70](2, 1)   = 8 * y * z3;
    hessXYZ[70](2, 2)   = 12 * y2 * z2;
    gggXYZ[70][1](1, 2) = 8 * z3;
    gggXYZ[70][1](2, 1) = 8 * z3;
    gggXYZ[70][1](2, 2) = 24 * y * z2;
    gggXYZ[70][2](1, 1) = 8 * z3;
    gggXYZ[70][2](1, 2) = 24 * y * z2;
    gggXYZ[70][2](2, 1) = 24 * y * z2;
    gggXYZ[70][2](2, 2) = 24 * y2 * z;
    XYZ[69]             = x2 * z4; // Z4X2
    gradXYZ[69][0]      = 2 * x * z4;
    gradXYZ[69][2]      = 4 * x2 * z3;
    hessXYZ[69](0, 0)   = 2 * z4;
    hessXYZ[69](0, 2)   = 8 * x * z3;
    hessXYZ[69](2, 0)   = 8 * x * z3;
    hessXYZ[69](2, 2)   = 12 * x2 * z2;
    gggXYZ[69][0](0, 2) = 8 * z3;
    gggXYZ[69][0](2, 0) = 8 * z3;
    gggXYZ[69][0](2, 2) = 24 * x * z2;
    gggXYZ[69][2](0, 0) = 8 * z3;
    gggXYZ[69][2](0, 2) = 24 * x * z2;
    gggXYZ[69][2](2, 0) = 24 * x * z2;
    gggXYZ[69][2](2, 2) = 24 * x2 * z;
    XYZ[68]             = y4 * z2; // Y4Z2
    gradXYZ[68][1]      = 4 * y3 * z2;
    gradXYZ[68][2]      = 2 * y4 * z;
    hessXYZ[68](1, 1)   = 12 * y2 * z2;
    hessXYZ[68](1, 2)   = 8 * y3 * z;
    hessXYZ[68](2, 1)   = 8 * y3 * z;
    hessXYZ[68](2, 2)   = 2 * y4;
    gggXYZ[68][1](1, 1) = 24 * y * z2;
    gggXYZ[68][1](1, 2) = 24 * y2 * z;
    gggXYZ[68][1](2, 1) = 24 * y2 * z;
    gggXYZ[68][1](2, 2) = 8 * y3;
    gggXYZ[68][2](1, 1) = 24 * y2 * z;
    gggXYZ[68][2](1, 2) = 8 * y3;
    gggXYZ[68][2](2, 1) = 8 * y3;
    XYZ[67]             = x2 * y4; // Y4X2
    gradXYZ[67][0]      = 2 * x * y4;
    gradXYZ[67][1]      = 4 * x2 * y3;
    hessXYZ[67](0, 0)   = 2 * y4;
    hessXYZ[67](0, 1)   = 8 * x * y3;
    hessXYZ[67](1, 0)   = 8 * x * y3;
    hessXYZ[67](1, 1)   = 12 * x2 * y2;
    gggXYZ[67][0](0, 1) = 8 * y3;
    gggXYZ[67][0](1, 0) = 8 * y3;
    gggXYZ[67][0](1, 1) = 24 * x * y2;
    gggXYZ[67][1](0, 0) = 8 * y3;
    gggXYZ[67][1](0, 1) = 24 * x * y2;
    gggXYZ[67][1](1, 0) = 24 * x * y2;
    gggXYZ[67][1](1, 1) = 24 * x2 * y;
    XYZ[66]             = x4 * z2; // X4Z2
    gradXYZ[66][0]      = 4 * x3 * z2;
    gradXYZ[66][2]      = 2 * x4 * z;
    hessXYZ[66](0, 0)   = 12 * x2 * z2;
    hessXYZ[66](0, 2)   = 8 * x3 * z;
    hessXYZ[66](2, 0)   = 8 * x3 * z;
    hessXYZ[66](2, 2)   = 2 * x4;
    gggXYZ[66][0](0, 0) = 24 * x * z2;
    gggXYZ[66][0](0, 2) = 24 * x2 * z;
    gggXYZ[66][0](2, 0) = 24 * x2 * z;
    gggXYZ[66][0](2, 2) = 8 * x3;
    gggXYZ[66][2](0, 0) = 24 * x2 * z;
    gggXYZ[66][2](0, 2) = 8 * x3;
    gggXYZ[66][2](2, 0) = 8 * x3;
    XYZ[65]             = x4 * y2; // X4Y2
    gradXYZ[65][0]      = 4 * x3 * y2;
    gradXYZ[65][1]      = 2 * x4 * y;
    hessXYZ[65](0, 0)   = 12 * x2 * y2;
    hessXYZ[65](0, 1)   = 8 * x3 * y;
    hessXYZ[65](1, 0)   = 8 * x3 * y;
    hessXYZ[65](1, 1)   = 2 * x4;
    gggXYZ[65][0](0, 0) = 24 * x * y2;
    gggXYZ[65][0](0, 1) = 24 * x2 * y;
    gggXYZ[65][0](1, 0) = 24 * x2 * y;
    gggXYZ[65][0](1, 1) = 8 * x3;
    gggXYZ[65][1](0, 0) = 24 * x2 * y;
    gggXYZ[65][1](0, 1) = 8 * x3;
    gggXYZ[65][1](1, 0) = 8 * x3;
    XYZ[64]             = y * z * z4; // Z5Y
    gradXYZ[64][1]      = z * z4;
    gradXYZ[64][2]      = 5 * y * z4;
    hessXYZ[64](1, 2)   = 5 * z4;
    hessXYZ[64](2, 1)   = 5 * z4;
    hessXYZ[64](2, 2)   = 20 * y * z3;
    gggXYZ[64][1](2, 2) = 20 * z3;
    gggXYZ[64][2](1, 2) = 20 * z3;
    gggXYZ[64][2](2, 1) = 20 * z3;
    gggXYZ[64][2](2, 2) = 60 * y * z2;
    XYZ[63]             = x * z * z4; // Z5X
    gradXYZ[63][0]      = z * z4;
    gradXYZ[63][2]      = 5 * x * z4;
    hessXYZ[63](0, 2)   = 5 * z4;
    hessXYZ[63](2, 0)   = 5 * z4;
    hessXYZ[63](2, 2)   = 20 * x * z3;
    gggXYZ[63][0](2, 2) = 20 * z3;
    gggXYZ[63][2](0, 2) = 20 * z3;
    gggXYZ[63][2](2, 0) = 20 * z3;
    gggXYZ[63][2](2, 2) = 60 * x * z2;
    XYZ[62]             = y * y4 * z; // Y5Z
    gradXYZ[62][1]      = 5 * y4 * z;
    gradXYZ[62][2]      = y * y4;
    hessXYZ[62](1, 1)   = 20 * y3 * z;
    hessXYZ[62](1, 2)   = 5 * y4;
    hessXYZ[62](2, 1)   = 5 * y4;
    gggXYZ[62][1](1, 1) = 60 * y2 * z;
    gggXYZ[62][1](1, 2) = 20 * y3;
    gggXYZ[62][1](2, 1) = 20 * y3;
    gggXYZ[62][2](1, 1) = 20 * y3;
    XYZ[61]             = x * y * y4; // Y5X
    gradXYZ[61][0]      = y * y4;
    gradXYZ[61][1]      = 5 * x * y4;
    hessXYZ[61](0, 1)   = 5 * y4;
    hessXYZ[61](1, 0)   = 5 * y4;
    hessXYZ[61](1, 1)   = 20 * x * y3;
    gggXYZ[61][0](1, 1) = 20 * y3;
    gggXYZ[61][1](0, 1) = 20 * y3;
    gggXYZ[61][1](1, 0) = 20 * y3;
    gggXYZ[61][1](1, 1) = 60 * x * y2;
    XYZ[60]             = x * x4 * z; // X5Z
    gradXYZ[60][0]      = 5 * x4 * z;
    gradXYZ[60][2]      = x * x4;
    hessXYZ[60](0, 0)   = 20 * x3 * z;
    hessXYZ[60](0, 2)   = 5 * x4;
    hessXYZ[60](2, 0)   = 5 * x4;
    gggXYZ[60][0](0, 0) = 60 * x2 * z;
    gggXYZ[60][0](0, 2) = 20 * x3;
    gggXYZ[60][0](2, 0) = 20 * x3;
    gggXYZ[60][2](0, 0) = 20 * x3;
    XYZ[59]             = x * x4 * y; // X5Y
    gradXYZ[59][0]      = 5 * x4 * y;
    gradXYZ[59][1]      = x * x4;
    hessXYZ[59](0, 0)   = 20 * x3 * y;
    hessXYZ[59](0, 1)   = 5 * x4;
    hessXYZ[59](1, 0)   = 5 * x4;
    gggXYZ[59][0](0, 0) = 60 * x2 * y;
    gggXYZ[59][0](0, 1) = 20 * x3;
    gggXYZ[59][0](1, 0) = 20 * x3;
    gggXYZ[59][1](0, 0) = 20 * x3;
    XYZ[58]             = z * z5; // Z6
    gradXYZ[58][2]      = 6 * z * z4;
    hessXYZ[58](2, 2)   = 30 * z4;
    gggXYZ[58][2](2, 2) = 120 * z3;
    XYZ[57]             = y * y5; // Y6
    gradXYZ[57][1]      = 6 * y * y4;
    hessXYZ[57](1, 1)   = 30 * y4;
    gggXYZ[57][1](1, 1) = 120 * y3;
    XYZ[56]             = x * x5; // X6
    gradXYZ[56][0]      = 6 * x * x4;
    hessXYZ[56](0, 0)   = 30 * x4;
    gggXYZ[56][0](0, 0) = 120 * x3;
  case 5:
    XYZ[55]             = x * y2 * z2; // YYZZX
    gradXYZ[55][0]      = y2 * z2;
    gradXYZ[55][1]      = 2 * x * y * z2;
    gradXYZ[55][2]      = 2 * x * y2 * z;
    hessXYZ[55](0, 1)   = 2 * y * z2;
    hessXYZ[55](0, 2)   = 2 * y2 * z;
    hessXYZ[55](1, 0)   = 2 * y * z2;
    hessXYZ[55](1, 1)   = 2 * x * z2;
    hessXYZ[55](1, 2)   = 4 * x * y * z;
    hessXYZ[55](2, 0)   = 2 * y2 * z;
    hessXYZ[55](2, 1)   = 4 * x * y * z;
    hessXYZ[55](2, 2)   = 2 * x * y2;
    gggXYZ[55][0](1, 1) = 2 * z2;
    gggXYZ[55][0](1, 2) = 4 * y * z;
    gggXYZ[55][0](2, 1) = 4 * y * z;
    gggXYZ[55][0](2, 2) = 2 * y2;
    gggXYZ[55][1](0, 1) = 2 * z2;
    gggXYZ[55][1](0, 2) = 4 * y * z;
    gggXYZ[55][1](1, 0) = 2 * z2;
    gggXYZ[55][1](1, 2) = 4 * x * z;
    gggXYZ[55][1](2, 0) = 4 * y * z;
    gggXYZ[55][1](2, 1) = 4 * x * z;
    gggXYZ[55][1](2, 2) = 4 * x * y;
    gggXYZ[55][2](0, 1) = 4 * y * z;
    gggXYZ[55][2](0, 2) = 2 * y2;
    gggXYZ[55][2](1, 0) = 4 * y * z;
    gggXYZ[55][2](1, 1) = 4 * x * z;
    gggXYZ[55][2](1, 2) = 4 * x * y;
    gggXYZ[55][2](2, 0) = 2 * y2;
    gggXYZ[55][2](2, 1) = 4 * x * y;
    XYZ[54]             = x2 * y * z2; // XXZZY
    gradXYZ[54][0]      = 2 * x * y * z2;
    gradXYZ[54][1]      = x2 * z2;
    gradXYZ[54][2]      = 2 * x2 * y * z;
    hessXYZ[54](0, 0)   = 2 * y * z2;
    hessXYZ[54](0, 1)   = 2 * x * z2;
    hessXYZ[54](0, 2)   = 4 * x * y * z;
    hessXYZ[54](1, 0)   = 2 * x * z2;
    hessXYZ[54](1, 2)   = 2 * x2 * z;
    hessXYZ[54](2, 0)   = 4 * x * y * z;
    hessXYZ[54](2, 1)   = 2 * x2 * z;
    hessXYZ[54](2, 2)   = 2 * x2 * y;
    gggXYZ[54][0](0, 1) = 2 * z2;
    gggXYZ[54][0](0, 2) = 4 * y * z;
    gggXYZ[54][0](1, 0) = 2 * z2;
    gggXYZ[54][0](1, 2) = 4 * x * z;
    gggXYZ[54][0](2, 0) = 4 * y * z;
    gggXYZ[54][0](2, 1) = 4 * x * z;
    gggXYZ[54][0](2, 2) = 4 * x * y;
    gggXYZ[54][1](0, 0) = 2 * z2;
    gggXYZ[54][1](0, 2) = 4 * x * z;
    gggXYZ[54][1](2, 0) = 4 * x * z;
    gggXYZ[54][1](2, 2) = 2 * x2;
    gggXYZ[54][2](0, 0) = 4 * y * z;
    gggXYZ[54][2](0, 1) = 4 * x * z;
    gggXYZ[54][2](0, 2) = 4 * x * y;
    gggXYZ[54][2](1, 0) = 4 * x * z;
    gggXYZ[54][2](1, 2) = 2 * x2;
    gggXYZ[54][2](2, 0) = 4 * x * y;
    gggXYZ[54][2](2, 1) = 2 * x2;
    XYZ[53]             = x2 * y2 * z; // XXYYZ
    gradXYZ[53][0]      = 2 * x * y2 * z;
    gradXYZ[53][1]      = 2 * x2 * y * z;
    gradXYZ[53][2]      = x2 * y2;
    hessXYZ[53](0, 0)   = 2 * y2 * z;
    hessXYZ[53](0, 1)   = 4 * x * y * z;
    hessXYZ[53](0, 2)   = 2 * x * y2;
    hessXYZ[53](1, 0)   = 4 * x * y * z;
    hessXYZ[53](1, 1)   = 2 * x2 * z;
    hessXYZ[53](1, 2)   = 2 * x2 * y;
    hessXYZ[53](2, 0)   = 2 * x * y2;
    hessXYZ[53](2, 1)   = 2 * x2 * y;
    gggXYZ[53][0](0, 1) = 4 * y * z;
    gggXYZ[53][0](0, 2) = 2 * y2;
    gggXYZ[53][0](1, 0) = 4 * y * z;
    gggXYZ[53][0](1, 1) = 4 * x * z;
    gggXYZ[53][0](1, 2) = 4 * x * y;
    gggXYZ[53][0](2, 0) = 2 * y2;
    gggXYZ[53][0](2, 1) = 4 * x * y;
    gggXYZ[53][1](0, 0) = 4 * y * z;
    gggXYZ[53][1](0, 1) = 4 * x * z;
    gggXYZ[53][1](0, 2) = 4 * x * y;
    gggXYZ[53][1](1, 0) = 4 * x * z;
    gggXYZ[53][1](1, 2) = 2 * x2;
    gggXYZ[53][1](2, 0) = 4 * x * y;
    gggXYZ[53][1](2, 1) = 2 * x2;
    gggXYZ[53][2](0, 0) = 2 * y2;
    gggXYZ[53][2](0, 1) = 4 * x * y;
    gggXYZ[53][2](1, 0) = 4 * x * y;
    gggXYZ[53][2](1, 1) = 2 * x2;
    XYZ[52]             = x * y * z3; // ZZZXY
    gradXYZ[52][0]      = y * z3;
    gradXYZ[52][1]      = x * z3;
    gradXYZ[52][2]      = 3 * x * y * z2;
    hessXYZ[52](0, 1)   = z3;
    hessXYZ[52](0, 2)   = 3 * y * z2;
    hessXYZ[52](1, 0)   = z3;
    hessXYZ[52](1, 2)   = 3 * x * z2;
    hessXYZ[52](2, 0)   = 3 * y * z2;
    hessXYZ[52](2, 1)   = 3 * x * z2;
    hessXYZ[52](2, 2)   = 6 * x * y * z;
    gggXYZ[52][0](1, 2) = 3 * z2;
    gggXYZ[52][0](2, 1) = 3 * z2;
    gggXYZ[52][0](2, 2) = 6 * y * z;
    gggXYZ[52][1](0, 2) = 3 * z2;
    gggXYZ[52][1](2, 0) = 3 * z2;
    gggXYZ[52][1](2, 2) = 6 * x * z;
    gggXYZ[52][2](0, 1) = 3 * z2;
    gggXYZ[52][2](0, 2) = 6 * y * z;
    gggXYZ[52][2](1, 0) = 3 * z2;
    gggXYZ[52][2](1, 2) = 6 * x * z;
    gggXYZ[52][2](2, 0) = 6 * y * z;
    gggXYZ[52][2](2, 1) = 6 * x * z;
    gggXYZ[52][2](2, 2) = 6 * x * y;
    XYZ[51]             = x * y3 * z; // YYYXZ
    gradXYZ[51][0]      = y3 * z;
    gradXYZ[51][1]      = 3 * x * y2 * z;
    gradXYZ[51][2]      = x * y3;
    hessXYZ[51](0, 1)   = 3 * y2 * z;
    hessXYZ[51](0, 2)   = y3;
    hessXYZ[51](1, 0)   = 3 * y2 * z;
    hessXYZ[51](1, 1)   = 6 * x * y * z;
    hessXYZ[51](1, 2)   = 3 * x * y2;
    hessXYZ[51](2, 0)   = y3;
    hessXYZ[51](2, 1)   = 3 * x * y2;
    gggXYZ[51][0](1, 1) = 6 * y * z;
    gggXYZ[51][0](1, 2) = 3 * y2;
    gggXYZ[51][0](2, 1) = 3 * y2;
    gggXYZ[51][1](0, 1) = 6 * y * z;
    gggXYZ[51][1](0, 2) = 3 * y2;
    gggXYZ[51][1](1, 0) = 6 * y * z;
    gggXYZ[51][1](1, 1) = 6 * x * z;
    gggXYZ[51][1](1, 2) = 6 * x * y;
    gggXYZ[51][1](2, 0) = 3 * y2;
    gggXYZ[51][1](2, 1) = 6 * x * y;
    gggXYZ[51][2](0, 1) = 3 * y2;
    gggXYZ[51][2](1, 0) = 3 * y2;
    gggXYZ[51][2](1, 1) = 6 * x * y;
    XYZ[50]             = x3 * y * z; // XXXYZ
    gradXYZ[50][0]      = 3 * x2 * y * z;
    gradXYZ[50][1]      = x3 * z;
    gradXYZ[50][2]      = x3 * y;
    hessXYZ[50](0, 0)   = 6 * x * y * z;
    hessXYZ[50](0, 1)   = 3 * x2 * z;
    hessXYZ[50](0, 2)   = 3 * x2 * y;
    hessXYZ[50](1, 0)   = 3 * x2 * z;
    hessXYZ[50](1, 2)   = x3;
    hessXYZ[50](2, 0)   = 3 * x2 * y;
    hessXYZ[50](2, 1)   = x3;
    gggXYZ[50][0](0, 0) = 6 * y * z;
    gggXYZ[50][0](0, 1) = 6 * x * z;
    gggXYZ[50][0](0, 2) = 6 * x * y;
    gggXYZ[50][0](1, 0) = 6 * x * z;
    gggXYZ[50][0](1, 2) = 3 * x2;
    gggXYZ[50][0](2, 0) = 6 * x * y;
    gggXYZ[50][0](2, 1) = 3 * x2;
    gggXYZ[50][1](0, 0) = 6 * x * z;
    gggXYZ[50][1](0, 2) = 3 * x2;
    gggXYZ[50][1](2, 0) = 3 * x2;
    gggXYZ[50][2](0, 0) = 6 * x * y;
    gggXYZ[50][2](0, 1) = 3 * x2;
    gggXYZ[50][2](1, 0) = 3 * x2;
    XYZ[49]             = y2 * z3; // ZZZYY
    gradXYZ[49][1]      = 2 * y * z3;
    gradXYZ[49][2]      = 3 * y2 * z2;
    hessXYZ[49](1, 1)   = 2 * z3;
    hessXYZ[49](1, 2)   = 6 * y * z2;
    hessXYZ[49](2, 1)   = 6 * y * z2;
    hessXYZ[49](2, 2)   = 6 * y2 * z;
    gggXYZ[49][1](1, 2) = 6 * z2;
    gggXYZ[49][1](2, 1) = 6 * z2;
    gggXYZ[49][1](2, 2) = 12 * y * z;
    gggXYZ[49][2](1, 1) = 6 * z2;
    gggXYZ[49][2](1, 2) = 12 * y * z;
    gggXYZ[49][2](2, 1) = 12 * y * z;
    gggXYZ[49][2](2, 2) = 6 * y2;
    XYZ[48]             = x2 * z3; // ZZZXX
    gradXYZ[48][0]      = 2 * x * z3;
    gradXYZ[48][2]      = 3 * x2 * z2;
    hessXYZ[48](0, 0)   = 2 * z3;
    hessXYZ[48](0, 2)   = 6 * x * z2;
    hessXYZ[48](2, 0)   = 6 * x * z2;
    hessXYZ[48](2, 2)   = 6 * x2 * z;
    gggXYZ[48][0](0, 2) = 6 * z2;
    gggXYZ[48][0](2, 0) = 6 * z2;
    gggXYZ[48][0](2, 2) = 12 * x * z;
    gggXYZ[48][2](0, 0) = 6 * z2;
    gggXYZ[48][2](0, 2) = 12 * x * z;
    gggXYZ[48][2](2, 0) = 12 * x * z;
    gggXYZ[48][2](2, 2) = 6 * x2;
    XYZ[47]             = y3 * z2; // YYYZZ
    gradXYZ[47][1]      = 3 * y2 * z2;
    gradXYZ[47][2]      = 2 * y3 * z;
    hessXYZ[47](1, 1)   = 6 * y * z2;
    hessXYZ[47](1, 2)   = 6 * y2 * z;
    hessXYZ[47](2, 1)   = 6 * y2 * z;
    hessXYZ[47](2, 2)   = 2 * y3;
    gggXYZ[47][1](1, 1) = 6 * z2;
    gggXYZ[47][1](1, 2) = 12 * y * z;
    gggXYZ[47][1](2, 1) = 12 * y * z;
    gggXYZ[47][1](2, 2) = 6 * y2;
    gggXYZ[47][2](1, 1) = 12 * y * z;
    gggXYZ[47][2](1, 2) = 6 * y2;
    gggXYZ[47][2](2, 1) = 6 * y2;
    XYZ[46]             = x2 * y3; // YYYXX
    gradXYZ[46][0]      = 2 * x * y3;
    gradXYZ[46][1]      = 3 * x2 * y2;
    hessXYZ[46](0, 0)   = 2 * y3;
    hessXYZ[46](0, 1)   = 6 * x * y2;
    hessXYZ[46](1, 0)   = 6 * x * y2;
    hessXYZ[46](1, 1)   = 6 * x2 * y;
    gggXYZ[46][0](0, 1) = 6 * y2;
    gggXYZ[46][0](1, 0) = 6 * y2;
    gggXYZ[46][0](1, 1) = 12 * x * y;
    gggXYZ[46][1](0, 0) = 6 * y2;
    gggXYZ[46][1](0, 1) = 12 * x * y;
    gggXYZ[46][1](1, 0) = 12 * x * y;
    gggXYZ[46][1](1, 1) = 6 * x2;
    XYZ[45]             = x3 * z2; // XXXZZ
    gradXYZ[45][0]      = 3 * x2 * z2;
    gradXYZ[45][2]      = 2 * x3 * z;
    hessXYZ[45](0, 0)   = 6 * x * z2;
    hessXYZ[45](0, 2)   = 6 * x2 * z;
    hessXYZ[45](2, 0)   = 6 * x2 * z;
    hessXYZ[45](2, 2)   = 2 * x3;
    gggXYZ[45][0](0, 0) = 6 * z2;
    gggXYZ[45][0](0, 2) = 12 * x * z;
    gggXYZ[45][0](2, 0) = 12 * x * z;
    gggXYZ[45][0](2, 2) = 6 * x2;
    gggXYZ[45][2](0, 0) = 12 * x * z;
    gggXYZ[45][2](0, 2) = 6 * x2;
    gggXYZ[45][2](2, 0) = 6 * x2;
    XYZ[44]             = x3 * y2; // XXXYY
    gradXYZ[44][0]      = 3 * x2 * y2;
    gradXYZ[44][1]      = 2 * x3 * y;
    hessXYZ[44](0, 0)   = 6 * x * y2;
    hessXYZ[44](0, 1)   = 6 * x2 * y;
    hessXYZ[44](1, 0)   = 6 * x2 * y;
    hessXYZ[44](1, 1)   = 2 * x3;
    gggXYZ[44][0](0, 0) = 6 * y2;
    gggXYZ[44][0](0, 1) = 12 * x * y;
    gggXYZ[44][0](1, 0) = 12 * x * y;
    gggXYZ[44][0](1, 1) = 6 * x2;
    gggXYZ[44][1](0, 0) = 12 * x * y;
    gggXYZ[44][1](0, 1) = 6 * x2;
    gggXYZ[44][1](1, 0) = 6 * x2;
    XYZ[43]             = y * z4; // ZZZZY
    gradXYZ[43][1]      = z4;
    gradXYZ[43][2]      = 4 * y * z3;
    hessXYZ[43](1, 2)   = 4 * z3;
    hessXYZ[43](2, 1)   = 4 * z3;
    hessXYZ[43](2, 2)   = 12 * y * z2;
    gggXYZ[43][1](2, 2) = 12 * z2;
    gggXYZ[43][2](1, 2) = 12 * z2;
    gggXYZ[43][2](2, 1) = 12 * z2;
    gggXYZ[43][2](2, 2) = 24 * y * z;
    XYZ[42]             = x * z4; // ZZZZX
    gradXYZ[42][0]      = z4;
    gradXYZ[42][2]      = 4 * x * z3;
    hessXYZ[42](0, 2)   = 4 * z3;
    hessXYZ[42](2, 0)   = 4 * z3;
    hessXYZ[42](2, 2)   = 12 * x * z2;
    gggXYZ[42][0](2, 2) = 12 * z2;
    gggXYZ[42][2](0, 2) = 12 * z2;
    gggXYZ[42][2](2, 0) = 12 * z2;
    gggXYZ[42][2](2, 2) = 24 * x * z;
    XYZ[41]             = y4 * z; // YYYYZ
    gradXYZ[41][1]      = 4 * y3 * z;
    gradXYZ[41][2]      = y4;
    hessXYZ[41](1, 1)   = 12 * y2 * z;
    hessXYZ[41](1, 2)   = 4 * y3;
    hessXYZ[41](2, 1)   = 4 * y3;
    gggXYZ[41][1](1, 1) = 24 * y * z;
    gggXYZ[41][1](1, 2) = 12 * y2;
    gggXYZ[41][1](2, 1) = 12 * y2;
    gggXYZ[41][2](1, 1) = 12 * y2;
    XYZ[40]             = x * y4; // YYYYX
    gradXYZ[40][0]      = y4;
    gradXYZ[40][1]      = 4 * x * y3;
    hessXYZ[40](0, 1)   = 4 * y3;
    hessXYZ[40](1, 0)   = 4 * y3;
    hessXYZ[40](1, 1)   = 12 * x * y2;
    gggXYZ[40][0](1, 1) = 12 * y2;
    gggXYZ[40][1](0, 1) = 12 * y2;
    gggXYZ[40][1](1, 0) = 12 * y2;
    gggXYZ[40][1](1, 1) = 24 * x * y;
    XYZ[39]             = x4 * z; // XXXXZ
    gradXYZ[39][0]      = 4 * x3 * z;
    gradXYZ[39][2]      = x4;
    hessXYZ[39](0, 0)   = 12 * x2 * z;
    hessXYZ[39](0, 2)   = 4 * x3;
    hessXYZ[39](2, 0)   = 4 * x3;
    gggXYZ[39][0](0, 0) = 24 * x * z;
    gggXYZ[39][0](0, 2) = 12 * x2;
    gggXYZ[39][0](2, 0) = 12 * x2;
    gggXYZ[39][2](0, 0) = 12 * x2;
    XYZ[38]             = x4 * y; // XXXXY
    gradXYZ[38][0]      = 4 * x3 * y;
    gradXYZ[38][1]      = x4;
    hessXYZ[38](0, 0)   = 12 * x2 * y;
    hessXYZ[38](0, 1)   = 4 * x3;
    hessXYZ[38](1, 0)   = 4 * x3;
    gggXYZ[38][0](0, 0) = 24 * x * y;
    gggXYZ[38][0](0, 1) = 12 * x2;
    gggXYZ[38][0](1, 0) = 12 * x2;
    gggXYZ[38][1](0, 0) = 12 * x2;
    XYZ[37]             = z * z4; // ZZZZZ
    gradXYZ[37][2]      = 5 * z4;
    hessXYZ[37](2, 2)   = 20 * z3;
    gggXYZ[37][2](2, 2) = 60 * z2;
    XYZ[36]             = y * y4; // YYYYY
    gradXYZ[36][1]      = 5 * y4;
    hessXYZ[36](1, 1)   = 20 * y3;
    gggXYZ[36][1](1, 1) = 60 * y2;
    XYZ[35]             = x * x4; // XXXXX
    gradXYZ[35][0]      = 5 * x4;
    hessXYZ[35](0, 0)   = 20 * x3;
    gggXYZ[35][0](0, 0) = 60 * x2;
  case 4:
    XYZ[34]             = x * y * z2; // ZZXY
    gradXYZ[34][0]      = y * z2;
    gradXYZ[34][1]      = x * z2;
    gradXYZ[34][2]      = 2 * x * y * z;
    hessXYZ[34](0, 1)   = z2;
    hessXYZ[34](0, 2)   = 2 * y * z;
    hessXYZ[34](1, 0)   = z2;
    hessXYZ[34](1, 2)   = 2 * x * z;
    hessXYZ[34](2, 0)   = 2 * y * z;
    hessXYZ[34](2, 1)   = 2 * x * z;
    hessXYZ[34](2, 2)   = 2 * x * y;
    gggXYZ[34][0](1, 2) = 2 * z;
    gggXYZ[34][0](2, 1) = 2 * z;
    gggXYZ[34][0](2, 2) = 2 * y;
    gggXYZ[34][1](0, 2) = 2 * z;
    gggXYZ[34][1](2, 0) = 2 * z;
    gggXYZ[34][1](2, 2) = 2 * x;
    gggXYZ[34][2](0, 1) = 2 * z;
    gggXYZ[34][2](0, 2) = 2 * y;
    gggXYZ[34][2](1, 0) = 2 * z;
    gggXYZ[34][2](1, 2) = 2 * x;
    gggXYZ[34][2](2, 0) = 2 * y;
    gggXYZ[34][2](2, 1) = 2 * x;
    XYZ[33]             = x * y2 * z; // YYXZ
    gradXYZ[33][0]      = y2 * z;
    gradXYZ[33][1]      = 2 * x * y * z;
    gradXYZ[33][2]      = x * y2;
    hessXYZ[33](0, 1)   = 2 * y * z;
    hessXYZ[33](0, 2)   = y2;
    hessXYZ[33](1, 0)   = 2 * y * z;
    hessXYZ[33](1, 1)   = 2 * x * z;
    hessXYZ[33](1, 2)   = 2 * x * y;
    hessXYZ[33](2, 0)   = y2;
    hessXYZ[33](2, 1)   = 2 * x * y;
    gggXYZ[33][0](1, 1) = 2 * z;
    gggXYZ[33][0](1, 2) = 2 * y;
    gggXYZ[33][0](2, 1) = 2 * y;
    gggXYZ[33][1](0, 1) = 2 * z;
    gggXYZ[33][1](0, 2) = 2 * y;
    gggXYZ[33][1](1, 0) = 2 * z;
    gggXYZ[33][1](1, 2) = 2 * x;
    gggXYZ[33][1](2, 0) = 2 * y;
    gggXYZ[33][1](2, 1) = 2 * x;
    gggXYZ[33][2](0, 1) = 2 * y;
    gggXYZ[33][2](1, 0) = 2 * y;
    gggXYZ[33][2](1, 1) = 2 * x;
    XYZ[32]             = x2 * y * z; // XXYZ
    gradXYZ[32][0]      = 2 * x * y * z;
    gradXYZ[32][1]      = x2 * z;
    gradXYZ[32][2]      = x2 * y;
    hessXYZ[32](0, 0)   = 2 * y * z;
    hessXYZ[32](0, 1)   = 2 * x * z;
    hessXYZ[32](0, 2)   = 2 * x * y;
    hessXYZ[32](1, 0)   = 2 * x * z;
    hessXYZ[32](1, 2)   = x2;
    hessXYZ[32](2, 0)   = 2 * x * y;
    hessXYZ[32](2, 1)   = x2;
    gggXYZ[32][0](0, 1) = 2 * z;
    gggXYZ[32][0](0, 2) = 2 * y;
    gggXYZ[32][0](1, 0) = 2 * z;
    gggXYZ[32][0](1, 2) = 2 * x;
    gggXYZ[32][0](2, 0) = 2 * y;
    gggXYZ[32][0](2, 1) = 2 * x;
    gggXYZ[32][1](0, 0) = 2 * z;
    gggXYZ[32][1](0, 2) = 2 * x;
    gggXYZ[32][1](2, 0) = 2 * x;
    gggXYZ[32][2](0, 0) = 2 * y;
    gggXYZ[32][2](0, 1) = 2 * x;
    gggXYZ[32][2](1, 0) = 2 * x;
    XYZ[31]             = y2 * z2; // YYZZ
    gradXYZ[31][1]      = 2 * y * z2;
    gradXYZ[31][2]      = 2 * y2 * z;
    hessXYZ[31](1, 1)   = 2 * z2;
    hessXYZ[31](1, 2)   = 4 * y * z;
    hessXYZ[31](2, 1)   = 4 * y * z;
    hessXYZ[31](2, 2)   = 2 * y2;
    gggXYZ[31][1](1, 2) = 4 * z;
    gggXYZ[31][1](2, 1) = 4 * z;
    gggXYZ[31][1](2, 2) = 4 * y;
    gggXYZ[31][2](1, 1) = 4 * z;
    gggXYZ[31][2](1, 2) = 4 * y;
    gggXYZ[31][2](2, 1) = 4 * y;
    XYZ[30]             = x2 * z2; // XXZZ
    gradXYZ[30][0]      = 2 * x * z2;
    gradXYZ[30][2]      = 2 * x2 * z;
    hessXYZ[30](0, 0)   = 2 * z2;
    hessXYZ[30](0, 2)   = 4 * x * z;
    hessXYZ[30](2, 0)   = 4 * x * z;
    hessXYZ[30](2, 2)   = 2 * x2;
    gggXYZ[30][0](0, 2) = 4 * z;
    gggXYZ[30][0](2, 0) = 4 * z;
    gggXYZ[30][0](2, 2) = 4 * x;
    gggXYZ[30][2](0, 0) = 4 * z;
    gggXYZ[30][2](0, 2) = 4 * x;
    gggXYZ[30][2](2, 0) = 4 * x;
    XYZ[29]             = x2 * y2; // XXYY
    gradXYZ[29][0]      = 2 * x * y2;
    gradXYZ[29][1]      = 2 * x2 * y;
    hessXYZ[29](0, 0)   = 2 * y2;
    hessXYZ[29](0, 1)   = 4 * x * y;
    hessXYZ[29](1, 0)   = 4 * x * y;
    hessXYZ[29](1, 1)   = 2 * x2;
    gggXYZ[29][0](0, 1) = 4 * y;
    gggXYZ[29][0](1, 0) = 4 * y;
    gggXYZ[29][0](1, 1) = 4 * x;
    gggXYZ[29][1](0, 0) = 4 * y;
    gggXYZ[29][1](0, 1) = 4 * x;
    gggXYZ[29][1](1, 0) = 4 * x;
    XYZ[28]             = y * z3; // ZZZY
    gradXYZ[28][1]      = z3;
    gradXYZ[28][2]      = 3 * y * z2;
    hessXYZ[28](1, 2)   = 3 * z2;
    hessXYZ[28](2, 1)   = 3 * z2;
    hessXYZ[28](2, 2)   = 6 * y * z;
    gggXYZ[28][1](2, 2) = 6 * z;
    gggXYZ[28][2](1, 2) = 6 * z;
    gggXYZ[28][2](2, 1) = 6 * z;
    gggXYZ[28][2](2, 2) = 6 * y;
    XYZ[27]             = x * z3; // ZZZX
    gradXYZ[27][0]      = z3;
    gradXYZ[27][2]      = 3 * x * z2;
    hessXYZ[27](0, 2)   = 3 * z2;
    hessXYZ[27](2, 0)   = 3 * z2;
    hessXYZ[27](2, 2)   = 6 * x * z;
    gggXYZ[27][0](2, 2) = 6 * z;
    gggXYZ[27][2](0, 2) = 6 * z;
    gggXYZ[27][2](2, 0) = 6 * z;
    gggXYZ[27][2](2, 2) = 6 * x;
    XYZ[26]             = y3 * z; // YYYZ
    gradXYZ[26][1]      = 3 * y2 * z;
    gradXYZ[26][2]      = y3;
    hessXYZ[26](1, 1)   = 6 * y * z;
    hessXYZ[26](1, 2)   = 3 * y2;
    hessXYZ[26](2, 1)   = 3 * y2;
    gggXYZ[26][1](1, 1) = 6 * z;
    gggXYZ[26][1](1, 2) = 6 * y;
    gggXYZ[26][1](2, 1) = 6 * y;
    gggXYZ[26][2](1, 1) = 6 * y;
    XYZ[25]             = x * y3; // YYYX
    gradXYZ[25][0]      = y3;
    gradXYZ[25][1]      = 3 * x * y2;
    hessXYZ[25](0, 1)   = 3 * y2;
    hessXYZ[25](1, 0)   = 3 * y2;
    hessXYZ[25](1, 1)   = 6 * x * y;
    gggXYZ[25][0](1, 1) = 6 * y;
    gggXYZ[25][1](0, 1) = 6 * y;
    gggXYZ[25][1](1, 0) = 6 * y;
    gggXYZ[25][1](1, 1) = 6 * x;
    XYZ[24]             = x3 * z; // XXXZ
    gradXYZ[24][0]      = 3 * x2 * z;
    gradXYZ[24][2]      = x3;
    hessXYZ[24](0, 0)   = 6 * x * z;
    hessXYZ[24](0, 2)   = 3 * x2;
    hessXYZ[24](2, 0)   = 3 * x2;
    gggXYZ[24][0](0, 0) = 6 * z;
    gggXYZ[24][0](0, 2) = 6 * x;
    gggXYZ[24][0](2, 0) = 6 * x;
    gggXYZ[24][2](0, 0) = 6 * x;
    XYZ[23]             = x3 * y; // XXXY
    gradXYZ[23][0]      = 3 * x2 * y;
    gradXYZ[23][1]      = x3;
    hessXYZ[23](0, 0)   = 6 * x * y;
    hessXYZ[23](0, 1)   = 3 * x2;
    hessXYZ[23](1, 0)   = 3 * x2;
    gggXYZ[23][0](0, 0) = 6 * y;
    gggXYZ[23][0](0, 1) = 6 * x;
    gggXYZ[23][0](1, 0) = 6 * x;
    gggXYZ[23][1](0, 0) = 6 * x;
    XYZ[22]             = z4; // ZZZZ
    gradXYZ[22][2]      = 4 * z3;
    hessXYZ[22](2, 2)   = 12 * z2;
    gggXYZ[22][2](2, 2) = 24 * z;
    XYZ[21]             = y4; // YYYY
    gradXYZ[21][1]      = 4 * y3;
    hessXYZ[21](1, 1)   = 12 * y2;
    gggXYZ[21][1](1, 1) = 24 * y;
    XYZ[20]             = x4; // XXXX
    gradXYZ[20][0]      = 4 * x3;
    hessXYZ[20](0, 0)   = 12 * x2;
    gggXYZ[20][0](0, 0) = 24 * x;
  case 3:
    XYZ[19]             = x * y * z; // XYZ
    gradXYZ[19][0]      = y * z;
    gradXYZ[19][1]      = x * z;
    gradXYZ[19][2]      = x * y;
    hessXYZ[19](0, 1)   = z;
    hessXYZ[19](0, 2)   = y;
    hessXYZ[19](1, 0)   = z;
    hessXYZ[19](1, 2)   = x;
    hessXYZ[19](2, 0)   = y;
    hessXYZ[19](2, 1)   = x;
    gggXYZ[19][0](1, 2) = 1;
    gggXYZ[19][0](2, 1) = 1;
    gggXYZ[19][1](0, 2) = 1;
    gggXYZ[19][1](2, 0) = 1;
    gggXYZ[19][2](0, 1) = 1;
    gggXYZ[19][2](1, 0) = 1;
    XYZ[18]             = y * z2; // ZZY
    gradXYZ[18][1]      = z2;
    gradXYZ[18][2]      = 2 * y * z;
    hessXYZ[18](1, 2)   = 2 * z;
    hessXYZ[18](2, 1)   = 2 * z;
    hessXYZ[18](2, 2)   = 2 * y;
    gggXYZ[18][1](2, 2) = 2;
    gggXYZ[18][2](1, 2) = 2;
    gggXYZ[18][2](2, 1) = 2;
    XYZ[17]             = x * z2; // ZZX
    gradXYZ[17][0]      = z2;
    gradXYZ[17][2]      = 2 * x * z;
    hessXYZ[17](0, 2)   = 2 * z;
    hessXYZ[17](2, 0)   = 2 * z;
    hessXYZ[17](2, 2)   = 2 * x;
    gggXYZ[17][0](2, 2) = 2;
    gggXYZ[17][2](0, 2) = 2;
    gggXYZ[17][2](2, 0) = 2;
    XYZ[16]             = y2 * z; // YYZ
    gradXYZ[16][1]      = 2 * y * z;
    gradXYZ[16][2]      = y2;
    hessXYZ[16](1, 1)   = 2 * z;
    hessXYZ[16](1, 2)   = 2 * y;
    hessXYZ[16](2, 1)   = 2 * y;
    gggXYZ[16][1](1, 2) = 2;
    gggXYZ[16][1](2, 1) = 2;
    gggXYZ[16][2](1, 1) = 2;
    XYZ[15]             = x * y2; // YYX
    gradXYZ[15][0]      = y2;
    gradXYZ[15][1]      = 2 * x * y;
    hessXYZ[15](0, 1)   = 2 * y;
    hessXYZ[15](1, 0)   = 2 * y;
    hessXYZ[15](1, 1)   = 2 * x;
    gggXYZ[15][0](1, 1) = 2;
    gggXYZ[15][1](0, 1) = 2;
    gggXYZ[15][1](1, 0) = 2;
    XYZ[14]             = x2 * z; // XXZ
    gradXYZ[14][0]      = 2 * x * z;
    gradXYZ[14][2]      = x2;
    hessXYZ[14](0, 0)   = 2 * z;
    hessXYZ[14](0, 2)   = 2 * x;
    hessXYZ[14](2, 0)   = 2 * x;
    gggXYZ[14][0](0, 2) = 2;
    gggXYZ[14][0](2, 0) = 2;
    gggXYZ[14][2](0, 0) = 2;
    XYZ[13]             = x2 * y; // XXY
    gradXYZ[13][0]      = 2 * x * y;
    gradXYZ[13][1]      = x2;
    hessXYZ[13](0, 0)   = 2 * y;
    hessXYZ[13](0, 1)   = 2 * x;
    hessXYZ[13](1, 0)   = 2 * x;
    gggXYZ[13][0](0, 1) = 2;
    gggXYZ[13][0](1, 0) = 2;
    gggXYZ[13][1](0, 0) = 2;
    XYZ[12]             = z3; // ZZZ
    gradXYZ[12][2]      = 3 * z2;
    hessXYZ[12](2, 2)   = 6 * z;
    gggXYZ[12][2](2, 2) = 6;
    XYZ[11]             = y3; // YYY
    gradXYZ[11][1]      = 3 * y2;
    hessXYZ[11](1, 1)   = 6 * y;
    gggXYZ[11][1](1, 1) = 6;
    XYZ[10]             = x3; // XXX
    gradXYZ[10][0]      = 3 * x2;
    hessXYZ[10](0, 0)   = 6 * x;
    gggXYZ[10][0](0, 0) = 6;
  case 2:
    XYZ[9]           = y * z; // YZ
    gradXYZ[9][1]    = z;
    gradXYZ[9][2]    = y;
    hessXYZ[9](1, 2) = 1;
    hessXYZ[9](2, 1) = 1;
    XYZ[8]           = x * z; // XZ
    gradXYZ[8][0]    = z;
    gradXYZ[8][2]    = x;
    hessXYZ[8](0, 2) = 1;
    hessXYZ[8](2, 0) = 1;
    XYZ[7]           = x * y; // XY
    gradXYZ[7][0]    = y;
    gradXYZ[7][1]    = x;
    hessXYZ[7](0, 1) = 1;
    hessXYZ[7](1, 0) = 1;
    XYZ[6]           = z2; // ZZ
    gradXYZ[6][2]    = 2 * z;
    hessXYZ[6](2, 2) = 2;
    XYZ[5]           = y2; // YY
    gradXYZ[5][1]    = 2 * y;
    hessXYZ[5](1, 1) = 2;
    XYZ[4]           = x2; // XX
    gradXYZ[4][0]    = 2 * x;
    hessXYZ[4](0, 0) = 2;
  case 1:
    XYZ[3]        = z; // Z
    gradXYZ[3][2] = 1;
    XYZ[2]        = y; // Y
    gradXYZ[2][1] = 1;
    XYZ[1]        = x; // X
    gradXYZ[1][0] = 1;
  case 0:
    XYZ[0] = 1; // S
  }

  for (int i = 0; i < ntot; i++)
    XYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    gradXYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
    hessXYZ[i] *= NormFactor[i];
  for (int i = 0; i < ntot; i++)
  {
    gggXYZ[i][0] *= NormFactor[i];
    gggXYZ[i][1] *= NormFactor[i];
    gggXYZ[i][2] *= NormFactor[i];
  }
}


template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::evaluateThirdDerivOnly(const Point_t& p)
{
  value_type x = p[0], y = p[1], z = p[2];
  value_type x2 = x * x, y2 = y * y, z2 = z * z;
  value_type x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  int ntot = XYZ.size();
  for (int i = 0; i < ntot; i++)
  {
    gggXYZ[i][0] = 0.0;
    gggXYZ[i][1] = 0.0;
    gggXYZ[i][2] = 0.0;
  }

  switch (Lmax)
  {
  case 6:
    gggXYZ[83][0](0, 1) = 4 * y * z2;
    gggXYZ[83][0](0, 2) = 4 * y2 * z;
    gggXYZ[83][0](1, 0) = 4 * y * z2;
    gggXYZ[83][0](1, 1) = 4 * x * z2;
    gggXYZ[83][0](1, 2) = 8 * x * y * z;
    gggXYZ[83][0](2, 0) = 4 * y2 * z;
    gggXYZ[83][0](2, 1) = 8 * x * y * z;
    gggXYZ[83][0](2, 2) = 4 * x * y2;
    gggXYZ[83][1](0, 0) = 4 * y * z2;
    gggXYZ[83][1](0, 1) = 4 * x * z2;
    gggXYZ[83][1](0, 2) = 8 * x * y * z;
    gggXYZ[83][1](1, 0) = 4 * x * z2;
    gggXYZ[83][1](1, 2) = 4 * x2 * z;
    gggXYZ[83][1](2, 0) = 8 * x * y * z;
    gggXYZ[83][1](2, 1) = 4 * x2 * z;
    gggXYZ[83][1](2, 2) = 4 * x2 * y;
    gggXYZ[83][2](0, 0) = 4 * y2 * z;
    gggXYZ[83][2](0, 1) = 8 * x * y * z;
    gggXYZ[83][2](0, 2) = 4 * x * y2;
    gggXYZ[83][2](1, 0) = 8 * x * y * z;
    gggXYZ[83][2](1, 1) = 4 * x2 * z;
    gggXYZ[83][2](1, 2) = 4 * x2 * y;
    gggXYZ[83][2](2, 0) = 4 * x * y2;
    gggXYZ[83][2](2, 1) = 4 * x2 * y;
    gggXYZ[82][0](1, 1) = 2 * z3;
    gggXYZ[82][0](1, 2) = 6 * y * z2;
    gggXYZ[82][0](2, 1) = 6 * y * z2;
    gggXYZ[82][0](2, 2) = 6 * y2 * z;
    gggXYZ[82][1](0, 1) = 2 * z3;
    gggXYZ[82][1](0, 2) = 6 * y * z2;
    gggXYZ[82][1](1, 0) = 2 * z3;
    gggXYZ[82][1](1, 2) = 6 * x * z2;
    gggXYZ[82][1](2, 0) = 6 * y * z2;
    gggXYZ[82][1](2, 1) = 6 * x * z2;
    gggXYZ[82][1](2, 2) = 12 * x * y * z;
    gggXYZ[82][2](0, 1) = 6 * y * z2;
    gggXYZ[82][2](0, 2) = 6 * y2 * z;
    gggXYZ[82][2](1, 0) = 6 * y * z2;
    gggXYZ[82][2](1, 1) = 6 * x * z2;
    gggXYZ[82][2](1, 2) = 12 * x * y * z;
    gggXYZ[82][2](2, 0) = 6 * y2 * z;
    gggXYZ[82][2](2, 1) = 12 * x * y * z;
    gggXYZ[82][2](2, 2) = 6 * x * y2;
    gggXYZ[81][0](0, 1) = 2 * z3;
    gggXYZ[81][0](0, 2) = 6 * y * z2;
    gggXYZ[81][0](1, 0) = 2 * z3;
    gggXYZ[81][0](1, 2) = 6 * x * z2;
    gggXYZ[81][0](2, 0) = 6 * y * z2;
    gggXYZ[81][0](2, 1) = 6 * x * z2;
    gggXYZ[81][0](2, 2) = 12 * x * y * z;
    gggXYZ[81][1](0, 0) = 2 * z3;
    gggXYZ[81][1](0, 2) = 6 * x * z2;
    gggXYZ[81][1](2, 0) = 6 * x * z2;
    gggXYZ[81][1](2, 2) = 6 * x2 * z;
    gggXYZ[81][2](0, 0) = 6 * y * z2;
    gggXYZ[81][2](0, 1) = 6 * x * z2;
    gggXYZ[81][2](0, 2) = 12 * x * y * z;
    gggXYZ[81][2](1, 0) = 6 * x * z2;
    gggXYZ[81][2](1, 2) = 6 * x2 * z;
    gggXYZ[81][2](2, 0) = 12 * x * y * z;
    gggXYZ[81][2](2, 1) = 6 * x2 * z;
    gggXYZ[81][2](2, 2) = 6 * x2 * y;
    gggXYZ[80][0](1, 1) = 6 * y * z2;
    gggXYZ[80][0](1, 2) = 6 * y2 * z;
    gggXYZ[80][0](2, 1) = 6 * y2 * z;
    gggXYZ[80][0](2, 2) = 2 * y3;
    gggXYZ[80][1](0, 1) = 6 * y * z2;
    gggXYZ[80][1](0, 2) = 6 * y2 * z;
    gggXYZ[80][1](1, 0) = 6 * y * z2;
    gggXYZ[80][1](1, 1) = 6 * x * z2;
    gggXYZ[80][1](1, 2) = 12 * x * y * z;
    gggXYZ[80][1](2, 0) = 6 * y2 * z;
    gggXYZ[80][1](2, 1) = 12 * x * y * z;
    gggXYZ[80][1](2, 2) = 6 * x * y2;
    gggXYZ[80][2](0, 1) = 6 * y2 * z;
    gggXYZ[80][2](0, 2) = 2 * y3;
    gggXYZ[80][2](1, 0) = 6 * y2 * z;
    gggXYZ[80][2](1, 1) = 12 * x * y * z;
    gggXYZ[80][2](1, 2) = 6 * x * y2;
    gggXYZ[80][2](2, 0) = 2 * y3;
    gggXYZ[80][2](2, 1) = 6 * x * y2;
    gggXYZ[79][0](0, 1) = 6 * y2 * z;
    gggXYZ[79][0](0, 2) = 2 * y3;
    gggXYZ[79][0](1, 0) = 6 * y2 * z;
    gggXYZ[79][0](1, 1) = 12 * x * y * z;
    gggXYZ[79][0](1, 2) = 6 * x * y2;
    gggXYZ[79][0](2, 0) = 2 * y3;
    gggXYZ[79][0](2, 1) = 6 * x * y2;
    gggXYZ[79][1](0, 0) = 6 * y2 * z;
    gggXYZ[79][1](0, 1) = 12 * x * y * z;
    gggXYZ[79][1](0, 2) = 6 * x * y2;
    gggXYZ[79][1](1, 0) = 12 * x * y * z;
    gggXYZ[79][1](1, 1) = 6 * x2 * z;
    gggXYZ[79][1](1, 2) = 6 * x2 * y;
    gggXYZ[79][1](2, 0) = 6 * x * y2;
    gggXYZ[79][1](2, 1) = 6 * x2 * y;
    gggXYZ[79][2](0, 0) = 2 * y3;
    gggXYZ[79][2](0, 1) = 6 * x * y2;
    gggXYZ[79][2](1, 0) = 6 * x * y2;
    gggXYZ[79][2](1, 1) = 6 * x2 * y;
    gggXYZ[78][0](0, 0) = 6 * y * z2;
    gggXYZ[78][0](0, 1) = 6 * x * z2;
    gggXYZ[78][0](0, 2) = 12 * x * y * z;
    gggXYZ[78][0](1, 0) = 6 * x * z2;
    gggXYZ[78][0](1, 2) = 6 * x2 * z;
    gggXYZ[78][0](2, 0) = 12 * x * y * z;
    gggXYZ[78][0](2, 1) = 6 * x2 * z;
    gggXYZ[78][0](2, 2) = 6 * x2 * y;
    gggXYZ[78][1](0, 0) = 6 * x * z2;
    gggXYZ[78][1](0, 2) = 6 * x2 * z;
    gggXYZ[78][1](2, 0) = 6 * x2 * z;
    gggXYZ[78][1](2, 2) = 2 * x3;
    gggXYZ[78][2](0, 0) = 12 * x * y * z;
    gggXYZ[78][2](0, 1) = 6 * x2 * z;
    gggXYZ[78][2](0, 2) = 6 * x2 * y;
    gggXYZ[78][2](1, 0) = 6 * x2 * z;
    gggXYZ[78][2](1, 2) = 2 * x3;
    gggXYZ[78][2](2, 0) = 6 * x2 * y;
    gggXYZ[78][2](2, 1) = 2 * x3;
    gggXYZ[77][0](0, 0) = 6 * y2 * z;
    gggXYZ[77][0](0, 1) = 12 * x * y * z;
    gggXYZ[77][0](0, 2) = 6 * x * y2;
    gggXYZ[77][0](1, 0) = 12 * x * y * z;
    gggXYZ[77][0](1, 1) = 6 * x2 * z;
    gggXYZ[77][0](1, 2) = 6 * x2 * y;
    gggXYZ[77][0](2, 0) = 6 * x * y2;
    gggXYZ[77][0](2, 1) = 6 * x2 * y;
    gggXYZ[77][1](0, 0) = 12 * x * y * z;
    gggXYZ[77][1](0, 1) = 6 * x2 * z;
    gggXYZ[77][1](0, 2) = 6 * x2 * y;
    gggXYZ[77][1](1, 0) = 6 * x2 * z;
    gggXYZ[77][1](1, 2) = 2 * x3;
    gggXYZ[77][1](2, 0) = 6 * x2 * y;
    gggXYZ[77][1](2, 1) = 2 * x3;
    gggXYZ[77][2](0, 0) = 6 * x * y2;
    gggXYZ[77][2](0, 1) = 6 * x2 * y;
    gggXYZ[77][2](1, 0) = 6 * x2 * y;
    gggXYZ[77][2](1, 1) = 2 * x3;
    gggXYZ[76][1](1, 1) = 6 * z3;
    gggXYZ[76][1](1, 2) = 18 * y * z2;
    gggXYZ[76][1](2, 1) = 18 * y * z2;
    gggXYZ[76][1](2, 2) = 18 * y2 * z;
    gggXYZ[76][2](1, 1) = 18 * y * z2;
    gggXYZ[76][2](1, 2) = 18 * y2 * z;
    gggXYZ[76][2](2, 1) = 18 * y2 * z;
    gggXYZ[76][2](2, 2) = 6 * y3;
    gggXYZ[75][0](0, 0) = 6 * z3;
    gggXYZ[75][0](0, 2) = 18 * x * z2;
    gggXYZ[75][0](2, 0) = 18 * x * z2;
    gggXYZ[75][0](2, 2) = 18 * x2 * z;
    gggXYZ[75][2](0, 0) = 18 * x * z2;
    gggXYZ[75][2](0, 2) = 18 * x2 * z;
    gggXYZ[75][2](2, 0) = 18 * x2 * z;
    gggXYZ[75][2](2, 2) = 6 * x3;
    gggXYZ[74][0](0, 0) = 6 * y3;
    gggXYZ[74][0](0, 1) = 18 * x * y2;
    gggXYZ[74][0](1, 0) = 18 * x * y2;
    gggXYZ[74][0](1, 1) = 18 * x2 * y;
    gggXYZ[74][1](0, 0) = 18 * x * y2;
    gggXYZ[74][1](0, 1) = 18 * x2 * y;
    gggXYZ[74][1](1, 0) = 18 * x2 * y;
    gggXYZ[74][1](1, 1) = 6 * x3;
    gggXYZ[73][0](1, 2) = 4 * z3;
    gggXYZ[73][0](2, 1) = 4 * z3;
    gggXYZ[73][0](2, 2) = 12 * y * z2;
    gggXYZ[73][1](0, 2) = 4 * z3;
    gggXYZ[73][1](2, 0) = 4 * z3;
    gggXYZ[73][1](2, 2) = 12 * x * z2;
    gggXYZ[73][2](0, 1) = 4 * z3;
    gggXYZ[73][2](0, 2) = 12 * y * z2;
    gggXYZ[73][2](1, 0) = 4 * z3;
    gggXYZ[73][2](1, 2) = 12 * x * z2;
    gggXYZ[73][2](2, 0) = 12 * y * z2;
    gggXYZ[73][2](2, 1) = 12 * x * z2;
    gggXYZ[73][2](2, 2) = 24 * x * y * z;
    gggXYZ[72][0](1, 1) = 12 * y2 * z;
    gggXYZ[72][0](1, 2) = 4 * y3;
    gggXYZ[72][0](2, 1) = 4 * y3;
    gggXYZ[72][1](0, 1) = 12 * y2 * z;
    gggXYZ[72][1](0, 2) = 4 * y3;
    gggXYZ[72][1](1, 0) = 12 * y2 * z;
    gggXYZ[72][1](1, 1) = 24 * x * y * z;
    gggXYZ[72][1](1, 2) = 12 * x * y2;
    gggXYZ[72][1](2, 0) = 4 * y3;
    gggXYZ[72][1](2, 1) = 12 * x * y2;
    gggXYZ[72][2](0, 1) = 4 * y3;
    gggXYZ[72][2](1, 0) = 4 * y3;
    gggXYZ[72][2](1, 1) = 12 * x * y2;
    gggXYZ[71][0](0, 0) = 24 * x * y * z;
    gggXYZ[71][0](0, 1) = 12 * x2 * z;
    gggXYZ[71][0](0, 2) = 12 * x2 * y;
    gggXYZ[71][0](1, 0) = 12 * x2 * z;
    gggXYZ[71][0](1, 2) = 4 * x3;
    gggXYZ[71][0](2, 0) = 12 * x2 * y;
    gggXYZ[71][0](2, 1) = 4 * x3;
    gggXYZ[71][1](0, 0) = 12 * x2 * z;
    gggXYZ[71][1](0, 2) = 4 * x3;
    gggXYZ[71][1](2, 0) = 4 * x3;
    gggXYZ[71][2](0, 0) = 12 * x2 * y;
    gggXYZ[71][2](0, 1) = 4 * x3;
    gggXYZ[71][2](1, 0) = 4 * x3;
    gggXYZ[70][1](1, 2) = 8 * z3;
    gggXYZ[70][1](2, 1) = 8 * z3;
    gggXYZ[70][1](2, 2) = 24 * y * z2;
    gggXYZ[70][2](1, 1) = 8 * z3;
    gggXYZ[70][2](1, 2) = 24 * y * z2;
    gggXYZ[70][2](2, 1) = 24 * y * z2;
    gggXYZ[70][2](2, 2) = 24 * y2 * z;
    gggXYZ[69][0](0, 2) = 8 * z3;
    gggXYZ[69][0](2, 0) = 8 * z3;
    gggXYZ[69][0](2, 2) = 24 * x * z2;
    gggXYZ[69][2](0, 0) = 8 * z3;
    gggXYZ[69][2](0, 2) = 24 * x * z2;
    gggXYZ[69][2](2, 0) = 24 * x * z2;
    gggXYZ[69][2](2, 2) = 24 * x2 * z;
    gggXYZ[68][1](1, 1) = 24 * y * z2;
    gggXYZ[68][1](1, 2) = 24 * y2 * z;
    gggXYZ[68][1](2, 1) = 24 * y2 * z;
    gggXYZ[68][1](2, 2) = 8 * y3;
    gggXYZ[68][2](1, 1) = 24 * y2 * z;
    gggXYZ[68][2](1, 2) = 8 * y3;
    gggXYZ[68][2](2, 1) = 8 * y3;
    gggXYZ[67][0](0, 1) = 8 * y3;
    gggXYZ[67][0](1, 0) = 8 * y3;
    gggXYZ[67][0](1, 1) = 24 * x * y2;
    gggXYZ[67][1](0, 0) = 8 * y3;
    gggXYZ[67][1](0, 1) = 24 * x * y2;
    gggXYZ[67][1](1, 0) = 24 * x * y2;
    gggXYZ[67][1](1, 1) = 24 * x2 * y;
    gggXYZ[66][0](0, 0) = 24 * x * z2;
    gggXYZ[66][0](0, 2) = 24 * x2 * z;
    gggXYZ[66][0](2, 0) = 24 * x2 * z;
    gggXYZ[66][0](2, 2) = 8 * x3;
    gggXYZ[66][2](0, 0) = 24 * x2 * z;
    gggXYZ[66][2](0, 2) = 8 * x3;
    gggXYZ[66][2](2, 0) = 8 * x3;
    gggXYZ[65][0](0, 0) = 24 * x * y2;
    gggXYZ[65][0](0, 1) = 24 * x2 * y;
    gggXYZ[65][0](1, 0) = 24 * x2 * y;
    gggXYZ[65][0](1, 1) = 8 * x3;
    gggXYZ[65][1](0, 0) = 24 * x2 * y;
    gggXYZ[65][1](0, 1) = 8 * x3;
    gggXYZ[65][1](1, 0) = 8 * x3;
    gggXYZ[64][1](2, 2) = 20 * z3;
    gggXYZ[64][2](1, 2) = 20 * z3;
    gggXYZ[64][2](2, 1) = 20 * z3;
    gggXYZ[64][2](2, 2) = 60 * y * z2;
    gggXYZ[63][0](2, 2) = 20 * z3;
    gggXYZ[63][2](0, 2) = 20 * z3;
    gggXYZ[63][2](2, 0) = 20 * z3;
    gggXYZ[63][2](2, 2) = 60 * x * z2;
    gggXYZ[62][1](1, 1) = 60 * y2 * z;
    gggXYZ[62][1](1, 2) = 20 * y3;
    gggXYZ[62][1](2, 1) = 20 * y3;
    gggXYZ[62][2](1, 1) = 20 * y3;
    gggXYZ[61][0](1, 1) = 20 * y3;
    gggXYZ[61][1](0, 1) = 20 * y3;
    gggXYZ[61][1](1, 0) = 20 * y3;
    gggXYZ[61][1](1, 1) = 60 * x * y2;
    gggXYZ[60][0](0, 0) = 60 * x2 * z;
    gggXYZ[60][0](0, 2) = 20 * x3;
    gggXYZ[60][0](2, 0) = 20 * x3;
    gggXYZ[60][2](0, 0) = 20 * x3;
    gggXYZ[59][0](0, 0) = 60 * x2 * y;
    gggXYZ[59][0](0, 1) = 20 * x3;
    gggXYZ[59][0](1, 0) = 20 * x3;
    gggXYZ[59][1](0, 0) = 20 * x3;
    gggXYZ[58][2](2, 2) = 120 * z3;
    gggXYZ[57][1](1, 1) = 120 * y3;
    gggXYZ[56][0](0, 0) = 120 * x3;
  case 5:
    gggXYZ[55][0](1, 1) = 2 * z2;
    gggXYZ[55][0](1, 2) = 4 * y * z;
    gggXYZ[55][0](2, 1) = 4 * y * z;
    gggXYZ[55][0](2, 2) = 2 * y2;
    gggXYZ[55][1](0, 1) = 2 * z2;
    gggXYZ[55][1](0, 2) = 4 * y * z;
    gggXYZ[55][1](1, 0) = 2 * z2;
    gggXYZ[55][1](1, 2) = 4 * x * z;
    gggXYZ[55][1](2, 0) = 4 * y * z;
    gggXYZ[55][1](2, 1) = 4 * x * z;
    gggXYZ[55][1](2, 2) = 4 * x * y;
    gggXYZ[55][2](0, 1) = 4 * y * z;
    gggXYZ[55][2](0, 2) = 2 * y2;
    gggXYZ[55][2](1, 0) = 4 * y * z;
    gggXYZ[55][2](1, 1) = 4 * x * z;
    gggXYZ[55][2](1, 2) = 4 * x * y;
    gggXYZ[55][2](2, 0) = 2 * y2;
    gggXYZ[55][2](2, 1) = 4 * x * y;
    gggXYZ[54][0](0, 1) = 2 * z2;
    gggXYZ[54][0](0, 2) = 4 * y * z;
    gggXYZ[54][0](1, 0) = 2 * z2;
    gggXYZ[54][0](1, 2) = 4 * x * z;
    gggXYZ[54][0](2, 0) = 4 * y * z;
    gggXYZ[54][0](2, 1) = 4 * x * z;
    gggXYZ[54][0](2, 2) = 4 * x * y;
    gggXYZ[54][1](0, 0) = 2 * z2;
    gggXYZ[54][1](0, 2) = 4 * x * z;
    gggXYZ[54][1](2, 0) = 4 * x * z;
    gggXYZ[54][1](2, 2) = 2 * x2;
    gggXYZ[54][2](0, 0) = 4 * y * z;
    gggXYZ[54][2](0, 1) = 4 * x * z;
    gggXYZ[54][2](0, 2) = 4 * x * y;
    gggXYZ[54][2](1, 0) = 4 * x * z;
    gggXYZ[54][2](1, 2) = 2 * x2;
    gggXYZ[54][2](2, 0) = 4 * x * y;
    gggXYZ[54][2](2, 1) = 2 * x2;
    gggXYZ[53][0](0, 1) = 4 * y * z;
    gggXYZ[53][0](0, 2) = 2 * y2;
    gggXYZ[53][0](1, 0) = 4 * y * z;
    gggXYZ[53][0](1, 1) = 4 * x * z;
    gggXYZ[53][0](1, 2) = 4 * x * y;
    gggXYZ[53][0](2, 0) = 2 * y2;
    gggXYZ[53][0](2, 1) = 4 * x * y;
    gggXYZ[53][1](0, 0) = 4 * y * z;
    gggXYZ[53][1](0, 1) = 4 * x * z;
    gggXYZ[53][1](0, 2) = 4 * x * y;
    gggXYZ[53][1](1, 0) = 4 * x * z;
    gggXYZ[53][1](1, 2) = 2 * x2;
    gggXYZ[53][1](2, 0) = 4 * x * y;
    gggXYZ[53][1](2, 1) = 2 * x2;
    gggXYZ[53][2](0, 0) = 2 * y2;
    gggXYZ[53][2](0, 1) = 4 * x * y;
    gggXYZ[53][2](1, 0) = 4 * x * y;
    gggXYZ[53][2](1, 1) = 2 * x2;
    gggXYZ[52][0](1, 2) = 3 * z2;
    gggXYZ[52][0](2, 1) = 3 * z2;
    gggXYZ[52][0](2, 2) = 6 * y * z;
    gggXYZ[52][1](0, 2) = 3 * z2;
    gggXYZ[52][1](2, 0) = 3 * z2;
    gggXYZ[52][1](2, 2) = 6 * x * z;
    gggXYZ[52][2](0, 1) = 3 * z2;
    gggXYZ[52][2](0, 2) = 6 * y * z;
    gggXYZ[52][2](1, 0) = 3 * z2;
    gggXYZ[52][2](1, 2) = 6 * x * z;
    gggXYZ[52][2](2, 0) = 6 * y * z;
    gggXYZ[52][2](2, 1) = 6 * x * z;
    gggXYZ[52][2](2, 2) = 6 * x * y;
    gggXYZ[51][0](1, 1) = 6 * y * z;
    gggXYZ[51][0](1, 2) = 3 * y2;
    gggXYZ[51][0](2, 1) = 3 * y2;
    gggXYZ[51][1](0, 1) = 6 * y * z;
    gggXYZ[51][1](0, 2) = 3 * y2;
    gggXYZ[51][1](1, 0) = 6 * y * z;
    gggXYZ[51][1](1, 1) = 6 * x * z;
    gggXYZ[51][1](1, 2) = 6 * x * y;
    gggXYZ[51][1](2, 0) = 3 * y2;
    gggXYZ[51][1](2, 1) = 6 * x * y;
    gggXYZ[51][2](0, 1) = 3 * y2;
    gggXYZ[51][2](1, 0) = 3 * y2;
    gggXYZ[51][2](1, 1) = 6 * x * y;
    gggXYZ[50][0](0, 0) = 6 * y * z;
    gggXYZ[50][0](0, 1) = 6 * x * z;
    gggXYZ[50][0](0, 2) = 6 * x * y;
    gggXYZ[50][0](1, 0) = 6 * x * z;
    gggXYZ[50][0](1, 2) = 3 * x2;
    gggXYZ[50][0](2, 0) = 6 * x * y;
    gggXYZ[50][0](2, 1) = 3 * x2;
    gggXYZ[50][1](0, 0) = 6 * x * z;
    gggXYZ[50][1](0, 2) = 3 * x2;
    gggXYZ[50][1](2, 0) = 3 * x2;
    gggXYZ[50][2](0, 0) = 6 * x * y;
    gggXYZ[50][2](0, 1) = 3 * x2;
    gggXYZ[50][2](1, 0) = 3 * x2;
    gggXYZ[49][1](1, 2) = 6 * z2;
    gggXYZ[49][1](2, 1) = 6 * z2;
    gggXYZ[49][1](2, 2) = 12 * y * z;
    gggXYZ[49][2](1, 1) = 6 * z2;
    gggXYZ[49][2](1, 2) = 12 * y * z;
    gggXYZ[49][2](2, 1) = 12 * y * z;
    gggXYZ[49][2](2, 2) = 6 * y2;
    gggXYZ[48][0](0, 2) = 6 * z2;
    gggXYZ[48][0](2, 0) = 6 * z2;
    gggXYZ[48][0](2, 2) = 12 * x * z;
    gggXYZ[48][2](0, 0) = 6 * z2;
    gggXYZ[48][2](0, 2) = 12 * x * z;
    gggXYZ[48][2](2, 0) = 12 * x * z;
    gggXYZ[48][2](2, 2) = 6 * x2;
    gggXYZ[47][1](1, 1) = 6 * z2;
    gggXYZ[47][1](1, 2) = 12 * y * z;
    gggXYZ[47][1](2, 1) = 12 * y * z;
    gggXYZ[47][1](2, 2) = 6 * y2;
    gggXYZ[47][2](1, 1) = 12 * y * z;
    gggXYZ[47][2](1, 2) = 6 * y2;
    gggXYZ[47][2](2, 1) = 6 * y2;
    gggXYZ[46][0](0, 1) = 6 * y2;
    gggXYZ[46][0](1, 0) = 6 * y2;
    gggXYZ[46][0](1, 1) = 12 * x * y;
    gggXYZ[46][1](0, 0) = 6 * y2;
    gggXYZ[46][1](0, 1) = 12 * x * y;
    gggXYZ[46][1](1, 0) = 12 * x * y;
    gggXYZ[46][1](1, 1) = 6 * x2;
    gggXYZ[45][0](0, 0) = 6 * z2;
    gggXYZ[45][0](0, 2) = 12 * x * z;
    gggXYZ[45][0](2, 0) = 12 * x * z;
    gggXYZ[45][0](2, 2) = 6 * x2;
    gggXYZ[45][2](0, 0) = 12 * x * z;
    gggXYZ[45][2](0, 2) = 6 * x2;
    gggXYZ[45][2](2, 0) = 6 * x2;
    gggXYZ[44][0](0, 0) = 6 * y2;
    gggXYZ[44][0](0, 1) = 12 * x * y;
    gggXYZ[44][0](1, 0) = 12 * x * y;
    gggXYZ[44][0](1, 1) = 6 * x2;
    gggXYZ[44][1](0, 0) = 12 * x * y;
    gggXYZ[44][1](0, 1) = 6 * x2;
    gggXYZ[44][1](1, 0) = 6 * x2;
    gggXYZ[43][1](2, 2) = 12 * z2;
    gggXYZ[43][2](1, 2) = 12 * z2;
    gggXYZ[43][2](2, 1) = 12 * z2;
    gggXYZ[43][2](2, 2) = 24 * y * z;
    gggXYZ[42][0](2, 2) = 12 * z2;
    gggXYZ[42][2](0, 2) = 12 * z2;
    gggXYZ[42][2](2, 0) = 12 * z2;
    gggXYZ[42][2](2, 2) = 24 * x * z;
    gggXYZ[41][1](1, 1) = 24 * y * z;
    gggXYZ[41][1](1, 2) = 12 * y2;
    gggXYZ[41][1](2, 1) = 12 * y2;
    gggXYZ[41][2](1, 1) = 12 * y2;
    gggXYZ[40][0](1, 1) = 12 * y2;
    gggXYZ[40][1](0, 1) = 12 * y2;
    gggXYZ[40][1](1, 0) = 12 * y2;
    gggXYZ[40][1](1, 1) = 24 * x * y;
    gggXYZ[39][0](0, 0) = 24 * x * z;
    gggXYZ[39][0](0, 2) = 12 * x2;
    gggXYZ[39][0](2, 0) = 12 * x2;
    gggXYZ[39][2](0, 0) = 12 * x2;
    gggXYZ[38][0](0, 0) = 24 * x * y;
    gggXYZ[38][0](0, 1) = 12 * x2;
    gggXYZ[38][0](1, 0) = 12 * x2;
    gggXYZ[38][1](0, 0) = 12 * x2;
    gggXYZ[37][2](2, 2) = 60 * z2;
    gggXYZ[36][1](1, 1) = 60 * y2;
    gggXYZ[35][0](0, 0) = 60 * x2;
  case 4:
    gggXYZ[34][0](1, 2) = 2 * z;
    gggXYZ[34][0](2, 1) = 2 * z;
    gggXYZ[34][0](2, 2) = 2 * y;
    gggXYZ[34][1](0, 2) = 2 * z;
    gggXYZ[34][1](2, 0) = 2 * z;
    gggXYZ[34][1](2, 2) = 2 * x;
    gggXYZ[34][2](0, 1) = 2 * z;
    gggXYZ[34][2](0, 2) = 2 * y;
    gggXYZ[34][2](1, 0) = 2 * z;
    gggXYZ[34][2](1, 2) = 2 * x;
    gggXYZ[34][2](2, 0) = 2 * y;
    gggXYZ[34][2](2, 1) = 2 * x;
    gggXYZ[33][0](1, 1) = 2 * z;
    gggXYZ[33][0](1, 2) = 2 * y;
    gggXYZ[33][0](2, 1) = 2 * y;
    gggXYZ[33][1](0, 1) = 2 * z;
    gggXYZ[33][1](0, 2) = 2 * y;
    gggXYZ[33][1](1, 0) = 2 * z;
    gggXYZ[33][1](1, 2) = 2 * x;
    gggXYZ[33][1](2, 0) = 2 * y;
    gggXYZ[33][1](2, 1) = 2 * x;
    gggXYZ[33][2](0, 1) = 2 * y;
    gggXYZ[33][2](1, 0) = 2 * y;
    gggXYZ[33][2](1, 1) = 2 * x;
    gggXYZ[32][0](0, 1) = 2 * z;
    gggXYZ[32][0](0, 2) = 2 * y;
    gggXYZ[32][0](1, 0) = 2 * z;
    gggXYZ[32][0](1, 2) = 2 * x;
    gggXYZ[32][0](2, 0) = 2 * y;
    gggXYZ[32][0](2, 1) = 2 * x;
    gggXYZ[32][1](0, 0) = 2 * z;
    gggXYZ[32][1](0, 2) = 2 * x;
    gggXYZ[32][1](2, 0) = 2 * x;
    gggXYZ[32][2](0, 0) = 2 * y;
    gggXYZ[32][2](0, 1) = 2 * x;
    gggXYZ[32][2](1, 0) = 2 * x;
    gggXYZ[31][1](1, 2) = 4 * z;
    gggXYZ[31][1](2, 1) = 4 * z;
    gggXYZ[31][1](2, 2) = 4 * y;
    gggXYZ[31][2](1, 1) = 4 * z;
    gggXYZ[31][2](1, 2) = 4 * y;
    gggXYZ[31][2](2, 1) = 4 * y;
    gggXYZ[30][0](0, 2) = 4 * z;
    gggXYZ[30][0](2, 0) = 4 * z;
    gggXYZ[30][0](2, 2) = 4 * x;
    gggXYZ[30][2](0, 0) = 4 * z;
    gggXYZ[30][2](0, 2) = 4 * x;
    gggXYZ[30][2](2, 0) = 4 * x;
    gggXYZ[29][0](0, 1) = 4 * y;
    gggXYZ[29][0](1, 0) = 4 * y;
    gggXYZ[29][0](1, 1) = 4 * x;
    gggXYZ[29][1](0, 0) = 4 * y;
    gggXYZ[29][1](0, 1) = 4 * x;
    gggXYZ[29][1](1, 0) = 4 * x;
    gggXYZ[28][1](2, 2) = 6 * z;
    gggXYZ[28][2](1, 2) = 6 * z;
    gggXYZ[28][2](2, 1) = 6 * z;
    gggXYZ[28][2](2, 2) = 6 * y;
    gggXYZ[27][0](2, 2) = 6 * z;
    gggXYZ[27][2](0, 2) = 6 * z;
    gggXYZ[27][2](2, 0) = 6 * z;
    gggXYZ[27][2](2, 2) = 6 * x;
    gggXYZ[26][1](1, 1) = 6 * z;
    gggXYZ[26][1](1, 2) = 6 * y;
    gggXYZ[26][1](2, 1) = 6 * y;
    gggXYZ[26][2](1, 1) = 6 * y;
    gggXYZ[25][0](1, 1) = 6 * y;
    gggXYZ[25][1](0, 1) = 6 * y;
    gggXYZ[25][1](1, 0) = 6 * y;
    gggXYZ[25][1](1, 1) = 6 * x;
    gggXYZ[24][0](0, 0) = 6 * z;
    gggXYZ[24][0](0, 2) = 6 * x;
    gggXYZ[24][0](2, 0) = 6 * x;
    gggXYZ[24][2](0, 0) = 6 * x;
    gggXYZ[23][0](0, 0) = 6 * y;
    gggXYZ[23][0](0, 1) = 6 * x;
    gggXYZ[23][0](1, 0) = 6 * x;
    gggXYZ[23][1](0, 0) = 6 * x;
    gggXYZ[22][2](2, 2) = 24 * z;
    gggXYZ[21][1](1, 1) = 24 * y;
    gggXYZ[20][0](0, 0) = 24 * x;
  case 3:
    gggXYZ[19][0](1, 2) = 1;
    gggXYZ[19][0](2, 1) = 1;
    gggXYZ[19][1](0, 2) = 1;
    gggXYZ[19][1](2, 0) = 1;
    gggXYZ[19][2](0, 1) = 1;
    gggXYZ[19][2](1, 0) = 1;
    gggXYZ[18][1](2, 2) = 2;
    gggXYZ[18][2](1, 2) = 2;
    gggXYZ[18][2](2, 1) = 2;
    gggXYZ[17][0](2, 2) = 2;
    gggXYZ[17][2](0, 2) = 2;
    gggXYZ[17][2](2, 0) = 2;
    gggXYZ[16][1](1, 2) = 2;
    gggXYZ[16][1](2, 1) = 2;
    gggXYZ[16][2](1, 1) = 2;
    gggXYZ[15][0](1, 1) = 2;
    gggXYZ[15][1](0, 1) = 2;
    gggXYZ[15][1](1, 0) = 2;
    gggXYZ[14][0](0, 2) = 2;
    gggXYZ[14][0](2, 0) = 2;
    gggXYZ[14][2](0, 0) = 2;
    gggXYZ[13][0](0, 1) = 2;
    gggXYZ[13][0](1, 0) = 2;
    gggXYZ[13][1](0, 0) = 2;
    gggXYZ[12][2](2, 2) = 6;
    gggXYZ[11][1](1, 1) = 6;
    gggXYZ[10][0](0, 0) = 6;
  case 2:
  case 1:; // empty statement
  case 0:; // empty statement
  }

  for (int i = 0; i < ntot; i++)
  {
    gggXYZ[i][0] *= NormFactor[i];
    gggXYZ[i][1] *= NormFactor[i];
    gggXYZ[i][2] *= NormFactor[i];
  }
}


// generated from read_order.py
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T, Point_t, Tensor_t, GGG_t>::getABC(int n, int& a, int& b, int& c)
{
  // following Gamess notation
  switch (n)
  {
  // S
  case 0: // S
    a = 0;
    b = 0;
    c = 0;
    break;
  // P
  case 1: // X
    a = 1;
    b = 0;
    c = 0;
    break;
  case 2: // Y
    a = 0;
    b = 1;
    c = 0;
    break;
  case 3: // Z
    a = 0;
    b = 0;
    c = 1;
    break;
  // D
  case 4: // XX
    a = 2;
    b = 0;
    c = 0;
    break;
  case 5: // YY
    a = 0;
    b = 2;
    c = 0;
    break;
  case 6: // ZZ
    a = 0;
    b = 0;
    c = 2;
    break;
  case 7: // XY
    a = 1;
    b = 1;
    c = 0;
    break;
  case 8: // XZ
    a = 1;
    b = 0;
    c = 1;
    break;
  case 9: // YZ
    a = 0;
    b = 1;
    c = 1;
    break;
  // F
  case 10: // XXX
    a = 3;
    b = 0;
    c = 0;
    break;
  case 11: // YYY
    a = 0;
    b = 3;
    c = 0;
    break;
  case 12: // ZZZ
    a = 0;
    b = 0;
    c = 3;
    break;
  case 13: // XXY
    a = 2;
    b = 1;
    c = 0;
    break;
  case 14: // XXZ
    a = 2;
    b = 0;
    c = 1;
    break;
  case 15: // YYX
    a = 1;
    b = 2;
    c = 0;
    break;
  case 16: // YYZ
    a = 0;
    b = 2;
    c = 1;
    break;
  case 17: // ZZX
    a = 1;
    b = 0;
    c = 2;
    break;
  case 18: // ZZY
    a = 0;
    b = 1;
    c = 2;
    break;
  case 19: // XYZ
    a = 1;
    b = 1;
    c = 1;
    break;
  // G
  case 20: // XXXX
    a = 4;
    b = 0;
    c = 0;
    break;
  case 21: // YYYY
    a = 0;
    b = 4;
    c = 0;
    break;
  case 22: // ZZZZ
    a = 0;
    b = 0;
    c = 4;
    break;
  case 23: // XXXY
    a = 3;
    b = 1;
    c = 0;
    break;
  case 24: // XXXZ
    a = 3;
    b = 0;
    c = 1;
    break;
  case 25: // YYYX
    a = 1;
    b = 3;
    c = 0;
    break;
  case 26: // YYYZ
    a = 0;
    b = 3;
    c = 1;
    break;
  case 27: // ZZZX
    a = 1;
    b = 0;
    c = 3;
    break;
  case 28: // ZZZY
    a = 0;
    b = 1;
    c = 3;
    break;
  case 29: // XXYY
    a = 2;
    b = 2;
    c = 0;
    break;
  case 30: // XXZZ
    a = 2;
    b = 0;
    c = 2;
    break;
  case 31: // YYZZ
    a = 0;
    b = 2;
    c = 2;
    break;
  case 32: // XXYZ
    a = 2;
    b = 1;
    c = 1;
    break;
  case 33: // YYXZ
    a = 1;
    b = 2;
    c = 1;
    break;
  case 34: // ZZXY
    a = 1;
    b = 1;
    c = 2;
    break;
  // H
  case 35: // XXXXX
    a = 5;
    b = 0;
    c = 0;
    break;
  case 36: // YYYYY
    a = 0;
    b = 5;
    c = 0;
    break;
  case 37: // ZZZZZ
    a = 0;
    b = 0;
    c = 5;
    break;
  case 38: // XXXXY
    a = 4;
    b = 1;
    c = 0;
    break;
  case 39: // XXXXZ
    a = 4;
    b = 0;
    c = 1;
    break;
  case 40: // YYYYX
    a = 1;
    b = 4;
    c = 0;
    break;
  case 41: // YYYYZ
    a = 0;
    b = 4;
    c = 1;
    break;
  case 42: // ZZZZX
    a = 1;
    b = 0;
    c = 4;
    break;
  case 43: // ZZZZY
    a = 0;
    b = 1;
    c = 4;
    break;
  case 44: // XXXYY
    a = 3;
    b = 2;
    c = 0;
    break;
  case 45: // XXXZZ
    a = 3;
    b = 0;
    c = 2;
    break;
  case 46: // YYYXX
    a = 2;
    b = 3;
    c = 0;
    break;
  case 47: // YYYZZ
    a = 0;
    b = 3;
    c = 2;
    break;
  case 48: // ZZZXX
    a = 2;
    b = 0;
    c = 3;
    break;
  case 49: // ZZZYY
    a = 0;
    b = 2;
    c = 3;
    break;
  case 50: // XXXYZ
    a = 3;
    b = 1;
    c = 1;
    break;
  case 51: // YYYXZ
    a = 1;
    b = 3;
    c = 1;
    break;
  case 52: // ZZZXY
    a = 1;
    b = 1;
    c = 3;
    break;
  case 53: // XXYYZ
    a = 2;
    b = 2;
    c = 1;
    break;
  case 54: // XXZZY
    a = 2;
    b = 1;
    c = 2;
    break;
  case 55: // YYZZX
    a = 1;
    b = 2;
    c = 2;
    break;
  // I
  case 56: // X6
    a = 6;
    b = 0;
    c = 0;
    break;
  case 57: // Y6
    a = 0;
    b = 6;
    c = 0;
    break;
  case 58: // Z6
    a = 0;
    b = 0;
    c = 6;
    break;
  case 59: // X5Y
    a = 5;
    b = 1;
    c = 0;
    break;
  case 60: // X5Z
    a = 5;
    b = 0;
    c = 1;
    break;
  case 61: // Y5X
    a = 1;
    b = 5;
    c = 0;
    break;
  case 62: // Y5Z
    a = 0;
    b = 5;
    c = 1;
    break;
  case 63: // Z5X
    a = 1;
    b = 0;
    c = 5;
    break;
  case 64: // Z5Y
    a = 0;
    b = 1;
    c = 5;
    break;
  case 65: // X4Y2
    a = 4;
    b = 2;
    c = 0;
    break;
  case 66: // X4Z2
    a = 4;
    b = 0;
    c = 2;
    break;
  case 67: // Y4X2
    a = 2;
    b = 4;
    c = 0;
    break;
  case 68: // Y4Z2
    a = 0;
    b = 4;
    c = 2;
    break;
  case 69: // Z4X2
    a = 2;
    b = 0;
    c = 4;
    break;
  case 70: // Z4Y2
    a = 0;
    b = 2;
    c = 4;
    break;
  case 71: // X4YZ
    a = 4;
    b = 1;
    c = 1;
    break;
  case 72: // Y4XZ
    a = 1;
    b = 4;
    c = 1;
    break;
  case 73: // Z4XY
    a = 1;
    b = 1;
    c = 4;
    break;
  case 74: // X3Y3
    a = 3;
    b = 3;
    c = 0;
    break;
  case 75: // X3Z3
    a = 3;
    b = 0;
    c = 3;
    break;
  case 76: // Y3Z3
    a = 0;
    b = 3;
    c = 3;
    break;
  case 77: // X3Y2Z
    a = 3;
    b = 2;
    c = 1;
    break;
  case 78: // X3Z2Y
    a = 3;
    b = 1;
    c = 2;
    break;
  case 79: // Y3X2Z
    a = 2;
    b = 3;
    c = 1;
    break;
  case 80: // Y3Z2X
    a = 1;
    b = 3;
    c = 2;
    break;
  case 81: // Z3X2Y
    a = 2;
    b = 1;
    c = 3;
    break;
  case 82: // Z3Y2X
    a = 1;
    b = 2;
    c = 3;
    break;
  case 83: // X2Y2Z2
    a = 2;
    b = 2;
    c = 2;
    break;

  default:
    std::cerr << "CartesianTensor::getABC() - Incorrect index." << std::endl;
    APP_ABORT("");
    break;
  }
}

#endif
