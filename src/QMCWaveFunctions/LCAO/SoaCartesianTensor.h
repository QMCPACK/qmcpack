//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/Numerics/codegen/gen_cartesian_tensor.py and SoaCartesianTensor.h.in

 Edit SoaCartesianTensor.h.in, rerun gen_cartesian_tensor.py, and copy the generated file here.
*/


#ifndef QMCPLUSPLUS_SOA_CARTESIAN_TENSOR_H
#define QMCPLUSPLUS_SOA_CARTESIAN_TENSOR_H

#include <stdexcept>
#include "OhmmsSoA/VectorSoaContainer.h"

namespace qmcplusplus
{
/** CartesianTensor according to Gamess order
 * @tparam T, value_type, e.g. double
 *
 * Original implementation Numerics/CartesianTensor.h
 * Modified to use SoA for cXYZ and used by SoaAtomicBasisSet
 * Array ordered as [S,X,Y,Z,XX,YY,ZZ,XY,XZ,YZ,...]
 *    (following Gamess order)
 */
template<class T>
struct SoaCartesianTensor
{
  using value_type = T;
  using ggg_type   = TinyVector<Tensor<T, 3>, 3>;

  ///maximum angular momentum
  size_t Lmax;
  ///normalization factor
  aligned_vector<T> NormFactor;
  ///composite V,Gx,Gy,Gz,[L | H00, H01, H02, H11, H12, H12]
  //   {GH000, GH001, GH002, GH011, GH012, GH022, GH111, GH112, GH122, GH222
  VectorSoaContainer<T, 20> cXYZ;

  /** constructor
   * @param l_max maximum angular momentum
   *
   * Evaluate all the constants and prefactors.
  */
  explicit SoaCartesianTensor(const int l_max, bool addsign = false);

  ///compute Ylm
  void evaluate_bare(T x, T y, T z, T* XYZ) const;

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* XYZ) const
  {
    evaluate_bare(x, y, z, XYZ);
    for (size_t i = 0, nl = cXYZ.size(); i < nl; i++)
      XYZ[i] *= NormFactor[i];
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* XYZ)
  {
    evaluate_bare(x, y, z, XYZ);
    for (size_t i = 0, nl = cXYZ.size(); i < nl; i++)
      XYZ[i] *= NormFactor[i];
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z) { evaluateV(x, y, z, cXYZ.data(0)); }

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGL(T x, T y, T z);

  void evaluateVGH(T x, T y, T z);

  void evaluateVGHGH(T x, T y, T z);
  ///returns dummy: this is not used
  inline int index(int l, int m) const { return (l * (l + 1)) + m; }

  /** return the starting address of the component
   *
   * component=0(V), 1(dx), 2(dy), 3(dz), 4(Lap)
   * See comment at VectorSoAContainer<T,20> cXYZ declaration.
   */
  inline const T* operator[](size_t component) const { return cXYZ.data(component); }

  inline size_t size() const { return cXYZ.size(); }

  inline int lmax() const { return Lmax; }

  inline void getABC(int n, int& a, int& b, int& c);

  int DFactorial(int num) { return (num < 2) ? 1 : num * DFactorial(num - 2); }
};

template<class T>
SoaCartesianTensor<T>::SoaCartesianTensor(const int l_max, bool addsign) : Lmax(l_max)
{
  if (Lmax < 0 || Lmax > 6)
    throw std::runtime_error("CartesianTensor can't handle Lmax > 6 or Lmax < 0.\n");

  int ntot = 0;
  for (int i = 0; i <= Lmax; i++)
    ntot += (i + 1) * (i + 2) / 2;
  cXYZ.resize(ntot);
  NormFactor.resize(ntot, 1);
  int p = 0;
  int a = 0, b = 0, c = 0;
  const double pi = 4.0 * std::atan(1.0);
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
      double L = static_cast<T>(l);
      double NormL =
          std::pow(2, L + 1) * std::sqrt(2.0 / static_cast<double>(DFactorial(2 * l + 1))) * std::pow(2.0 / pi, 0.25);
      NormFactor[p++] = static_cast<T>(
          std::pow(2.0 / pi, 0.75) * std::pow(4.0, 0.5 * (a + b + c)) *
          std::sqrt(1.0 /
                    static_cast<double>((DFactorial(2 * a - 1) * DFactorial(2 * b - 1) * DFactorial(2 * c - 1)))) /
          NormL);
    }
  }
}


template<class T>
void SoaCartesianTensor<T>::evaluate_bare(T x, T y, T z, T* restrict XYZ) const
{
  const T x2 = x * x, y2 = y * y, z2 = z * z;
  const T x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  const T x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  const T x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;
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
}


template<class T>
void SoaCartesianTensor<T>::evaluateVGL(T x, T y, T z)
{
  constexpr T czero(0);
  cXYZ = czero;

  const T x2 = x * x, y2 = y * y, z2 = z * z;
  const T x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  const T x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  const T x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;
  T* restrict XYZ = cXYZ.data(0);
  T* restrict gr0 = cXYZ.data(1);
  T* restrict gr1 = cXYZ.data(2);
  T* restrict gr2 = cXYZ.data(3);
  T* restrict lap = cXYZ.data(4);

  switch (Lmax)
  {
  case 6:
    XYZ[83] = x2 * y2 * z2; // X2Y2Z2
    gr0[83] = 2 * x * y2 * z2;
    gr1[83] = 2 * x2 * y * z2;
    gr2[83] = 2 * x2 * y2 * z;
    lap[83] = 2 * x2 * y2 + 2 * x2 * z2 + 2 * y2 * z2;
    XYZ[82] = x * y2 * z3; // Z3Y2X
    gr0[82] = y2 * z3;
    gr1[82] = 2 * x * y * z3;
    gr2[82] = 3 * x * y2 * z2;
    lap[82] = 6 * x * y2 * z + 2 * x * z3;
    XYZ[81] = x2 * y * z3; // Z3X2Y
    gr0[81] = 2 * x * y * z3;
    gr1[81] = x2 * z3;
    gr2[81] = 3 * x2 * y * z2;
    lap[81] = 6 * x2 * y * z + 2 * y * z3;
    XYZ[80] = x * y3 * z2; // Y3Z2X
    gr0[80] = y3 * z2;
    gr1[80] = 3 * x * y2 * z2;
    gr2[80] = 2 * x * y3 * z;
    lap[80] = 6 * x * y * z2 + 2 * x * y3;
    XYZ[79] = x2 * y3 * z; // Y3X2Z
    gr0[79] = 2 * x * y3 * z;
    gr1[79] = 3 * x2 * y2 * z;
    gr2[79] = x2 * y3;
    lap[79] = 6 * x2 * y * z + 2 * y3 * z;
    XYZ[78] = x3 * y * z2; // X3Z2Y
    gr0[78] = 3 * x2 * y * z2;
    gr1[78] = x3 * z2;
    gr2[78] = 2 * x3 * y * z;
    lap[78] = 6 * x * y * z2 + 2 * x3 * y;
    XYZ[77] = x3 * y2 * z; // X3Y2Z
    gr0[77] = 3 * x2 * y2 * z;
    gr1[77] = 2 * x3 * y * z;
    gr2[77] = x3 * y2;
    lap[77] = 6 * x * y2 * z + 2 * x3 * z;
    XYZ[76] = y3 * z3; // Y3Z3
    gr1[76] = 3 * y2 * z3;
    gr2[76] = 3 * y3 * z2;
    lap[76] = 6 * y * z3 + 6 * y3 * z;
    XYZ[75] = x3 * z3; // X3Z3
    gr0[75] = 3 * x2 * z3;
    gr2[75] = 3 * x3 * z2;
    lap[75] = 6 * x * z3 + 6 * x3 * z;
    XYZ[74] = x3 * y3; // X3Y3
    gr0[74] = 3 * x2 * y3;
    gr1[74] = 3 * x3 * y2;
    lap[74] = 6 * x * y3 + 6 * x3 * y;
    XYZ[73] = x * y * z4; // Z4XY
    gr0[73] = y * z4;
    gr1[73] = x * z4;
    gr2[73] = 4 * x * y * z3;
    lap[73] = 12 * x * y * z2;
    XYZ[72] = x * y4 * z; // Y4XZ
    gr0[72] = y4 * z;
    gr1[72] = 4 * x * y3 * z;
    gr2[72] = x * y4;
    lap[72] = 12 * x * y2 * z;
    XYZ[71] = x4 * y * z; // X4YZ
    gr0[71] = 4 * x3 * y * z;
    gr1[71] = x4 * z;
    gr2[71] = x4 * y;
    lap[71] = 12 * x2 * y * z;
    XYZ[70] = y2 * z4; // Z4Y2
    gr1[70] = 2 * y * z4;
    gr2[70] = 4 * y2 * z3;
    lap[70] = 12 * y2 * z2 + 2 * z4;
    XYZ[69] = x2 * z4; // Z4X2
    gr0[69] = 2 * x * z4;
    gr2[69] = 4 * x2 * z3;
    lap[69] = 12 * x2 * z2 + 2 * z4;
    XYZ[68] = y4 * z2; // Y4Z2
    gr1[68] = 4 * y3 * z2;
    gr2[68] = 2 * y4 * z;
    lap[68] = 12 * y2 * z2 + 2 * y4;
    XYZ[67] = x2 * y4; // Y4X2
    gr0[67] = 2 * x * y4;
    gr1[67] = 4 * x2 * y3;
    lap[67] = 12 * x2 * y2 + 2 * y4;
    XYZ[66] = x4 * z2; // X4Z2
    gr0[66] = 4 * x3 * z2;
    gr2[66] = 2 * x4 * z;
    lap[66] = 12 * x2 * z2 + 2 * x4;
    XYZ[65] = x4 * y2; // X4Y2
    gr0[65] = 4 * x3 * y2;
    gr1[65] = 2 * x4 * y;
    lap[65] = 12 * x2 * y2 + 2 * x4;
    XYZ[64] = y * z * z4; // Z5Y
    gr1[64] = z * z4;
    gr2[64] = 5 * y * z4;
    lap[64] = 20 * y * z3;
    XYZ[63] = x * z * z4; // Z5X
    gr0[63] = z * z4;
    gr2[63] = 5 * x * z4;
    lap[63] = 20 * x * z3;
    XYZ[62] = y * y4 * z; // Y5Z
    gr1[62] = 5 * y4 * z;
    gr2[62] = y * y4;
    lap[62] = 20 * y3 * z;
    XYZ[61] = x * y * y4; // Y5X
    gr0[61] = y * y4;
    gr1[61] = 5 * x * y4;
    lap[61] = 20 * x * y3;
    XYZ[60] = x * x4 * z; // X5Z
    gr0[60] = 5 * x4 * z;
    gr2[60] = x * x4;
    lap[60] = 20 * x3 * z;
    XYZ[59] = x * x4 * y; // X5Y
    gr0[59] = 5 * x4 * y;
    gr1[59] = x * x4;
    lap[59] = 20 * x3 * y;
    XYZ[58] = z * z5; // Z6
    gr2[58] = 6 * z * z4;
    lap[58] = 30 * z4;
    XYZ[57] = y * y5; // Y6
    gr1[57] = 6 * y * y4;
    lap[57] = 30 * y4;
    XYZ[56] = x * x5; // X6
    gr0[56] = 6 * x * x4;
    lap[56] = 30 * x4;
  case 5:
    XYZ[55] = x * y2 * z2; // YYZZX
    gr0[55] = y2 * z2;
    gr1[55] = 2 * x * y * z2;
    gr2[55] = 2 * x * y2 * z;
    lap[55] = 2 * x * y2 + 2 * x * z2;
    XYZ[54] = x2 * y * z2; // XXZZY
    gr0[54] = 2 * x * y * z2;
    gr1[54] = x2 * z2;
    gr2[54] = 2 * x2 * y * z;
    lap[54] = 2 * x2 * y + 2 * y * z2;
    XYZ[53] = x2 * y2 * z; // XXYYZ
    gr0[53] = 2 * x * y2 * z;
    gr1[53] = 2 * x2 * y * z;
    gr2[53] = x2 * y2;
    lap[53] = 2 * x2 * z + 2 * y2 * z;
    XYZ[52] = x * y * z3; // ZZZXY
    gr0[52] = y * z3;
    gr1[52] = x * z3;
    gr2[52] = 3 * x * y * z2;
    lap[52] = 6 * x * y * z;
    XYZ[51] = x * y3 * z; // YYYXZ
    gr0[51] = y3 * z;
    gr1[51] = 3 * x * y2 * z;
    gr2[51] = x * y3;
    lap[51] = 6 * x * y * z;
    XYZ[50] = x3 * y * z; // XXXYZ
    gr0[50] = 3 * x2 * y * z;
    gr1[50] = x3 * z;
    gr2[50] = x3 * y;
    lap[50] = 6 * x * y * z;
    XYZ[49] = y2 * z3; // ZZZYY
    gr1[49] = 2 * y * z3;
    gr2[49] = 3 * y2 * z2;
    lap[49] = 6 * y2 * z + 2 * z3;
    XYZ[48] = x2 * z3; // ZZZXX
    gr0[48] = 2 * x * z3;
    gr2[48] = 3 * x2 * z2;
    lap[48] = 6 * x2 * z + 2 * z3;
    XYZ[47] = y3 * z2; // YYYZZ
    gr1[47] = 3 * y2 * z2;
    gr2[47] = 2 * y3 * z;
    lap[47] = 6 * y * z2 + 2 * y3;
    XYZ[46] = x2 * y3; // YYYXX
    gr0[46] = 2 * x * y3;
    gr1[46] = 3 * x2 * y2;
    lap[46] = 6 * x2 * y + 2 * y3;
    XYZ[45] = x3 * z2; // XXXZZ
    gr0[45] = 3 * x2 * z2;
    gr2[45] = 2 * x3 * z;
    lap[45] = 6 * x * z2 + 2 * x3;
    XYZ[44] = x3 * y2; // XXXYY
    gr0[44] = 3 * x2 * y2;
    gr1[44] = 2 * x3 * y;
    lap[44] = 6 * x * y2 + 2 * x3;
    XYZ[43] = y * z4; // ZZZZY
    gr1[43] = z4;
    gr2[43] = 4 * y * z3;
    lap[43] = 12 * y * z2;
    XYZ[42] = x * z4; // ZZZZX
    gr0[42] = z4;
    gr2[42] = 4 * x * z3;
    lap[42] = 12 * x * z2;
    XYZ[41] = y4 * z; // YYYYZ
    gr1[41] = 4 * y3 * z;
    gr2[41] = y4;
    lap[41] = 12 * y2 * z;
    XYZ[40] = x * y4; // YYYYX
    gr0[40] = y4;
    gr1[40] = 4 * x * y3;
    lap[40] = 12 * x * y2;
    XYZ[39] = x4 * z; // XXXXZ
    gr0[39] = 4 * x3 * z;
    gr2[39] = x4;
    lap[39] = 12 * x2 * z;
    XYZ[38] = x4 * y; // XXXXY
    gr0[38] = 4 * x3 * y;
    gr1[38] = x4;
    lap[38] = 12 * x2 * y;
    XYZ[37] = z * z4; // ZZZZZ
    gr2[37] = 5 * z4;
    lap[37] = 20 * z3;
    XYZ[36] = y * y4; // YYYYY
    gr1[36] = 5 * y4;
    lap[36] = 20 * y3;
    XYZ[35] = x * x4; // XXXXX
    gr0[35] = 5 * x4;
    lap[35] = 20 * x3;
  case 4:
    XYZ[34] = x * y * z2; // ZZXY
    gr0[34] = y * z2;
    gr1[34] = x * z2;
    gr2[34] = 2 * x * y * z;
    lap[34] = 2 * x * y;
    XYZ[33] = x * y2 * z; // YYXZ
    gr0[33] = y2 * z;
    gr1[33] = 2 * x * y * z;
    gr2[33] = x * y2;
    lap[33] = 2 * x * z;
    XYZ[32] = x2 * y * z; // XXYZ
    gr0[32] = 2 * x * y * z;
    gr1[32] = x2 * z;
    gr2[32] = x2 * y;
    lap[32] = 2 * y * z;
    XYZ[31] = y2 * z2; // YYZZ
    gr1[31] = 2 * y * z2;
    gr2[31] = 2 * y2 * z;
    lap[31] = 2 * y2 + 2 * z2;
    XYZ[30] = x2 * z2; // XXZZ
    gr0[30] = 2 * x * z2;
    gr2[30] = 2 * x2 * z;
    lap[30] = 2 * x2 + 2 * z2;
    XYZ[29] = x2 * y2; // XXYY
    gr0[29] = 2 * x * y2;
    gr1[29] = 2 * x2 * y;
    lap[29] = 2 * x2 + 2 * y2;
    XYZ[28] = y * z3; // ZZZY
    gr1[28] = z3;
    gr2[28] = 3 * y * z2;
    lap[28] = 6 * y * z;
    XYZ[27] = x * z3; // ZZZX
    gr0[27] = z3;
    gr2[27] = 3 * x * z2;
    lap[27] = 6 * x * z;
    XYZ[26] = y3 * z; // YYYZ
    gr1[26] = 3 * y2 * z;
    gr2[26] = y3;
    lap[26] = 6 * y * z;
    XYZ[25] = x * y3; // YYYX
    gr0[25] = y3;
    gr1[25] = 3 * x * y2;
    lap[25] = 6 * x * y;
    XYZ[24] = x3 * z; // XXXZ
    gr0[24] = 3 * x2 * z;
    gr2[24] = x3;
    lap[24] = 6 * x * z;
    XYZ[23] = x3 * y; // XXXY
    gr0[23] = 3 * x2 * y;
    gr1[23] = x3;
    lap[23] = 6 * x * y;
    XYZ[22] = z4; // ZZZZ
    gr2[22] = 4 * z3;
    lap[22] = 12 * z2;
    XYZ[21] = y4; // YYYY
    gr1[21] = 4 * y3;
    lap[21] = 12 * y2;
    XYZ[20] = x4; // XXXX
    gr0[20] = 4 * x3;
    lap[20] = 12 * x2;
  case 3:
    XYZ[19] = x * y * z; // XYZ
    gr0[19] = y * z;
    gr1[19] = x * z;
    gr2[19] = x * y;
    XYZ[18] = y * z2; // ZZY
    gr1[18] = z2;
    gr2[18] = 2 * y * z;
    lap[18] = 2 * y;
    XYZ[17] = x * z2; // ZZX
    gr0[17] = z2;
    gr2[17] = 2 * x * z;
    lap[17] = 2 * x;
    XYZ[16] = y2 * z; // YYZ
    gr1[16] = 2 * y * z;
    gr2[16] = y2;
    lap[16] = 2 * z;
    XYZ[15] = x * y2; // YYX
    gr0[15] = y2;
    gr1[15] = 2 * x * y;
    lap[15] = 2 * x;
    XYZ[14] = x2 * z; // XXZ
    gr0[14] = 2 * x * z;
    gr2[14] = x2;
    lap[14] = 2 * z;
    XYZ[13] = x2 * y; // XXY
    gr0[13] = 2 * x * y;
    gr1[13] = x2;
    lap[13] = 2 * y;
    XYZ[12] = z3; // ZZZ
    gr2[12] = 3 * z2;
    lap[12] = 6 * z;
    XYZ[11] = y3; // YYY
    gr1[11] = 3 * y2;
    lap[11] = 6 * y;
    XYZ[10] = x3; // XXX
    gr0[10] = 3 * x2;
    lap[10] = 6 * x;
  case 2:
    XYZ[9] = y * z; // YZ
    gr1[9] = z;
    gr2[9] = y;
    XYZ[8] = x * z; // XZ
    gr0[8] = z;
    gr2[8] = x;
    XYZ[7] = x * y; // XY
    gr0[7] = y;
    gr1[7] = x;
    XYZ[6] = z2; // ZZ
    gr2[6] = 2 * z;
    lap[6] = 2;
    XYZ[5] = y2; // YY
    gr1[5] = 2 * y;
    lap[5] = 2;
    XYZ[4] = x2; // XX
    gr0[4] = 2 * x;
    lap[4] = 2;
  case 1:
    XYZ[3] = z; // Z
    gr2[3] = 1;
    XYZ[2] = y; // Y
    gr1[2] = 1;
    XYZ[1] = x; // X
    gr0[1] = 1;
  case 0:
    XYZ[0] = 1; // S
  }

  const size_t ntot = NormFactor.size();
  for (size_t i = 0; i < ntot; i++)
  {
    XYZ[i] *= NormFactor[i];
    gr0[i] *= NormFactor[i];
    gr1[i] *= NormFactor[i];
    gr2[i] *= NormFactor[i];
    lap[i] *= NormFactor[i];
  }
}


template<class T>
void SoaCartesianTensor<T>::evaluateVGH(T x, T y, T z)
{
  constexpr T czero(0);
  cXYZ = czero;

  const T x2 = x * x, y2 = y * y, z2 = z * z;
  const T x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  const T x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  const T x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;

  T* restrict XYZ = cXYZ.data(0);
  T* restrict gr0 = cXYZ.data(1);
  T* restrict gr1 = cXYZ.data(2);
  T* restrict gr2 = cXYZ.data(3);
  T* restrict h00 = cXYZ.data(4);
  T* restrict h01 = cXYZ.data(5);
  T* restrict h02 = cXYZ.data(6);
  T* restrict h11 = cXYZ.data(7);
  T* restrict h12 = cXYZ.data(8);
  T* restrict h22 = cXYZ.data(9);


  switch (Lmax)
  {
  case 6:
    XYZ[83] = x2 * y2 * z2; // X2Y2Z2
    gr0[83] = 2 * x * y2 * z2;
    gr1[83] = 2 * x2 * y * z2;
    gr2[83] = 2 * x2 * y2 * z;
    h00[83] = 2 * y2 * z2;
    h01[83] = 4 * x * y * z2;
    h02[83] = 4 * x * y2 * z;
    h11[83] = 2 * x2 * z2;
    h12[83] = 4 * x2 * y * z;
    h22[83] = 2 * x2 * y2;
    XYZ[82] = x * y2 * z3; // Z3Y2X
    gr0[82] = y2 * z3;
    gr1[82] = 2 * x * y * z3;
    gr2[82] = 3 * x * y2 * z2;
    h01[82] = 2 * y * z3;
    h02[82] = 3 * y2 * z2;
    h11[82] = 2 * x * z3;
    h12[82] = 6 * x * y * z2;
    h22[82] = 6 * x * y2 * z;
    XYZ[81] = x2 * y * z3; // Z3X2Y
    gr0[81] = 2 * x * y * z3;
    gr1[81] = x2 * z3;
    gr2[81] = 3 * x2 * y * z2;
    h00[81] = 2 * y * z3;
    h01[81] = 2 * x * z3;
    h02[81] = 6 * x * y * z2;
    h12[81] = 3 * x2 * z2;
    h22[81] = 6 * x2 * y * z;
    XYZ[80] = x * y3 * z2; // Y3Z2X
    gr0[80] = y3 * z2;
    gr1[80] = 3 * x * y2 * z2;
    gr2[80] = 2 * x * y3 * z;
    h01[80] = 3 * y2 * z2;
    h02[80] = 2 * y3 * z;
    h11[80] = 6 * x * y * z2;
    h12[80] = 6 * x * y2 * z;
    h22[80] = 2 * x * y3;
    XYZ[79] = x2 * y3 * z; // Y3X2Z
    gr0[79] = 2 * x * y3 * z;
    gr1[79] = 3 * x2 * y2 * z;
    gr2[79] = x2 * y3;
    h00[79] = 2 * y3 * z;
    h01[79] = 6 * x * y2 * z;
    h02[79] = 2 * x * y3;
    h11[79] = 6 * x2 * y * z;
    h12[79] = 3 * x2 * y2;
    XYZ[78] = x3 * y * z2; // X3Z2Y
    gr0[78] = 3 * x2 * y * z2;
    gr1[78] = x3 * z2;
    gr2[78] = 2 * x3 * y * z;
    h00[78] = 6 * x * y * z2;
    h01[78] = 3 * x2 * z2;
    h02[78] = 6 * x2 * y * z;
    h12[78] = 2 * x3 * z;
    h22[78] = 2 * x3 * y;
    XYZ[77] = x3 * y2 * z; // X3Y2Z
    gr0[77] = 3 * x2 * y2 * z;
    gr1[77] = 2 * x3 * y * z;
    gr2[77] = x3 * y2;
    h00[77] = 6 * x * y2 * z;
    h01[77] = 6 * x2 * y * z;
    h02[77] = 3 * x2 * y2;
    h11[77] = 2 * x3 * z;
    h12[77] = 2 * x3 * y;
    XYZ[76] = y3 * z3; // Y3Z3
    gr1[76] = 3 * y2 * z3;
    gr2[76] = 3 * y3 * z2;
    h11[76] = 6 * y * z3;
    h12[76] = 9 * y2 * z2;
    h22[76] = 6 * y3 * z;
    XYZ[75] = x3 * z3; // X3Z3
    gr0[75] = 3 * x2 * z3;
    gr2[75] = 3 * x3 * z2;
    h00[75] = 6 * x * z3;
    h02[75] = 9 * x2 * z2;
    h22[75] = 6 * x3 * z;
    XYZ[74] = x3 * y3; // X3Y3
    gr0[74] = 3 * x2 * y3;
    gr1[74] = 3 * x3 * y2;
    h00[74] = 6 * x * y3;
    h01[74] = 9 * x2 * y2;
    h11[74] = 6 * x3 * y;
    XYZ[73] = x * y * z4; // Z4XY
    gr0[73] = y * z4;
    gr1[73] = x * z4;
    gr2[73] = 4 * x * y * z3;
    h01[73] = z4;
    h02[73] = 4 * y * z3;
    h12[73] = 4 * x * z3;
    h22[73] = 12 * x * y * z2;
    XYZ[72] = x * y4 * z; // Y4XZ
    gr0[72] = y4 * z;
    gr1[72] = 4 * x * y3 * z;
    gr2[72] = x * y4;
    h01[72] = 4 * y3 * z;
    h02[72] = y4;
    h11[72] = 12 * x * y2 * z;
    h12[72] = 4 * x * y3;
    XYZ[71] = x4 * y * z; // X4YZ
    gr0[71] = 4 * x3 * y * z;
    gr1[71] = x4 * z;
    gr2[71] = x4 * y;
    h00[71] = 12 * x2 * y * z;
    h01[71] = 4 * x3 * z;
    h02[71] = 4 * x3 * y;
    h12[71] = x4;
    XYZ[70] = y2 * z4; // Z4Y2
    gr1[70] = 2 * y * z4;
    gr2[70] = 4 * y2 * z3;
    h11[70] = 2 * z4;
    h12[70] = 8 * y * z3;
    h22[70] = 12 * y2 * z2;
    XYZ[69] = x2 * z4; // Z4X2
    gr0[69] = 2 * x * z4;
    gr2[69] = 4 * x2 * z3;
    h00[69] = 2 * z4;
    h02[69] = 8 * x * z3;
    h22[69] = 12 * x2 * z2;
    XYZ[68] = y4 * z2; // Y4Z2
    gr1[68] = 4 * y3 * z2;
    gr2[68] = 2 * y4 * z;
    h11[68] = 12 * y2 * z2;
    h12[68] = 8 * y3 * z;
    h22[68] = 2 * y4;
    XYZ[67] = x2 * y4; // Y4X2
    gr0[67] = 2 * x * y4;
    gr1[67] = 4 * x2 * y3;
    h00[67] = 2 * y4;
    h01[67] = 8 * x * y3;
    h11[67] = 12 * x2 * y2;
    XYZ[66] = x4 * z2; // X4Z2
    gr0[66] = 4 * x3 * z2;
    gr2[66] = 2 * x4 * z;
    h00[66] = 12 * x2 * z2;
    h02[66] = 8 * x3 * z;
    h22[66] = 2 * x4;
    XYZ[65] = x4 * y2; // X4Y2
    gr0[65] = 4 * x3 * y2;
    gr1[65] = 2 * x4 * y;
    h00[65] = 12 * x2 * y2;
    h01[65] = 8 * x3 * y;
    h11[65] = 2 * x4;
    XYZ[64] = y * z * z4; // Z5Y
    gr1[64] = z * z4;
    gr2[64] = 5 * y * z4;
    h12[64] = 5 * z4;
    h22[64] = 20 * y * z3;
    XYZ[63] = x * z * z4; // Z5X
    gr0[63] = z * z4;
    gr2[63] = 5 * x * z4;
    h02[63] = 5 * z4;
    h22[63] = 20 * x * z3;
    XYZ[62] = y * y4 * z; // Y5Z
    gr1[62] = 5 * y4 * z;
    gr2[62] = y * y4;
    h11[62] = 20 * y3 * z;
    h12[62] = 5 * y4;
    XYZ[61] = x * y * y4; // Y5X
    gr0[61] = y * y4;
    gr1[61] = 5 * x * y4;
    h01[61] = 5 * y4;
    h11[61] = 20 * x * y3;
    XYZ[60] = x * x4 * z; // X5Z
    gr0[60] = 5 * x4 * z;
    gr2[60] = x * x4;
    h00[60] = 20 * x3 * z;
    h02[60] = 5 * x4;
    XYZ[59] = x * x4 * y; // X5Y
    gr0[59] = 5 * x4 * y;
    gr1[59] = x * x4;
    h00[59] = 20 * x3 * y;
    h01[59] = 5 * x4;
    XYZ[58] = z * z5; // Z6
    gr2[58] = 6 * z * z4;
    h22[58] = 30 * z4;
    XYZ[57] = y * y5; // Y6
    gr1[57] = 6 * y * y4;
    h11[57] = 30 * y4;
    XYZ[56] = x * x5; // X6
    gr0[56] = 6 * x * x4;
    h00[56] = 30 * x4;
  case 5:
    XYZ[55] = x * y2 * z2; // YYZZX
    gr0[55] = y2 * z2;
    gr1[55] = 2 * x * y * z2;
    gr2[55] = 2 * x * y2 * z;
    h01[55] = 2 * y * z2;
    h02[55] = 2 * y2 * z;
    h11[55] = 2 * x * z2;
    h12[55] = 4 * x * y * z;
    h22[55] = 2 * x * y2;
    XYZ[54] = x2 * y * z2; // XXZZY
    gr0[54] = 2 * x * y * z2;
    gr1[54] = x2 * z2;
    gr2[54] = 2 * x2 * y * z;
    h00[54] = 2 * y * z2;
    h01[54] = 2 * x * z2;
    h02[54] = 4 * x * y * z;
    h12[54] = 2 * x2 * z;
    h22[54] = 2 * x2 * y;
    XYZ[53] = x2 * y2 * z; // XXYYZ
    gr0[53] = 2 * x * y2 * z;
    gr1[53] = 2 * x2 * y * z;
    gr2[53] = x2 * y2;
    h00[53] = 2 * y2 * z;
    h01[53] = 4 * x * y * z;
    h02[53] = 2 * x * y2;
    h11[53] = 2 * x2 * z;
    h12[53] = 2 * x2 * y;
    XYZ[52] = x * y * z3; // ZZZXY
    gr0[52] = y * z3;
    gr1[52] = x * z3;
    gr2[52] = 3 * x * y * z2;
    h01[52] = z3;
    h02[52] = 3 * y * z2;
    h12[52] = 3 * x * z2;
    h22[52] = 6 * x * y * z;
    XYZ[51] = x * y3 * z; // YYYXZ
    gr0[51] = y3 * z;
    gr1[51] = 3 * x * y2 * z;
    gr2[51] = x * y3;
    h01[51] = 3 * y2 * z;
    h02[51] = y3;
    h11[51] = 6 * x * y * z;
    h12[51] = 3 * x * y2;
    XYZ[50] = x3 * y * z; // XXXYZ
    gr0[50] = 3 * x2 * y * z;
    gr1[50] = x3 * z;
    gr2[50] = x3 * y;
    h00[50] = 6 * x * y * z;
    h01[50] = 3 * x2 * z;
    h02[50] = 3 * x2 * y;
    h12[50] = x3;
    XYZ[49] = y2 * z3; // ZZZYY
    gr1[49] = 2 * y * z3;
    gr2[49] = 3 * y2 * z2;
    h11[49] = 2 * z3;
    h12[49] = 6 * y * z2;
    h22[49] = 6 * y2 * z;
    XYZ[48] = x2 * z3; // ZZZXX
    gr0[48] = 2 * x * z3;
    gr2[48] = 3 * x2 * z2;
    h00[48] = 2 * z3;
    h02[48] = 6 * x * z2;
    h22[48] = 6 * x2 * z;
    XYZ[47] = y3 * z2; // YYYZZ
    gr1[47] = 3 * y2 * z2;
    gr2[47] = 2 * y3 * z;
    h11[47] = 6 * y * z2;
    h12[47] = 6 * y2 * z;
    h22[47] = 2 * y3;
    XYZ[46] = x2 * y3; // YYYXX
    gr0[46] = 2 * x * y3;
    gr1[46] = 3 * x2 * y2;
    h00[46] = 2 * y3;
    h01[46] = 6 * x * y2;
    h11[46] = 6 * x2 * y;
    XYZ[45] = x3 * z2; // XXXZZ
    gr0[45] = 3 * x2 * z2;
    gr2[45] = 2 * x3 * z;
    h00[45] = 6 * x * z2;
    h02[45] = 6 * x2 * z;
    h22[45] = 2 * x3;
    XYZ[44] = x3 * y2; // XXXYY
    gr0[44] = 3 * x2 * y2;
    gr1[44] = 2 * x3 * y;
    h00[44] = 6 * x * y2;
    h01[44] = 6 * x2 * y;
    h11[44] = 2 * x3;
    XYZ[43] = y * z4; // ZZZZY
    gr1[43] = z4;
    gr2[43] = 4 * y * z3;
    h12[43] = 4 * z3;
    h22[43] = 12 * y * z2;
    XYZ[42] = x * z4; // ZZZZX
    gr0[42] = z4;
    gr2[42] = 4 * x * z3;
    h02[42] = 4 * z3;
    h22[42] = 12 * x * z2;
    XYZ[41] = y4 * z; // YYYYZ
    gr1[41] = 4 * y3 * z;
    gr2[41] = y4;
    h11[41] = 12 * y2 * z;
    h12[41] = 4 * y3;
    XYZ[40] = x * y4; // YYYYX
    gr0[40] = y4;
    gr1[40] = 4 * x * y3;
    h01[40] = 4 * y3;
    h11[40] = 12 * x * y2;
    XYZ[39] = x4 * z; // XXXXZ
    gr0[39] = 4 * x3 * z;
    gr2[39] = x4;
    h00[39] = 12 * x2 * z;
    h02[39] = 4 * x3;
    XYZ[38] = x4 * y; // XXXXY
    gr0[38] = 4 * x3 * y;
    gr1[38] = x4;
    h00[38] = 12 * x2 * y;
    h01[38] = 4 * x3;
    XYZ[37] = z * z4; // ZZZZZ
    gr2[37] = 5 * z4;
    h22[37] = 20 * z3;
    XYZ[36] = y * y4; // YYYYY
    gr1[36] = 5 * y4;
    h11[36] = 20 * y3;
    XYZ[35] = x * x4; // XXXXX
    gr0[35] = 5 * x4;
    h00[35] = 20 * x3;
  case 4:
    XYZ[34] = x * y * z2; // ZZXY
    gr0[34] = y * z2;
    gr1[34] = x * z2;
    gr2[34] = 2 * x * y * z;
    h01[34] = z2;
    h02[34] = 2 * y * z;
    h12[34] = 2 * x * z;
    h22[34] = 2 * x * y;
    XYZ[33] = x * y2 * z; // YYXZ
    gr0[33] = y2 * z;
    gr1[33] = 2 * x * y * z;
    gr2[33] = x * y2;
    h01[33] = 2 * y * z;
    h02[33] = y2;
    h11[33] = 2 * x * z;
    h12[33] = 2 * x * y;
    XYZ[32] = x2 * y * z; // XXYZ
    gr0[32] = 2 * x * y * z;
    gr1[32] = x2 * z;
    gr2[32] = x2 * y;
    h00[32] = 2 * y * z;
    h01[32] = 2 * x * z;
    h02[32] = 2 * x * y;
    h12[32] = x2;
    XYZ[31] = y2 * z2; // YYZZ
    gr1[31] = 2 * y * z2;
    gr2[31] = 2 * y2 * z;
    h11[31] = 2 * z2;
    h12[31] = 4 * y * z;
    h22[31] = 2 * y2;
    XYZ[30] = x2 * z2; // XXZZ
    gr0[30] = 2 * x * z2;
    gr2[30] = 2 * x2 * z;
    h00[30] = 2 * z2;
    h02[30] = 4 * x * z;
    h22[30] = 2 * x2;
    XYZ[29] = x2 * y2; // XXYY
    gr0[29] = 2 * x * y2;
    gr1[29] = 2 * x2 * y;
    h00[29] = 2 * y2;
    h01[29] = 4 * x * y;
    h11[29] = 2 * x2;
    XYZ[28] = y * z3; // ZZZY
    gr1[28] = z3;
    gr2[28] = 3 * y * z2;
    h12[28] = 3 * z2;
    h22[28] = 6 * y * z;
    XYZ[27] = x * z3; // ZZZX
    gr0[27] = z3;
    gr2[27] = 3 * x * z2;
    h02[27] = 3 * z2;
    h22[27] = 6 * x * z;
    XYZ[26] = y3 * z; // YYYZ
    gr1[26] = 3 * y2 * z;
    gr2[26] = y3;
    h11[26] = 6 * y * z;
    h12[26] = 3 * y2;
    XYZ[25] = x * y3; // YYYX
    gr0[25] = y3;
    gr1[25] = 3 * x * y2;
    h01[25] = 3 * y2;
    h11[25] = 6 * x * y;
    XYZ[24] = x3 * z; // XXXZ
    gr0[24] = 3 * x2 * z;
    gr2[24] = x3;
    h00[24] = 6 * x * z;
    h02[24] = 3 * x2;
    XYZ[23] = x3 * y; // XXXY
    gr0[23] = 3 * x2 * y;
    gr1[23] = x3;
    h00[23] = 6 * x * y;
    h01[23] = 3 * x2;
    XYZ[22] = z4; // ZZZZ
    gr2[22] = 4 * z3;
    h22[22] = 12 * z2;
    XYZ[21] = y4; // YYYY
    gr1[21] = 4 * y3;
    h11[21] = 12 * y2;
    XYZ[20] = x4; // XXXX
    gr0[20] = 4 * x3;
    h00[20] = 12 * x2;
  case 3:
    XYZ[19] = x * y * z; // XYZ
    gr0[19] = y * z;
    gr1[19] = x * z;
    gr2[19] = x * y;
    h01[19] = z;
    h02[19] = y;
    h12[19] = x;
    XYZ[18] = y * z2; // ZZY
    gr1[18] = z2;
    gr2[18] = 2 * y * z;
    h12[18] = 2 * z;
    h22[18] = 2 * y;
    XYZ[17] = x * z2; // ZZX
    gr0[17] = z2;
    gr2[17] = 2 * x * z;
    h02[17] = 2 * z;
    h22[17] = 2 * x;
    XYZ[16] = y2 * z; // YYZ
    gr1[16] = 2 * y * z;
    gr2[16] = y2;
    h11[16] = 2 * z;
    h12[16] = 2 * y;
    XYZ[15] = x * y2; // YYX
    gr0[15] = y2;
    gr1[15] = 2 * x * y;
    h01[15] = 2 * y;
    h11[15] = 2 * x;
    XYZ[14] = x2 * z; // XXZ
    gr0[14] = 2 * x * z;
    gr2[14] = x2;
    h00[14] = 2 * z;
    h02[14] = 2 * x;
    XYZ[13] = x2 * y; // XXY
    gr0[13] = 2 * x * y;
    gr1[13] = x2;
    h00[13] = 2 * y;
    h01[13] = 2 * x;
    XYZ[12] = z3; // ZZZ
    gr2[12] = 3 * z2;
    h22[12] = 6 * z;
    XYZ[11] = y3; // YYY
    gr1[11] = 3 * y2;
    h11[11] = 6 * y;
    XYZ[10] = x3; // XXX
    gr0[10] = 3 * x2;
    h00[10] = 6 * x;
  case 2:
    XYZ[9] = y * z; // YZ
    gr1[9] = z;
    gr2[9] = y;
    h12[9] = 1;
    XYZ[8] = x * z; // XZ
    gr0[8] = z;
    gr2[8] = x;
    h02[8] = 1;
    XYZ[7] = x * y; // XY
    gr0[7] = y;
    gr1[7] = x;
    h01[7] = 1;
    XYZ[6] = z2; // ZZ
    gr2[6] = 2 * z;
    h22[6] = 2;
    XYZ[5] = y2; // YY
    gr1[5] = 2 * y;
    h11[5] = 2;
    XYZ[4] = x2; // XX
    gr0[4] = 2 * x;
    h00[4] = 2;
  case 1:
    XYZ[3] = z; // Z
    gr2[3] = 1;
    XYZ[2] = y; // Y
    gr1[2] = 1;
    XYZ[1] = x; // X
    gr0[1] = 1;
  case 0:
    XYZ[0] = 1; // S
  }

  const size_t ntot = cXYZ.size();
  for (size_t i = 0; i < ntot; ++i)
  {
    XYZ[i] *= NormFactor[i];
    gr0[i] *= NormFactor[i];
    gr1[i] *= NormFactor[i];
    gr2[i] *= NormFactor[i];
    h00[i] *= NormFactor[i];
    h01[i] *= NormFactor[i];
    h02[i] *= NormFactor[i];
    h11[i] *= NormFactor[i];
    h12[i] *= NormFactor[i];
    h22[i] *= NormFactor[i];
  }
}


template<class T>
void SoaCartesianTensor<T>::evaluateVGHGH(T x, T y, T z)
{
  constexpr T czero(0);
  cXYZ = czero;

  const T x2 = x * x, y2 = y * y, z2 = z * z;
  const T x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
  const T x4 = x3 * x, y4 = y3 * y, z4 = z3 * z;
  const T x5 = x4 * x, y5 = y4 * y, z5 = z4 * z;

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

  switch (Lmax)
  {
  case 6:
    XYZ[83]   = x2 * y2 * z2; // X2Y2Z2
    gr0[83]   = 2 * x * y2 * z2;
    gr1[83]   = 2 * x2 * y * z2;
    gr2[83]   = 2 * x2 * y2 * z;
    h00[83]   = 2 * y2 * z2;
    h01[83]   = 4 * x * y * z2;
    h02[83]   = 4 * x * y2 * z;
    h11[83]   = 2 * x2 * z2;
    h12[83]   = 4 * x2 * y * z;
    h22[83]   = 2 * x2 * y2;
    gh001[83] = 4 * y * z2;
    gh002[83] = 4 * y2 * z;
    gh011[83] = 4 * x * z2;
    gh012[83] = 8 * x * y * z;
    gh022[83] = 4 * x * y2;
    gh112[83] = 4 * x2 * z;
    gh122[83] = 4 * x2 * y;
    XYZ[82]   = x * y2 * z3; // Z3Y2X
    gr0[82]   = y2 * z3;
    gr1[82]   = 2 * x * y * z3;
    gr2[82]   = 3 * x * y2 * z2;
    h01[82]   = 2 * y * z3;
    h02[82]   = 3 * y2 * z2;
    h11[82]   = 2 * x * z3;
    h12[82]   = 6 * x * y * z2;
    h22[82]   = 6 * x * y2 * z;
    gh011[82] = 2 * z3;
    gh012[82] = 6 * y * z2;
    gh022[82] = 6 * y2 * z;
    gh112[82] = 6 * x * z2;
    gh122[82] = 12 * x * y * z;
    gh222[82] = 6 * x * y2;
    XYZ[81]   = x2 * y * z3; // Z3X2Y
    gr0[81]   = 2 * x * y * z3;
    gr1[81]   = x2 * z3;
    gr2[81]   = 3 * x2 * y * z2;
    h00[81]   = 2 * y * z3;
    h01[81]   = 2 * x * z3;
    h02[81]   = 6 * x * y * z2;
    h12[81]   = 3 * x2 * z2;
    h22[81]   = 6 * x2 * y * z;
    gh001[81] = 2 * z3;
    gh002[81] = 6 * y * z2;
    gh012[81] = 6 * x * z2;
    gh022[81] = 12 * x * y * z;
    gh122[81] = 6 * x2 * z;
    gh222[81] = 6 * x2 * y;
    XYZ[80]   = x * y3 * z2; // Y3Z2X
    gr0[80]   = y3 * z2;
    gr1[80]   = 3 * x * y2 * z2;
    gr2[80]   = 2 * x * y3 * z;
    h01[80]   = 3 * y2 * z2;
    h02[80]   = 2 * y3 * z;
    h11[80]   = 6 * x * y * z2;
    h12[80]   = 6 * x * y2 * z;
    h22[80]   = 2 * x * y3;
    gh011[80] = 6 * y * z2;
    gh012[80] = 6 * y2 * z;
    gh022[80] = 2 * y3;
    gh111[80] = 6 * x * z2;
    gh112[80] = 12 * x * y * z;
    gh122[80] = 6 * x * y2;
    XYZ[79]   = x2 * y3 * z; // Y3X2Z
    gr0[79]   = 2 * x * y3 * z;
    gr1[79]   = 3 * x2 * y2 * z;
    gr2[79]   = x2 * y3;
    h00[79]   = 2 * y3 * z;
    h01[79]   = 6 * x * y2 * z;
    h02[79]   = 2 * x * y3;
    h11[79]   = 6 * x2 * y * z;
    h12[79]   = 3 * x2 * y2;
    gh001[79] = 6 * y2 * z;
    gh002[79] = 2 * y3;
    gh011[79] = 12 * x * y * z;
    gh012[79] = 6 * x * y2;
    gh111[79] = 6 * x2 * z;
    gh112[79] = 6 * x2 * y;
    XYZ[78]   = x3 * y * z2; // X3Z2Y
    gr0[78]   = 3 * x2 * y * z2;
    gr1[78]   = x3 * z2;
    gr2[78]   = 2 * x3 * y * z;
    h00[78]   = 6 * x * y * z2;
    h01[78]   = 3 * x2 * z2;
    h02[78]   = 6 * x2 * y * z;
    h12[78]   = 2 * x3 * z;
    h22[78]   = 2 * x3 * y;
    gh000[78] = 6 * y * z2;
    gh001[78] = 6 * x * z2;
    gh002[78] = 12 * x * y * z;
    gh012[78] = 6 * x2 * z;
    gh022[78] = 6 * x2 * y;
    gh122[78] = 2 * x3;
    XYZ[77]   = x3 * y2 * z; // X3Y2Z
    gr0[77]   = 3 * x2 * y2 * z;
    gr1[77]   = 2 * x3 * y * z;
    gr2[77]   = x3 * y2;
    h00[77]   = 6 * x * y2 * z;
    h01[77]   = 6 * x2 * y * z;
    h02[77]   = 3 * x2 * y2;
    h11[77]   = 2 * x3 * z;
    h12[77]   = 2 * x3 * y;
    gh000[77] = 6 * y2 * z;
    gh001[77] = 12 * x * y * z;
    gh002[77] = 6 * x * y2;
    gh011[77] = 6 * x2 * z;
    gh012[77] = 6 * x2 * y;
    gh112[77] = 2 * x3;
    XYZ[76]   = y3 * z3; // Y3Z3
    gr1[76]   = 3 * y2 * z3;
    gr2[76]   = 3 * y3 * z2;
    h11[76]   = 6 * y * z3;
    h12[76]   = 9 * y2 * z2;
    h22[76]   = 6 * y3 * z;
    gh111[76] = 6 * z3;
    gh112[76] = 18 * y * z2;
    gh122[76] = 18 * y2 * z;
    gh222[76] = 6 * y3;
    XYZ[75]   = x3 * z3; // X3Z3
    gr0[75]   = 3 * x2 * z3;
    gr2[75]   = 3 * x3 * z2;
    h00[75]   = 6 * x * z3;
    h02[75]   = 9 * x2 * z2;
    h22[75]   = 6 * x3 * z;
    gh000[75] = 6 * z3;
    gh002[75] = 18 * x * z2;
    gh022[75] = 18 * x2 * z;
    gh222[75] = 6 * x3;
    XYZ[74]   = x3 * y3; // X3Y3
    gr0[74]   = 3 * x2 * y3;
    gr1[74]   = 3 * x3 * y2;
    h00[74]   = 6 * x * y3;
    h01[74]   = 9 * x2 * y2;
    h11[74]   = 6 * x3 * y;
    gh000[74] = 6 * y3;
    gh001[74] = 18 * x * y2;
    gh011[74] = 18 * x2 * y;
    gh111[74] = 6 * x3;
    XYZ[73]   = x * y * z4; // Z4XY
    gr0[73]   = y * z4;
    gr1[73]   = x * z4;
    gr2[73]   = 4 * x * y * z3;
    h01[73]   = z4;
    h02[73]   = 4 * y * z3;
    h12[73]   = 4 * x * z3;
    h22[73]   = 12 * x * y * z2;
    gh012[73] = 4 * z3;
    gh022[73] = 12 * y * z2;
    gh122[73] = 12 * x * z2;
    gh222[73] = 24 * x * y * z;
    XYZ[72]   = x * y4 * z; // Y4XZ
    gr0[72]   = y4 * z;
    gr1[72]   = 4 * x * y3 * z;
    gr2[72]   = x * y4;
    h01[72]   = 4 * y3 * z;
    h02[72]   = y4;
    h11[72]   = 12 * x * y2 * z;
    h12[72]   = 4 * x * y3;
    gh011[72] = 12 * y2 * z;
    gh012[72] = 4 * y3;
    gh111[72] = 24 * x * y * z;
    gh112[72] = 12 * x * y2;
    XYZ[71]   = x4 * y * z; // X4YZ
    gr0[71]   = 4 * x3 * y * z;
    gr1[71]   = x4 * z;
    gr2[71]   = x4 * y;
    h00[71]   = 12 * x2 * y * z;
    h01[71]   = 4 * x3 * z;
    h02[71]   = 4 * x3 * y;
    h12[71]   = x4;
    gh000[71] = 24 * x * y * z;
    gh001[71] = 12 * x2 * z;
    gh002[71] = 12 * x2 * y;
    gh012[71] = 4 * x3;
    XYZ[70]   = y2 * z4; // Z4Y2
    gr1[70]   = 2 * y * z4;
    gr2[70]   = 4 * y2 * z3;
    h11[70]   = 2 * z4;
    h12[70]   = 8 * y * z3;
    h22[70]   = 12 * y2 * z2;
    gh112[70] = 8 * z3;
    gh122[70] = 24 * y * z2;
    gh222[70] = 24 * y2 * z;
    XYZ[69]   = x2 * z4; // Z4X2
    gr0[69]   = 2 * x * z4;
    gr2[69]   = 4 * x2 * z3;
    h00[69]   = 2 * z4;
    h02[69]   = 8 * x * z3;
    h22[69]   = 12 * x2 * z2;
    gh002[69] = 8 * z3;
    gh022[69] = 24 * x * z2;
    gh222[69] = 24 * x2 * z;
    XYZ[68]   = y4 * z2; // Y4Z2
    gr1[68]   = 4 * y3 * z2;
    gr2[68]   = 2 * y4 * z;
    h11[68]   = 12 * y2 * z2;
    h12[68]   = 8 * y3 * z;
    h22[68]   = 2 * y4;
    gh111[68] = 24 * y * z2;
    gh112[68] = 24 * y2 * z;
    gh122[68] = 8 * y3;
    XYZ[67]   = x2 * y4; // Y4X2
    gr0[67]   = 2 * x * y4;
    gr1[67]   = 4 * x2 * y3;
    h00[67]   = 2 * y4;
    h01[67]   = 8 * x * y3;
    h11[67]   = 12 * x2 * y2;
    gh001[67] = 8 * y3;
    gh011[67] = 24 * x * y2;
    gh111[67] = 24 * x2 * y;
    XYZ[66]   = x4 * z2; // X4Z2
    gr0[66]   = 4 * x3 * z2;
    gr2[66]   = 2 * x4 * z;
    h00[66]   = 12 * x2 * z2;
    h02[66]   = 8 * x3 * z;
    h22[66]   = 2 * x4;
    gh000[66] = 24 * x * z2;
    gh002[66] = 24 * x2 * z;
    gh022[66] = 8 * x3;
    XYZ[65]   = x4 * y2; // X4Y2
    gr0[65]   = 4 * x3 * y2;
    gr1[65]   = 2 * x4 * y;
    h00[65]   = 12 * x2 * y2;
    h01[65]   = 8 * x3 * y;
    h11[65]   = 2 * x4;
    gh000[65] = 24 * x * y2;
    gh001[65] = 24 * x2 * y;
    gh011[65] = 8 * x3;
    XYZ[64]   = y * z * z4; // Z5Y
    gr1[64]   = z * z4;
    gr2[64]   = 5 * y * z4;
    h12[64]   = 5 * z4;
    h22[64]   = 20 * y * z3;
    gh122[64] = 20 * z3;
    gh222[64] = 60 * y * z2;
    XYZ[63]   = x * z * z4; // Z5X
    gr0[63]   = z * z4;
    gr2[63]   = 5 * x * z4;
    h02[63]   = 5 * z4;
    h22[63]   = 20 * x * z3;
    gh022[63] = 20 * z3;
    gh222[63] = 60 * x * z2;
    XYZ[62]   = y * y4 * z; // Y5Z
    gr1[62]   = 5 * y4 * z;
    gr2[62]   = y * y4;
    h11[62]   = 20 * y3 * z;
    h12[62]   = 5 * y4;
    gh111[62] = 60 * y2 * z;
    gh112[62] = 20 * y3;
    XYZ[61]   = x * y * y4; // Y5X
    gr0[61]   = y * y4;
    gr1[61]   = 5 * x * y4;
    h01[61]   = 5 * y4;
    h11[61]   = 20 * x * y3;
    gh011[61] = 20 * y3;
    gh111[61] = 60 * x * y2;
    XYZ[60]   = x * x4 * z; // X5Z
    gr0[60]   = 5 * x4 * z;
    gr2[60]   = x * x4;
    h00[60]   = 20 * x3 * z;
    h02[60]   = 5 * x4;
    gh000[60] = 60 * x2 * z;
    gh002[60] = 20 * x3;
    XYZ[59]   = x * x4 * y; // X5Y
    gr0[59]   = 5 * x4 * y;
    gr1[59]   = x * x4;
    h00[59]   = 20 * x3 * y;
    h01[59]   = 5 * x4;
    gh000[59] = 60 * x2 * y;
    gh001[59] = 20 * x3;
    XYZ[58]   = z * z5; // Z6
    gr2[58]   = 6 * z * z4;
    h22[58]   = 30 * z4;
    gh222[58] = 120 * z3;
    XYZ[57]   = y * y5; // Y6
    gr1[57]   = 6 * y * y4;
    h11[57]   = 30 * y4;
    gh111[57] = 120 * y3;
    XYZ[56]   = x * x5; // X6
    gr0[56]   = 6 * x * x4;
    h00[56]   = 30 * x4;
    gh000[56] = 120 * x3;
  case 5:
    XYZ[55]   = x * y2 * z2; // YYZZX
    gr0[55]   = y2 * z2;
    gr1[55]   = 2 * x * y * z2;
    gr2[55]   = 2 * x * y2 * z;
    h01[55]   = 2 * y * z2;
    h02[55]   = 2 * y2 * z;
    h11[55]   = 2 * x * z2;
    h12[55]   = 4 * x * y * z;
    h22[55]   = 2 * x * y2;
    gh011[55] = 2 * z2;
    gh012[55] = 4 * y * z;
    gh022[55] = 2 * y2;
    gh112[55] = 4 * x * z;
    gh122[55] = 4 * x * y;
    XYZ[54]   = x2 * y * z2; // XXZZY
    gr0[54]   = 2 * x * y * z2;
    gr1[54]   = x2 * z2;
    gr2[54]   = 2 * x2 * y * z;
    h00[54]   = 2 * y * z2;
    h01[54]   = 2 * x * z2;
    h02[54]   = 4 * x * y * z;
    h12[54]   = 2 * x2 * z;
    h22[54]   = 2 * x2 * y;
    gh001[54] = 2 * z2;
    gh002[54] = 4 * y * z;
    gh012[54] = 4 * x * z;
    gh022[54] = 4 * x * y;
    gh122[54] = 2 * x2;
    XYZ[53]   = x2 * y2 * z; // XXYYZ
    gr0[53]   = 2 * x * y2 * z;
    gr1[53]   = 2 * x2 * y * z;
    gr2[53]   = x2 * y2;
    h00[53]   = 2 * y2 * z;
    h01[53]   = 4 * x * y * z;
    h02[53]   = 2 * x * y2;
    h11[53]   = 2 * x2 * z;
    h12[53]   = 2 * x2 * y;
    gh001[53] = 4 * y * z;
    gh002[53] = 2 * y2;
    gh011[53] = 4 * x * z;
    gh012[53] = 4 * x * y;
    gh112[53] = 2 * x2;
    XYZ[52]   = x * y * z3; // ZZZXY
    gr0[52]   = y * z3;
    gr1[52]   = x * z3;
    gr2[52]   = 3 * x * y * z2;
    h01[52]   = z3;
    h02[52]   = 3 * y * z2;
    h12[52]   = 3 * x * z2;
    h22[52]   = 6 * x * y * z;
    gh012[52] = 3 * z2;
    gh022[52] = 6 * y * z;
    gh122[52] = 6 * x * z;
    gh222[52] = 6 * x * y;
    XYZ[51]   = x * y3 * z; // YYYXZ
    gr0[51]   = y3 * z;
    gr1[51]   = 3 * x * y2 * z;
    gr2[51]   = x * y3;
    h01[51]   = 3 * y2 * z;
    h02[51]   = y3;
    h11[51]   = 6 * x * y * z;
    h12[51]   = 3 * x * y2;
    gh011[51] = 6 * y * z;
    gh012[51] = 3 * y2;
    gh111[51] = 6 * x * z;
    gh112[51] = 6 * x * y;
    XYZ[50]   = x3 * y * z; // XXXYZ
    gr0[50]   = 3 * x2 * y * z;
    gr1[50]   = x3 * z;
    gr2[50]   = x3 * y;
    h00[50]   = 6 * x * y * z;
    h01[50]   = 3 * x2 * z;
    h02[50]   = 3 * x2 * y;
    h12[50]   = x3;
    gh000[50] = 6 * y * z;
    gh001[50] = 6 * x * z;
    gh002[50] = 6 * x * y;
    gh012[50] = 3 * x2;
    XYZ[49]   = y2 * z3; // ZZZYY
    gr1[49]   = 2 * y * z3;
    gr2[49]   = 3 * y2 * z2;
    h11[49]   = 2 * z3;
    h12[49]   = 6 * y * z2;
    h22[49]   = 6 * y2 * z;
    gh112[49] = 6 * z2;
    gh122[49] = 12 * y * z;
    gh222[49] = 6 * y2;
    XYZ[48]   = x2 * z3; // ZZZXX
    gr0[48]   = 2 * x * z3;
    gr2[48]   = 3 * x2 * z2;
    h00[48]   = 2 * z3;
    h02[48]   = 6 * x * z2;
    h22[48]   = 6 * x2 * z;
    gh002[48] = 6 * z2;
    gh022[48] = 12 * x * z;
    gh222[48] = 6 * x2;
    XYZ[47]   = y3 * z2; // YYYZZ
    gr1[47]   = 3 * y2 * z2;
    gr2[47]   = 2 * y3 * z;
    h11[47]   = 6 * y * z2;
    h12[47]   = 6 * y2 * z;
    h22[47]   = 2 * y3;
    gh111[47] = 6 * z2;
    gh112[47] = 12 * y * z;
    gh122[47] = 6 * y2;
    XYZ[46]   = x2 * y3; // YYYXX
    gr0[46]   = 2 * x * y3;
    gr1[46]   = 3 * x2 * y2;
    h00[46]   = 2 * y3;
    h01[46]   = 6 * x * y2;
    h11[46]   = 6 * x2 * y;
    gh001[46] = 6 * y2;
    gh011[46] = 12 * x * y;
    gh111[46] = 6 * x2;
    XYZ[45]   = x3 * z2; // XXXZZ
    gr0[45]   = 3 * x2 * z2;
    gr2[45]   = 2 * x3 * z;
    h00[45]   = 6 * x * z2;
    h02[45]   = 6 * x2 * z;
    h22[45]   = 2 * x3;
    gh000[45] = 6 * z2;
    gh002[45] = 12 * x * z;
    gh022[45] = 6 * x2;
    XYZ[44]   = x3 * y2; // XXXYY
    gr0[44]   = 3 * x2 * y2;
    gr1[44]   = 2 * x3 * y;
    h00[44]   = 6 * x * y2;
    h01[44]   = 6 * x2 * y;
    h11[44]   = 2 * x3;
    gh000[44] = 6 * y2;
    gh001[44] = 12 * x * y;
    gh011[44] = 6 * x2;
    XYZ[43]   = y * z4; // ZZZZY
    gr1[43]   = z4;
    gr2[43]   = 4 * y * z3;
    h12[43]   = 4 * z3;
    h22[43]   = 12 * y * z2;
    gh122[43] = 12 * z2;
    gh222[43] = 24 * y * z;
    XYZ[42]   = x * z4; // ZZZZX
    gr0[42]   = z4;
    gr2[42]   = 4 * x * z3;
    h02[42]   = 4 * z3;
    h22[42]   = 12 * x * z2;
    gh022[42] = 12 * z2;
    gh222[42] = 24 * x * z;
    XYZ[41]   = y4 * z; // YYYYZ
    gr1[41]   = 4 * y3 * z;
    gr2[41]   = y4;
    h11[41]   = 12 * y2 * z;
    h12[41]   = 4 * y3;
    gh111[41] = 24 * y * z;
    gh112[41] = 12 * y2;
    XYZ[40]   = x * y4; // YYYYX
    gr0[40]   = y4;
    gr1[40]   = 4 * x * y3;
    h01[40]   = 4 * y3;
    h11[40]   = 12 * x * y2;
    gh011[40] = 12 * y2;
    gh111[40] = 24 * x * y;
    XYZ[39]   = x4 * z; // XXXXZ
    gr0[39]   = 4 * x3 * z;
    gr2[39]   = x4;
    h00[39]   = 12 * x2 * z;
    h02[39]   = 4 * x3;
    gh000[39] = 24 * x * z;
    gh002[39] = 12 * x2;
    XYZ[38]   = x4 * y; // XXXXY
    gr0[38]   = 4 * x3 * y;
    gr1[38]   = x4;
    h00[38]   = 12 * x2 * y;
    h01[38]   = 4 * x3;
    gh000[38] = 24 * x * y;
    gh001[38] = 12 * x2;
    XYZ[37]   = z * z4; // ZZZZZ
    gr2[37]   = 5 * z4;
    h22[37]   = 20 * z3;
    gh222[37] = 60 * z2;
    XYZ[36]   = y * y4; // YYYYY
    gr1[36]   = 5 * y4;
    h11[36]   = 20 * y3;
    gh111[36] = 60 * y2;
    XYZ[35]   = x * x4; // XXXXX
    gr0[35]   = 5 * x4;
    h00[35]   = 20 * x3;
    gh000[35] = 60 * x2;
  case 4:
    XYZ[34]   = x * y * z2; // ZZXY
    gr0[34]   = y * z2;
    gr1[34]   = x * z2;
    gr2[34]   = 2 * x * y * z;
    h01[34]   = z2;
    h02[34]   = 2 * y * z;
    h12[34]   = 2 * x * z;
    h22[34]   = 2 * x * y;
    gh012[34] = 2 * z;
    gh022[34] = 2 * y;
    gh122[34] = 2 * x;
    XYZ[33]   = x * y2 * z; // YYXZ
    gr0[33]   = y2 * z;
    gr1[33]   = 2 * x * y * z;
    gr2[33]   = x * y2;
    h01[33]   = 2 * y * z;
    h02[33]   = y2;
    h11[33]   = 2 * x * z;
    h12[33]   = 2 * x * y;
    gh011[33] = 2 * z;
    gh012[33] = 2 * y;
    gh112[33] = 2 * x;
    XYZ[32]   = x2 * y * z; // XXYZ
    gr0[32]   = 2 * x * y * z;
    gr1[32]   = x2 * z;
    gr2[32]   = x2 * y;
    h00[32]   = 2 * y * z;
    h01[32]   = 2 * x * z;
    h02[32]   = 2 * x * y;
    h12[32]   = x2;
    gh001[32] = 2 * z;
    gh002[32] = 2 * y;
    gh012[32] = 2 * x;
    XYZ[31]   = y2 * z2; // YYZZ
    gr1[31]   = 2 * y * z2;
    gr2[31]   = 2 * y2 * z;
    h11[31]   = 2 * z2;
    h12[31]   = 4 * y * z;
    h22[31]   = 2 * y2;
    gh112[31] = 4 * z;
    gh122[31] = 4 * y;
    XYZ[30]   = x2 * z2; // XXZZ
    gr0[30]   = 2 * x * z2;
    gr2[30]   = 2 * x2 * z;
    h00[30]   = 2 * z2;
    h02[30]   = 4 * x * z;
    h22[30]   = 2 * x2;
    gh002[30] = 4 * z;
    gh022[30] = 4 * x;
    XYZ[29]   = x2 * y2; // XXYY
    gr0[29]   = 2 * x * y2;
    gr1[29]   = 2 * x2 * y;
    h00[29]   = 2 * y2;
    h01[29]   = 4 * x * y;
    h11[29]   = 2 * x2;
    gh001[29] = 4 * y;
    gh011[29] = 4 * x;
    XYZ[28]   = y * z3; // ZZZY
    gr1[28]   = z3;
    gr2[28]   = 3 * y * z2;
    h12[28]   = 3 * z2;
    h22[28]   = 6 * y * z;
    gh122[28] = 6 * z;
    gh222[28] = 6 * y;
    XYZ[27]   = x * z3; // ZZZX
    gr0[27]   = z3;
    gr2[27]   = 3 * x * z2;
    h02[27]   = 3 * z2;
    h22[27]   = 6 * x * z;
    gh022[27] = 6 * z;
    gh222[27] = 6 * x;
    XYZ[26]   = y3 * z; // YYYZ
    gr1[26]   = 3 * y2 * z;
    gr2[26]   = y3;
    h11[26]   = 6 * y * z;
    h12[26]   = 3 * y2;
    gh111[26] = 6 * z;
    gh112[26] = 6 * y;
    XYZ[25]   = x * y3; // YYYX
    gr0[25]   = y3;
    gr1[25]   = 3 * x * y2;
    h01[25]   = 3 * y2;
    h11[25]   = 6 * x * y;
    gh011[25] = 6 * y;
    gh111[25] = 6 * x;
    XYZ[24]   = x3 * z; // XXXZ
    gr0[24]   = 3 * x2 * z;
    gr2[24]   = x3;
    h00[24]   = 6 * x * z;
    h02[24]   = 3 * x2;
    gh000[24] = 6 * z;
    gh002[24] = 6 * x;
    XYZ[23]   = x3 * y; // XXXY
    gr0[23]   = 3 * x2 * y;
    gr1[23]   = x3;
    h00[23]   = 6 * x * y;
    h01[23]   = 3 * x2;
    gh000[23] = 6 * y;
    gh001[23] = 6 * x;
    XYZ[22]   = z4; // ZZZZ
    gr2[22]   = 4 * z3;
    h22[22]   = 12 * z2;
    gh222[22] = 24 * z;
    XYZ[21]   = y4; // YYYY
    gr1[21]   = 4 * y3;
    h11[21]   = 12 * y2;
    gh111[21] = 24 * y;
    XYZ[20]   = x4; // XXXX
    gr0[20]   = 4 * x3;
    h00[20]   = 12 * x2;
    gh000[20] = 24 * x;
  case 3:
    XYZ[19]   = x * y * z; // XYZ
    gr0[19]   = y * z;
    gr1[19]   = x * z;
    gr2[19]   = x * y;
    h01[19]   = z;
    h02[19]   = y;
    h12[19]   = x;
    gh012[19] = 1;
    XYZ[18]   = y * z2; // ZZY
    gr1[18]   = z2;
    gr2[18]   = 2 * y * z;
    h12[18]   = 2 * z;
    h22[18]   = 2 * y;
    gh122[18] = 2;
    XYZ[17]   = x * z2; // ZZX
    gr0[17]   = z2;
    gr2[17]   = 2 * x * z;
    h02[17]   = 2 * z;
    h22[17]   = 2 * x;
    gh022[17] = 2;
    XYZ[16]   = y2 * z; // YYZ
    gr1[16]   = 2 * y * z;
    gr2[16]   = y2;
    h11[16]   = 2 * z;
    h12[16]   = 2 * y;
    gh112[16] = 2;
    XYZ[15]   = x * y2; // YYX
    gr0[15]   = y2;
    gr1[15]   = 2 * x * y;
    h01[15]   = 2 * y;
    h11[15]   = 2 * x;
    gh011[15] = 2;
    XYZ[14]   = x2 * z; // XXZ
    gr0[14]   = 2 * x * z;
    gr2[14]   = x2;
    h00[14]   = 2 * z;
    h02[14]   = 2 * x;
    gh002[14] = 2;
    XYZ[13]   = x2 * y; // XXY
    gr0[13]   = 2 * x * y;
    gr1[13]   = x2;
    h00[13]   = 2 * y;
    h01[13]   = 2 * x;
    gh001[13] = 2;
    XYZ[12]   = z3; // ZZZ
    gr2[12]   = 3 * z2;
    h22[12]   = 6 * z;
    gh222[12] = 6;
    XYZ[11]   = y3; // YYY
    gr1[11]   = 3 * y2;
    h11[11]   = 6 * y;
    gh111[11] = 6;
    XYZ[10]   = x3; // XXX
    gr0[10]   = 3 * x2;
    h00[10]   = 6 * x;
    gh000[10] = 6;
  case 2:
    XYZ[9] = y * z; // YZ
    gr1[9] = z;
    gr2[9] = y;
    h12[9] = 1;
    XYZ[8] = x * z; // XZ
    gr0[8] = z;
    gr2[8] = x;
    h02[8] = 1;
    XYZ[7] = x * y; // XY
    gr0[7] = y;
    gr1[7] = x;
    h01[7] = 1;
    XYZ[6] = z2; // ZZ
    gr2[6] = 2 * z;
    h22[6] = 2;
    XYZ[5] = y2; // YY
    gr1[5] = 2 * y;
    h11[5] = 2;
    XYZ[4] = x2; // XX
    gr0[4] = 2 * x;
    h00[4] = 2;
  case 1:
    XYZ[3] = z; // Z
    gr2[3] = 1;
    XYZ[2] = y; // Y
    gr1[2] = 1;
    XYZ[1] = x; // X
    gr0[1] = 1;
  case 0:
    XYZ[0] = 1; // S
  }

  const size_t ntot = cXYZ.size();
  for (size_t i = 0; i < ntot; ++i)
  {
    XYZ[i] *= NormFactor[i];
    gr0[i] *= NormFactor[i];
    gr1[i] *= NormFactor[i];
    gr2[i] *= NormFactor[i];
    h00[i] *= NormFactor[i];
    h01[i] *= NormFactor[i];
    h02[i] *= NormFactor[i];
    h11[i] *= NormFactor[i];
    h12[i] *= NormFactor[i];
    h22[i] *= NormFactor[i];
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

// generated from read_order.py
template<class T>
void SoaCartesianTensor<T>::getABC(int n, int& a, int& b, int& c)
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
    throw std::runtime_error("CartesianTensor::getABC() - Incorrect index.\n");
    break;
  }
}

} //namespace qmcplusplus
#endif
