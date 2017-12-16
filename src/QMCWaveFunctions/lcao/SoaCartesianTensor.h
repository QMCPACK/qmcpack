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
    
    



#ifndef QMCPLUSPLUS_SOA_CARTESIAN_TENSOR_H
#define QMCPLUSPLUS_SOA_CARTESIAN_TENSOR_H

#include "OhmmsSoA/Container.h"

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
  typedef T value_type;
  typedef TinyVector<Tensor<T,3>,3> ggg_type;

  ///maximum angular momentum
  size_t Lmax;
  ///normalization factor
  aligned_vector<T> NormFactor;
  ///composite V,Gx,Gy,Gz,[L | H00, H01, H02, H11, H12, H12]
  VectorSoaContainer<T,10> cXYZ;
  ///third derivative: keep the TinyyVector<Tensor<T,3>,3> and AoS for now
  std::vector<ggg_type> gggXYZ;
  //std::vector<ggg_type> gggXYZ;
  //VectorSoAContainer<T,18> gggXYZ;

  /** constructor
   * @param l_max maximum angular momentum
   *
   * Evaluate all the constants and prefactors.
  */
  explicit SoaCartesianTensor(const int l_max, bool addsign=false);

  ///compute Ylm
  void evaluate_bare(T x, T y, T z, T* XYZ) const;

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* XYZ) const
  {
    evaluate_bare(x,y,z,XYZ);
    for (size_t i=0, nl=cXYZ.size(); i<nl; i++)
      XYZ[i]*= NormFactor[i];
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* XYZ) 
  {
    evaluate_bare(x,y,z,XYZ);
    for (size_t i=0, nl=cXYZ.size(); i<nl; i++)
      XYZ[i]*= NormFactor[i];
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z)
  {
    evaluateV(x,y,z,cXYZ.data(0));
  }

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGL(T x, T y, T z);

  void evaluateVGH(T x, T y, T z);

  ///returns dummy: this is not used
  inline int index(int l, int m) const
  {
    return (l*(l+1))+m;
  }

  /** return the starting address of the component
   *
   * component=0(V), 1(dx), 2(dy), 3(dz), 4(Lap)
   */
  inline const T* operator[](size_t component) const
  {
    return cXYZ.data(component);
  }

  inline size_t size() const
  {
    return cXYZ.size();
  }

  inline int lmax() const
  {
    return Lmax;
  }

  inline void getABC(int n, int& a, int& b, int& c);

  int DFactorial(int num)
  {
    return (num<2)? 1: num*DFactorial(num-2);
  }

};

template<class T>
SoaCartesianTensor<T>::SoaCartesianTensor(const int l_max, bool addsign) : Lmax(l_max)
{
  if(Lmax < 0 || Lmax > 4)
  {
    std::cerr <<"CartesianTensor can't handle Lmax > 4 or Lmax < 0.\n";
    APP_ABORT("");
  }
  int ntot = 0;
  for(int i=0; i<=Lmax; i++)
    ntot+=(i+1)*(i+2)/2;
  cXYZ.resize(ntot);
  NormFactor.resize(ntot,1);
  int p=0;
  int a,b,c;
  const double pi = 4.0*std::atan(1.0);
  for(int l=0; l<=Lmax; l++)
  {
    int n = (l+1)*(l+2)/2;
    for(int k=0; k<n; k++)
    {
      getABC(p,a,b,c);
// factor of (alpha^(l+3/2))^(1/2) goes into the radial function
// mmorales: HACK HACK HACK, to avoid modifyng the radial functions,
//           I add a term to the normalization to cancel the term
//           coming from the Spherical Harmonics
//           NormL = pow(2,L+1)*sqrt(2.0/static_cast<real_type>(DFactorial(2*l+1)))*pow(2.0/pi,0.25)
      double L = static_cast<T>(l);
      double NormL = std::pow(2,L+1)*std::sqrt(2.0/static_cast<double>(DFactorial(2*l+1)))*std::pow(2.0/pi,0.25);
      NormFactor[p++] = static_cast<T>(std::pow(2.0/pi,0.75)*std::pow(4.0,0.5*(a+b+c))
          *std::sqrt(1.0/static_cast<double>((DFactorial(2*a-1)*DFactorial(2*b-1)*DFactorial(2*c-1))))/NormL);
    }
  }
}

template<class T>
void SoaCartesianTensor<T>::evaluate_bare(T x, T y, T z, T* restrict XYZ) const
{
  const T x2=x*x, y2=y*y, z2=z*z;
  switch(Lmax)
  {
  case 4:
    XYZ[20]=x2*x2;
    XYZ[21]=y2*y2;
    XYZ[22]=z2*z2;
    XYZ[23]=x2*x*y;
    XYZ[24]=x2*x*z;
    XYZ[25]=y2*y*x;
    XYZ[26]=y2*y*z;
    XYZ[27]=z2*z*x;
    XYZ[28]=z2*z*y;
    XYZ[29]=x2*y2;
    XYZ[30]=x2*z2;
    XYZ[31]=y2*z2;
    XYZ[32]=x2*y*z;
    XYZ[33]=x*y2*z;
    XYZ[34]=x*y*z2;
  case 3:
    XYZ[10]=x2*x;
    XYZ[11]=y2*y;
    XYZ[12]=z2*z;
    XYZ[13]=x2*y;
    XYZ[14]=x2*z;
    XYZ[15]=y2*x;
    XYZ[16]=y2*z;
    XYZ[17]=z2*x;
    XYZ[18]=z2*y;
    XYZ[19]=x*y*z;
  case 2:
    XYZ[4]=x*x;
    XYZ[5]=y*y;
    XYZ[6]=z*z;
    XYZ[7]=x*y;
    XYZ[8]=x*z;
    XYZ[9]=y*z;
  case 1:
    XYZ[1]=x;
    XYZ[2]=y;
    XYZ[3]=z;
  case 0:
    XYZ[0]=1.0;
  }
  //for (int i=0; i<XYZ.size(); i++)
  //  XYZ[i]*= NormFactor[i];
}

template<class T>
void SoaCartesianTensor<T>::evaluateVGL(T x, T y, T z)
{
  CONSTEXPR T czero(0);
  CONSTEXPR T cone(1);
  CONSTEXPR T two(2);
  CONSTEXPR T three(3);
  CONSTEXPR T four(4);
  CONSTEXPR T six(6);
  CONSTEXPR T twelve(12);

  cXYZ=czero;
  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;
  T* restrict XYZ=cXYZ.data(0);
  T* restrict gr0=cXYZ.data(1);
  T* restrict gr1=cXYZ.data(2);
  T* restrict gr2=cXYZ.data(3);
  T* restrict lap=cXYZ.data(4);

  switch(Lmax)
  {
  case 4:
    XYZ[20]= x2*x2;
    gr0[20]= four*x3;
    lap[20]= twelve*x2;
    XYZ[21]= y2*y2;
    gr1[21]= four*y3;
    lap[21]= twelve*y2;
    XYZ[22]= z2*z2;
    gr2[22]= four*z3;
    lap[22]= twelve*z2;
    XYZ[23]= x3*y;
    gr0[23]= three*x2*y;
    gr1[23]= x3;
    lap[23]= six*x*y;
    XYZ[24]= x3*z;
    gr0[24]= three*x2*z;
    gr2[24]= x3;
    lap[24]= six*x*z;
    XYZ[25]= y3*x;
    gr1[25]= three*y2*x;
    gr0[25]= y3;
    lap[25]= six*x*y;
    XYZ[26]= y3*z;
    gr1[26]= three*y2*z;
    gr2[26]= y3;
    lap[26]= six*y*z;
    XYZ[27]= z3*x;
    gr2[27]= three*z2*x;
    gr0[27]= z3;
    lap[27]= six*x*z;
    XYZ[28]= z3*y;
    gr2[28]= three*z2*y;
    gr1[28]= z3;
    lap[28]= six*z*y;
    XYZ[29]= x2*y2;
    gr0[29]= two*x*y2;
    gr1[29]= two*x2*y;
    lap[29]= two*(x2+y2);
    XYZ[30]= x2*z2;
    gr0[30]= two*x*z2;
    gr2[30]= two*x2*z;
    lap[30]= two*(x2+z2);
    XYZ[31]= y2*z2;
    gr1[31]= two*y*z2;
    gr2[31]= two*y2*z;
    lap[31]= two*(y2+z2);
    XYZ[32]= x2*y*z;
    gr0[32]= two*x*y*z;
    gr1[32]= x2*z;
    gr2[32]= x2*y;
    lap[32]= two*y*z;
    XYZ[33]= x*y2*z;
    gr1[33]= two*x*y*z;
    gr0[33]= y2*z;
    gr2[33]= y2*x;
    lap[33]= two*x*z;
    XYZ[34]= x*y*z2;
    gr2[34]= two*x*y*z;
    gr1[34]= z2*x;
    gr0[34]= z2*y;
    lap[34]= two*x*y;
  case 3:
    XYZ[10]= x3;
    gr0[10]= three*x2;
    lap[10]= six*x;
    XYZ[11]= y3;
    gr1[11]= three*y2;
    lap[11]= six*y;
    XYZ[12]= z3;
    gr2[12]= three*z2;
    lap[12]= six*z;
    XYZ[13]= x2*y;
    gr0[13]= two*x*y;
    gr1[13]= x2;
    lap[13]= two*y;
    XYZ[14]= x2*z;
    gr0[14]= two*x*z;
    gr2[14]= x2;
    lap[14]= two*z;
    XYZ[15]= y2*x;
    gr1[15]= two*x*y;
    gr0[15]= y2;
    lap[15]= two*x;
    XYZ[16]= y2*z;
    gr1[16]= two*z*y;
    gr2[16]= y2;
    lap[16]= two*z;
    XYZ[17]= z2*x;
    gr2[17]= two*x*z;
    gr0[17]= z2;
    lap[17]= two*x;
    XYZ[18]= z2*y;
    gr2[18]= two*z*y;
    gr1[18]= z2;
    lap[18]= two*y;
    XYZ[19]= x*y*z;
    gr0[19]= y*z;
    gr1[19]= x*z;
    gr2[19]= x*y;
  case 2:
    XYZ[4]= x*x;
    gr0[4]= 2.0*x;
    lap[4]= 2.0;
    XYZ[5]= y*y;
    gr1[5]= 2.0*y;
    lap[5]= 2.0;
    XYZ[6]= z*z;
    gr2[6]= 2.0*z;
    lap[6]= 2.0;
    XYZ[7]= x*y;
    gr0[7]= y;
    gr1[7]= x;
    XYZ[8]= x*z;
    gr0[8]= z;
    gr2[8]= x;
    XYZ[9]= y*z;
    gr1[9]= z;
    gr2[9]= y;
  case 1:
    XYZ[1]= x;
    gr0[1]= cone;
    XYZ[2]= y;
    gr1[2]= cone;
    XYZ[3]= z;
    gr2[3]= cone;
  case 0:
    XYZ[0]=cone;
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

template<class T>
void SoaCartesianTensor<T>::evaluateVGH(T x, T y, T z)
{
  const T x2=x*x, y2=y*y, z2=z*z;
  const T x3=x2*x, y3=y2*y, z3=z2*z;

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
  case 4:
    XYZ[20] = x2*x2;
    gr0[20] = 4.0*x3;
    h00[20] = 12.0*x2;
    XYZ[21] = y2*y2;
    gr1[21] = 4.0*y3;
    h11[21] = 12.0*y2;
    XYZ[22] = z2*z2;
    gr2[22] = 4.0*z3;
    h22[22] = 12.0*z2;
    XYZ[23] = x3*y;
    gr0[23] = 3.0*x2*y;
    gr1[23] = x3;
    h00[23] = 6.0*x*y;
    h01[23] = 3*x2;
    XYZ[24] = x3*z;
    gr0[24] = 3.0*x2*z;
    gr2[24] = x3;
    h00[24] = 6.0*x*z;
    h02[24] = 3*x2;
    XYZ[25] = y3*x;
    gr1[25] = 3.0*y2*x;
    gr0[25] = y3;
    h11[25] = 6.0*x*y;
    h01[25] = 3*y2;
    XYZ[26] = y3*z;
    gr1[26] = 3.0*y2*z;
    gr2[26] = y3;
    h11[26] = 6.0*z*y;
    h12[26] = 3*y2;
    XYZ[27] = z3*x;
    gr2[27] = 3.0*z2*x;
    gr0[27] = z3;
    h11[27] = 6.0*z*x;
    h02[27] = 3*z2;
    XYZ[28] = z3*y;
    gr2[28] = 3.0*z2*y;
    gr1[28] = z3;
    h11[28] = 6.0*z*y;
    h12[28] = 3*z2;
    XYZ[29] = x2*y2;
    gr0[29] = 2.0*x*y2;
    gr1[29] = 2.0*x2*y;
    h00[29] = 2.0*y2;
    h11[29] = 2.0*x2;
    h01[29] = 4.0*x*y;
    XYZ[30] = x2*z2;
    gr0[30] = 2.0*x*z2;
    gr2[30] = 2.0*x2*z;
    h00[30] = 2.0*z2;
    h22[30] = 2.0*x2;
    h02[30] = 4.0*x*z;
    XYZ[31] = y2*z2;
    gr1[31] = 2.0*y*z2;
    gr2[31] = 2.0*y2*z;
    h22[31] = 2.0*y2;
    h11[31] = 2.0*z2;
    h12[31] = 4.0*z*y;
    XYZ[32] = x2*y*z;
    gr0[32] = 2.0*x*y*z;
    gr1[32] = x2*z;
    gr2[32] = x2*y;
    h00[32] = 2.0*y*z;
    h01[32] = 2.0*x*z;
    h02[32] = 2.0*x*y;
    h12[32] = x*x;
    XYZ[33]= x*y2*z;
    gr1[33] = 2.0*x*y*z;
    gr0[33] = y2*z;
    gr2[33] = y2*x;
    h11[33](1,1)= 2.0*x*z;
    h01[33] = 2.0*y*z;
    h02[33] = y*y;
    h12[33] = 2.0*x*y;
    XYZ[34] = x*y*z2;
    gr2[34] = 2.0*x*y*z;
    gr1[34] = z2*x;
    gr0[34] = z2*y;
    h22[34] = 2.0*y*x;
    h01[34] = z*z;
    h02[34] = 2.0*z*y;
    h12[34] = 2.0*z*x;
  case 3:
    XYZ[10]= x3;
    gr0[10]= 3.0*x2;
    h00[10]= 6.0*x;
    XYZ[11]= y3;
    gr1[11]= 3.0*y2;
    h11[11]= 6.0*y;
    XYZ[12]= z3;
    gr2[12]= 3.0*z2;
    h22[12]= 6.0*z;
    XYZ[13]= x2*y;
    gr0[13]= 2.0*x*y;
    gr1[13]= x2;
    h00[13]= 2.0*y;
    h01[13]= 2.0*x;
    XYZ[14]= x2*z;
    gr0[14]= 2.0*x*z;
    gr2[14]= x2;
    h00[14]= 2.0*z;
    h02[14]= 2.0*x;
    XYZ[15]= y2*x;
    gr1[15]= 2.0*x*y;
    gr0[15]= y2;
    h11[15]= 2.0*x;
    h01[15]= 2.0*y;
    XYZ[16]= y2*z;
    gr1[16]= 2.0*z*y;
    gr2[16]= y2;
    h11[16]= 2.0*z;
    h12[16]= 2.0*y;
    XYZ[17]= z2*x;
    gr2[17]= 2.0*x*z;
    gr0[17]= z2;
    h22[17]= 2.0*x;
    h02[17]= 2.0*z;
    XYZ[18]= z2*y;
    gr2[18]= 2.0*z*y;
    gr1[18]= z2;
    h22[18]= 2.0*y;
    h12[18]= 2.0*z;
    XYZ[19]= x*y*z;
    gr0[19]= y*z;
    gr1[19]= x*z;
    gr2[19]= x*y;
    h01[19]= z;
    h02[19]= y;
    h12[19]= x;
  case 2:
    XYZ[4]= x*x;
    gr0[4]= 2.0*x;
    h00[4]= 2.0;
    XYZ[5]= y*y;
    gr1[5]= 2.0*y;
    h11[5]= 2.0;
    XYZ[6]= z*z;
    gr2[6]= 2.0*z;
    h22[6]= 2.0;
    XYZ[7]= x*y;
    gr0[7]= y;
    gr1[7]= x;
    h01[7]= 1.0;
    XYZ[8]= x*z;
    gr0[8]= z;
    gr2[8]= x;
    h02[8]= 1.0;
    XYZ[9]= y*z;
    gr1[9]= z;
    gr2[9]= y;
    h12[9]= 1.0;
  case 1:
    XYZ[1]=x;
    gr0[1][0]= 1.0;
    XYZ[2]=y;
    gr1[2][1]= 1.0;
    XYZ[3]=z;
    gr2[3][2]= 1.0;
  case 0:
    XYZ[0]=1.0;
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

#if 0
template<class T>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateWithThirdDeriv(const Point_t& p)
{
  const value_type x=p[0], y=p[1], z=p[2];
  evaluateVGL(x,y,z);
  evaluateThirdDerivOnly(p);
  const value_type x2=x*x, y2=y*y, z2=z*z;
  const value_type x3=x2*x, y3=y2*y, z3=z2*z;

}

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
  case 4:
 ////// Autogenerated cartesian 3rd derivative tensor components /////
    (gggXYZ[20][0])(0,0) = 24.0*x;
 
    (gggXYZ[21][1])(1,1) = 24.0*y;
 
    (gggXYZ[22][2])(2,2) = 24.0*z;
 
    (gggXYZ[23][0])(0,0) = 6.0*y;
    (gggXYZ[23][1])(0,0) = (gggXYZ[23][0])(1,0) = (gggXYZ[23][0])(0,1) = 6.0*x;
 
    (gggXYZ[24][0])(0,0) = 6.0*z;
    (gggXYZ[24][0])(2,0) = (gggXYZ[24][0])(0,2) = (gggXYZ[24][2])(0,0) = 6.0*x;
 
    (gggXYZ[25][0])(1,1) = (gggXYZ[25][1])(1,0) = (gggXYZ[25][1])(0,1) = 6.0*y;
    (gggXYZ[25][1])(1,1) = 6.0*x;
 
    (gggXYZ[26][1])(1,1) = 6.0*z;
    (gggXYZ[26][1])(2,1) = (gggXYZ[26][2])(1,1) = (gggXYZ[26][1])(1,2) = 6.0*y;
 
    (gggXYZ[27][2])(0,2) = (gggXYZ[27][2])(2,0) = (gggXYZ[27][0])(2,2) = 6.0*z;
    (gggXYZ[27][2])(2,2) = 6.0*x;
 
    (gggXYZ[28][1])(2,2) = (gggXYZ[28][2])(2,1) = (gggXYZ[28][2])(1,2) = 6.0*z;
    (gggXYZ[28][2])(2,2) = 6.0*y;
 
    (gggXYZ[29][1])(0,0) = (gggXYZ[29][0])(1,0) = (gggXYZ[29][0])(0,1) = 4.0*y;
    (gggXYZ[29][0])(1,1) = (gggXYZ[29][1])(1,0) = (gggXYZ[29][1])(0,1) = 4.0*x;
 
    (gggXYZ[30][0])(2,0) = (gggXYZ[30][0])(0,2) = (gggXYZ[30][2])(0,0) = 4.0*z;
    (gggXYZ[30][2])(0,2) = (gggXYZ[30][2])(2,0) = (gggXYZ[30][0])(2,2) = 4.0*x;
 
    (gggXYZ[31][1])(2,1) = (gggXYZ[31][2])(1,1) = (gggXYZ[31][1])(1,2) = 4.0*z;
    (gggXYZ[31][1])(2,2) = (gggXYZ[31][2])(2,1) = (gggXYZ[31][2])(1,2) = 4.0*y;
 
    (gggXYZ[32][1])(0,0) = (gggXYZ[32][0])(1,0) = (gggXYZ[32][0])(0,1) = 2.0*z;
    (gggXYZ[32][0])(2,0) = (gggXYZ[32][0])(0,2) = (gggXYZ[32][2])(0,0) = 2.0*y;
    (gggXYZ[32][2])(1,0) = (gggXYZ[32][0])(1,2) = (gggXYZ[32][1])(0,2) = (gggXYZ[32][2])(0,1) = (gggXYZ[32][0])(2,1) = (gggXYZ[32][1])(2,0) = 2.0*x;
 
    (gggXYZ[33][0])(1,1) = (gggXYZ[33][1])(1,0) = (gggXYZ[33][1])(0,1) = 2.0*z;
    (gggXYZ[33][2])(1,0) = (gggXYZ[33][0])(1,2) = (gggXYZ[33][1])(0,2) = (gggXYZ[33][2])(0,1) = (gggXYZ[33][0])(2,1) = (gggXYZ[33][1])(2,0) = 2.0*y;
    (gggXYZ[33][1])(2,1) = (gggXYZ[33][2])(1,1) = (gggXYZ[33][1])(1,2) = 2.0*x;
 
    (gggXYZ[34][2])(1,0) = (gggXYZ[34][0])(1,2) = (gggXYZ[34][1])(0,2) = (gggXYZ[34][2])(0,1) = (gggXYZ[34][0])(2,1) = (gggXYZ[34][1])(2,0) = 2.0*z;
    (gggXYZ[34][2])(0,2) = (gggXYZ[34][2])(2,0) = (gggXYZ[34][0])(2,2) = 2.0*y;
    (gggXYZ[34][1])(2,2) = (gggXYZ[34][2])(2,1) = (gggXYZ[34][2])(1,2) = 2.0*x;

  case 3:
 ////// Autogenerated cartesian 3rd derivative tensor components /////
    (gggXYZ[10][0])(0,0) = 6.0;
 
    (gggXYZ[11][1])(1,1) = 6.0;
 
    (gggXYZ[12][2])(2,2) = 6.0;
 
    (gggXYZ[13][1])(0,0) = (gggXYZ[13][0])(1,0) = (gggXYZ[13][0])(0,1) = 2.0;
 
    (gggXYZ[14][0])(2,0) = (gggXYZ[14][0])(0,2) = (gggXYZ[14][2])(0,0) = 2.0;
 
    (gggXYZ[15][0])(1,1) = (gggXYZ[15][1])(1,0) = (gggXYZ[15][1])(0,1) = 2.0;
 
    (gggXYZ[16][1])(2,1) = (gggXYZ[16][2])(1,1) = (gggXYZ[16][1])(1,2) = 2.0;
 
    (gggXYZ[17][2])(0,2) = (gggXYZ[17][2])(2,0) = (gggXYZ[17][0])(2,2) = 2.0;
 
    (gggXYZ[18][1])(2,2) = (gggXYZ[18][2])(2,1) = (gggXYZ[18][2])(1,2) = 2.0;
 
    (gggXYZ[19][2])(1,0) = (gggXYZ[19][0])(1,2) = (gggXYZ[19][1])(0,2) = (gggXYZ[19][2])(0,1) = (gggXYZ[19][0])(2,1) = (gggXYZ[19][1])(2,0) = 1.0;
 
  }
  for (int i=0; i<ntot; i++)
  {
    gggXYZ[i][0] *= NormFactor[i];
    gggXYZ[i][1] *= NormFactor[i];
    gggXYZ[i][2] *= NormFactor[i];
  }
}
#endif

template<class T>
void SoaCartesianTensor<T>::getABC(int n, int& a, int& b, int& c)
{
// following Gamess notation
  switch(n)
  {
//S
  case 0: // S
    a=b=c=0;
    break;
//P
  case 1: // X
    a=1;
    b=c=0;
    break;
  case 2: // Y
    b=1;
    a=c=0;
    break;
  case 3: // Z
    c=1;
    a=b=0;
    break;
// D
  case 4: // XX
    a=2;
    b=c=0;
    break;
  case 5: // YY
    b=2;
    a=c=0;
    break;
  case 6: // ZZ
    c=2;
    a=b=0;
    break;
  case 7: // XY
    a=b=1;
    c=0;
    break;
  case 8: // XZ
    a=c=1;
    b=0;
    break;
  case 9: // YZ
    c=b=1;
    a=0;
    break;
//F
  case 10: // XXX
    a=3;
    b=c=0;
    break;
  case 11: // YYY
    b=3;
    a=c=0;
    break;
  case 12: // ZZZ
    c=3;
    a=b=0;
    break;
  case 13: // XXY
    a=2;
    b=1;
    c=0;
    break;
  case 14: // XXZ
    a=2;
    c=1;
    b=0;
    break;
  case 15: // YYX
    b=2;
    a=1;
    c=0;
    break;
  case 16: // YYZ
    b=2;
    c=1;
    a=0;
    break;
  case 17: // ZZX
    c=2;
    b=0;
    a=1;
    break;
  case 18: // ZZY
    c=2;
    b=1;
    a=0;
    break;
  case 19: // XYZ
    a=b=c=1;
    break;
// G
  case 20: // XXXX
    a=4;
    b=0;
    c=0;
    break;
  case 21: // YYYY
    a=0;
    b=4;
    c=0;
    break;
  case 22: // ZZZZ
    a=0;
    b=0;
    c=4;
    break;
  case 23: // XXXY
    a=3;
    b=1;
    c=0;
    break;
  case 24: // XXXZ
    a=3;
    b=0;
    c=1;
    break;
  case 25: // YYYX
    a=1;
    b=3;
    c=0;
    break;
  case 26: // YYYZ
    a=0;
    b=3;
    c=1;
    break;
  case 27: // ZZZX
    a=1;
    b=0;
    c=3;
    break;
  case 28: // ZZZY
    a=0;
    b=1;
    c=3;
    break;
  case 29: // XXYY
    a=2;
    b=2;
    c=0;
    break;
  case 30: // XXZZ
    a=2;
    b=0;
    c=2;
    break;
  case 31: // YYZZ
    a=0;
    b=2;
    c=2;
    break;
  case 32: // XXYZ
    a=2;
    b=1;
    c=1;
    break;
  case 33: // YYXZ
    a=1;
    b=2;
    c=1;
    break;
  case 34: // ZZXY
    a=1;
    b=1;
    c=2;
    break;
// no H or higher in Gamess
  default:
    std::cerr <<"CartesianTensor::getABC() - Incorrect index.\n";
    APP_ABORT("");
    break;
  }
}

} //namespace qmcplusplus
#endif
