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
template<class T, class Point_t, class Tensor_t = qmcplusplus::Tensor<T,3> , class GGG_t = qmcplusplus::TinyVector<Tensor_t, 3> >
class CartesianTensor
{
public :

  typedef T value_type;
  typedef Point_t pos_type;
  typedef Tensor_t hess_type;
  typedef GGG_t ggg_type;
  typedef CartesianTensor<T,Point_t,Tensor_t> This_t;

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

  inline value_type getYlm(int lm) const
  {
    return XYZ[lm];
  }

  inline Point_t getGradYlm(int lm) const
  {
    return gradXYZ[lm];
  }

  inline value_type getLaplYlm(int lm) const
  {
    return laplXYZ[lm];
  }

  inline Tensor_t getHessYlm(int lm) const
  {
    return hessXYZ[lm];
  }

  inline GGG_t getGGGYlm(int lm) const
  {
    return gggXYZ[lm];
  }

  inline int size() const
  {
    return XYZ.size();
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
  if(Lmax < 0 || Lmax > 4)
  {
    std::cerr <<"CartesianTensor can't handle Lmax > 4 or Lmax < 0.\n";
    APP_ABORT("");
  }
  int ntot = 0;
  for(int i=0; i<=Lmax; i++)
    ntot+=(i+1)*(i+2)/2;
  XYZ.resize(ntot);
  gradXYZ.resize(ntot);
  laplXYZ.resize(ntot);
  hessXYZ.resize(ntot);
  gggXYZ.resize(ntot);
  NormFactor.resize(ntot,1);
  int p=0;
  int a,b,c;
  const double pi = 4.0*atan(1.0);
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
      double L = static_cast<double>(l);
      double NormL = pow(2,L+1)*sqrt(2.0/static_cast<double>(DFactorial(2*l+1)))*pow(2.0/pi,0.25);
      NormFactor[p++] = pow(2.0/pi,0.75)*pow(4.0,0.5*(a+b+c))*std::sqrt(1.0/static_cast<double>((DFactorial(2*a-1)*DFactorial(2*b-1)*DFactorial(2*c-1))))/NormL;
    }
  }
}

template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluate(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
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
  for (int i=0; i<XYZ.size(); i++)
    XYZ[i]*= NormFactor[i];
}

template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateAll(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
    gradXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
    laplXYZ[i]=0.0;
  switch(Lmax)
  {
  case 4:
    XYZ[20]=x2*x2;
    gradXYZ[20][0]= 4.0*x3;
    laplXYZ[20]= 12.0*x2;
    XYZ[21]=y2*y2;
    gradXYZ[21][1]= 4.0*y3;
    laplXYZ[21]= 12.0*y2;
    XYZ[22]=z2*z2;
    gradXYZ[22][2]= 4.0*z3;
    laplXYZ[22]= 12.0*z2;
    XYZ[23]=x3*y;
    gradXYZ[23][0]= 3.0*x2*y;
    gradXYZ[23][1]= x3;
    laplXYZ[23]= 6.0*x*y;
    XYZ[24]=x3*z;
    gradXYZ[24][0]= 3.0*x2*z;
    gradXYZ[24][2]= x3;
    laplXYZ[24]= 6.0*x*z;
    XYZ[25]=y3*x;
    gradXYZ[25][1]= 3.0*y2*x;
    gradXYZ[25][0]= y3;
    laplXYZ[25]= 6.0*x*y;
    XYZ[26]=y3*z;
    gradXYZ[26][1]= 3.0*y2*z;
    gradXYZ[26][2]= y3;
    laplXYZ[26]= 6.0*y*z;
    XYZ[27]=z3*x;
    gradXYZ[27][2]= 3.0*z2*x;
    gradXYZ[27][0]= z3;
    laplXYZ[27]= 6.0*x*z;
    XYZ[28]=z3*y;
    gradXYZ[28][2]= 3.0*z2*y;
    gradXYZ[28][1]= z3;
    laplXYZ[28]= 6.0*z*y;
    XYZ[29]=x2*y2;
    gradXYZ[29][0]= 2.0*x*y2;
    gradXYZ[29][1]= 2.0*x2*y;
    laplXYZ[29]= 2.0*x2+2.0*y2;
    XYZ[30]=x2*z2;
    gradXYZ[30][0]= 2.0*x*z2;
    gradXYZ[30][2]= 2.0*x2*z;
    laplXYZ[30]= 2.0*x2+2.0*z2;
    XYZ[31]=y2*z2;
    gradXYZ[31][1]= 2.0*y*z2;
    gradXYZ[31][2]= 2.0*y2*z;
    laplXYZ[31]= 2.0*y2+2.0*z2;
    XYZ[32]=x2*y*z;
    gradXYZ[32][0]= 2.0*x*y*z;
    gradXYZ[32][1]= x2*z;
    gradXYZ[32][2]= x2*y;
    laplXYZ[32]= 2.0*y*z;
    XYZ[33]=x*y2*z;
    gradXYZ[33][1]= 2.0*x*y*z;
    gradXYZ[33][0]= y2*z;
    gradXYZ[33][2]= y2*x;
    laplXYZ[33]= 2.0*x*z;
    XYZ[34]=x*y*z2;
    gradXYZ[34][2]= 2.0*x*y*z;
    gradXYZ[34][1]= z2*x;
    gradXYZ[34][0]= z2*y;
    laplXYZ[34]= 2.0*x*y;
  case 3:
    XYZ[10]=x3;
    gradXYZ[10][0]= 3.0*x2;
    laplXYZ[10]= 6.0*x;
    XYZ[11]=y3;
    gradXYZ[11][1]= 3.0*y2;
    laplXYZ[11]= 6.0*y;
    XYZ[12]=z3;
    gradXYZ[12][2]= 3.0*z2;
    laplXYZ[12]= 6.0*z;
    XYZ[13]=x2*y;
    gradXYZ[13][0]= 2.0*x*y;
    gradXYZ[13][1]= x2;
    laplXYZ[13]= 2.0*y;
    XYZ[14]=x2*z;
    gradXYZ[14][0]= 2.0*x*z;
    gradXYZ[14][2]= x2;
    laplXYZ[14]= 2.0*z;
    XYZ[15]=y2*x;
    gradXYZ[15][1]= 2.0*x*y;
    gradXYZ[15][0]= y2;
    laplXYZ[15]= 2.0*x;
    XYZ[16]=y2*z;
    gradXYZ[16][1]= 2.0*z*y;
    gradXYZ[16][2]= y2;
    laplXYZ[16]= 2.0*z;
    XYZ[17]=z2*x;
    gradXYZ[17][2]= 2.0*x*z;
    gradXYZ[17][0]= z2;
    laplXYZ[17]= 2.0*x;
    XYZ[18]=z2*y;
    gradXYZ[18][2]= 2.0*z*y;
    gradXYZ[18][1]= z2;
    laplXYZ[18]= 2.0*y;
    XYZ[19]=x*y*z;
    gradXYZ[19][0]= y*z;
    gradXYZ[19][1]= x*z;
    gradXYZ[19][2]= x*y;
  case 2:
    XYZ[4]=x*x;
    gradXYZ[4][0]= 2.0*x;
    laplXYZ[4]= 2.0;
    XYZ[5]=y*y;
    gradXYZ[5][1]= 2.0*y;
    laplXYZ[5]= 2.0;
    XYZ[6]=z*z;
    gradXYZ[6][2]= 2.0*z;
    laplXYZ[6]= 2.0;
    XYZ[7]=x*y;
    gradXYZ[7][0]= y;
    gradXYZ[7][1]= x;
    XYZ[8]=x*z;
    gradXYZ[8][0]= z;
    gradXYZ[8][2]= x;
    XYZ[9]=y*z;
    gradXYZ[9][1]= z;
    gradXYZ[9][2]= y;
  case 1:
    XYZ[1]=x;
    gradXYZ[1][0]= 1.0;
    XYZ[2]=y;
    gradXYZ[2][1]= 1.0;
    XYZ[3]=z;
    gradXYZ[3][2]= 1.0;
  case 0:
    XYZ[0]=1.0;
  }
  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    laplXYZ[i]*= NormFactor[i];
}

template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateWithHessian(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
  int ntot=XYZ.size();
  for (int i=0; i<ntot; i++)
    gradXYZ[i]=0.0;
  for (int i=0; i<ntot; i++)
    hessXYZ[i]=0.0;
  switch(Lmax)
  {
  case 4:
    XYZ[20]=x2*x2;
    gradXYZ[20][0]= 4.0*x3;
    hessXYZ[20](0,0)= 12.0*x2;
    XYZ[21]=y2*y2;
    gradXYZ[21][1]= 4.0*y3;
    hessXYZ[21](1,1)= 12.0*y2;
    XYZ[22]=z2*z2;
    gradXYZ[22][2]= 4.0*z3;
    hessXYZ[22](2,2)= 12.0*z2;
    XYZ[23]=x3*y;
    gradXYZ[23][0]= 3.0*x2*y;
    gradXYZ[23][1]= x3;
    hessXYZ[23](0,0)= 6.0*x*y;
    hessXYZ[23](0,1)= hessXYZ[23](1,0) = 3*x2;
    XYZ[24]=x3*z;
    gradXYZ[24][0]= 3.0*x2*z;
    gradXYZ[24][2]= x3;
    hessXYZ[24](0,0)= 6.0*x*z;
    hessXYZ[24](0,2)= hessXYZ[24](2,0) = 3*x2;
    XYZ[25]=y3*x;
    gradXYZ[25][1]= 3.0*y2*x;
    gradXYZ[25][0]= y3;
    hessXYZ[25](1,1)= 6.0*x*y;
    hessXYZ[25](0,1)= hessXYZ[25](1,0) = 3*y2;
    XYZ[26]=y3*z;
    gradXYZ[26][1]= 3.0*y2*z;
    gradXYZ[26][2]= y3;
    hessXYZ[26](1,1)= 6.0*z*y;
    hessXYZ[26](2,1)= hessXYZ[26](1,2) = 3*y2;
    XYZ[27]=z3*x;
    gradXYZ[27][2]= 3.0*z2*x;
    gradXYZ[27][0]= z3;
    hessXYZ[27](2,2)= 6.0*z*x;
    hessXYZ[27](2,0)= hessXYZ[27](0,2) = 3*z2;
    XYZ[28]=z3*y;
    gradXYZ[28][2]= 3.0*z2*y;
    gradXYZ[28][1]= z3;
    hessXYZ[28](2,2)= 6.0*z*y;
    hessXYZ[28](2,1)= hessXYZ[28](1,2) = 3*z2;
    XYZ[29]=x2*y2;
    gradXYZ[29][0]= 2.0*x*y2;
    gradXYZ[29][1]= 2.0*x2*y;
    hessXYZ[29](0,0)= 2.0*y2;
    hessXYZ[29](1,1)= 2.0*x2;
    hessXYZ[29](0,1)= hessXYZ[29](1,0) = 4.0*x*y;
    XYZ[30]=x2*z2;
    gradXYZ[30][0]= 2.0*x*z2;
    gradXYZ[30][2]= 2.0*x2*z;
    hessXYZ[30](0,0)= 2.0*z2;
    hessXYZ[30](2,2)= 2.0*x2;
    hessXYZ[30](0,2)= hessXYZ[30](2,0) = 4.0*x*z;
    XYZ[31]=y2*z2;
    gradXYZ[31][1]= 2.0*y*z2;
    gradXYZ[31][2]= 2.0*y2*z;
    hessXYZ[31](2,2)= 2.0*y2;
    hessXYZ[31](1,1)= 2.0*z2;
    hessXYZ[31](1,2)= hessXYZ[31](2,1) = 4.0*z*y;
    XYZ[32]=x2*y*z;
    gradXYZ[32][0]= 2.0*x*y*z;
    gradXYZ[32][1]= x2*z;
    gradXYZ[32][2]= x2*y;
    hessXYZ[32](0,0)= 2.0*y*z;
    hessXYZ[32](0,1)= hessXYZ[32](1,0) = 2.0*x*z;
    hessXYZ[32](0,2)= hessXYZ[32](2,0) = 2.0*x*y;
    hessXYZ[32](1,2)= hessXYZ[32](2,1) = x*x;
    XYZ[33]=x*y2*z;
    gradXYZ[33][1]= 2.0*x*y*z;
    gradXYZ[33][0]= y2*z;
    gradXYZ[33][2]= y2*x;
    hessXYZ[33](1,1)= 2.0*x*z;
    hessXYZ[33](0,1)= hessXYZ[33](1,0) = 2.0*y*z;
    hessXYZ[33](0,2)= hessXYZ[33](2,0) = y*y;
    hessXYZ[33](1,2)= hessXYZ[33](2,1) = 2.0*x*y;
    XYZ[34]=x*y*z2;
    gradXYZ[34][2]= 2.0*x*y*z;
    gradXYZ[34][1]= z2*x;
    gradXYZ[34][0]= z2*y;
    hessXYZ[34](2,2)= 2.0*y*x;
    hessXYZ[34](0,1)= hessXYZ[34](1,0) = z*z;
    hessXYZ[34](0,2)= hessXYZ[34](2,0) = 2.0*z*y;
    hessXYZ[34](1,2)= hessXYZ[34](2,1) = 2.0*z*x;
  case 3:
    XYZ[10]=x3;
    gradXYZ[10][0]= 3.0*x2;
    hessXYZ[10](0,0)= 6.0*x;
    XYZ[11]=y3;
    gradXYZ[11][1]= 3.0*y2;
    hessXYZ[11](1,1)= 6.0*y;
    XYZ[12]=z3;
    gradXYZ[12][2]= 3.0*z2;
    hessXYZ[12](2,2)= 6.0*z;
    XYZ[13]=x2*y;
    gradXYZ[13][0]= 2.0*x*y;
    gradXYZ[13][1]= x2;
    hessXYZ[13](0,0)= 2.0*y;
    hessXYZ[13](0,1)=hessXYZ[13](1,0)= 2.0*x;
    XYZ[14]=x2*z;
    gradXYZ[14][0]= 2.0*x*z;
    gradXYZ[14][2]= x2;
    hessXYZ[14](0,0)= 2.0*z;
    hessXYZ[14](0,2)=hessXYZ[14](2,0)= 2.0*x;
    XYZ[15]=y2*x;
    gradXYZ[15][1]= 2.0*x*y;
    gradXYZ[15][0]= y2;
    hessXYZ[15](1,1)= 2.0*x;
    hessXYZ[15](0,1)=hessXYZ[15](1,0)= 2.0*y;
    XYZ[16]=y2*z;
    gradXYZ[16][1]= 2.0*z*y;
    gradXYZ[16][2]= y2;
    hessXYZ[16](1,1)= 2.0*z;
    hessXYZ[16](2,1)=hessXYZ[16](1,2)= 2.0*y;
    XYZ[17]=z2*x;
    gradXYZ[17][2]= 2.0*x*z;
    gradXYZ[17][0]= z2;
    hessXYZ[17](2,2)= 2.0*x;
    hessXYZ[17](0,2)=hessXYZ[17](2,0)= 2.0*z;
    XYZ[18]=z2*y;
    gradXYZ[18][2]= 2.0*z*y;
    gradXYZ[18][1]= z2;
    hessXYZ[18](2,2)= 2.0*y;
    hessXYZ[18](2,1)=hessXYZ[18](1,2)= 2.0*z;
    XYZ[19]=x*y*z;
    gradXYZ[19][0]= y*z;
    gradXYZ[19][1]= x*z;
    gradXYZ[19][2]= x*y;
    hessXYZ[19](0,1)=hessXYZ[19](1,0)=z;
    hessXYZ[19](0,2)=hessXYZ[19](2,0)=y;
    hessXYZ[19](2,1)=hessXYZ[19](1,2)=x;
  case 2:
    XYZ[4]=x*x;
    gradXYZ[4][0]= 2.0*x;
    hessXYZ[4](0,0)= 2.0;
    XYZ[5]=y*y;
    gradXYZ[5][1]= 2.0*y;
    hessXYZ[5](1,1)= 2.0;
    XYZ[6]=z*z;
    gradXYZ[6][2]= 2.0*z;
    hessXYZ[6](2,2)= 2.0;
    XYZ[7]=x*y;
    gradXYZ[7][0]= y;
    gradXYZ[7][1]= x;
    hessXYZ[7](0,1)=hessXYZ[7](1,0)= 1.0;
    XYZ[8]=x*z;
    gradXYZ[8][0]= z;
    gradXYZ[8][2]= x;
    hessXYZ[8](0,2)=hessXYZ[8](2,0)= 1.0;
    XYZ[9]=y*z;
    gradXYZ[9][1]= z;
    gradXYZ[9][2]= y;
    hessXYZ[9](2,1)=hessXYZ[9](1,2)= 1.0;
  case 1:
    XYZ[1]=x;
    gradXYZ[1][0]= 1.0;
    XYZ[2]=y;
    gradXYZ[2][1]= 1.0;
    XYZ[3]=z;
    gradXYZ[3][2]= 1.0;
  case 0:
    XYZ[0]=1.0;
  }
  for (int i=0; i<ntot; i++)
    XYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    gradXYZ[i]*= NormFactor[i];
  for (int i=0; i<ntot; i++)
    hessXYZ[i]*= NormFactor[i];
}

template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::evaluateWithThirdDeriv(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  value_type x2=x*x, y2=y*y, z2=z*z;
  value_type x3=x2*x, y3=y2*y, z3=z2*z;
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
  case 4:
    XYZ[20]=x2*x2;
    gradXYZ[20][0]= 4.0*x3;
    hessXYZ[20](0,0)= 12.0*x2;
    XYZ[21]=y2*y2;
    gradXYZ[21][1]= 4.0*y3;
    hessXYZ[21](1,1)= 12.0*y2;
    XYZ[22]=z2*z2;
    gradXYZ[22][2]= 4.0*z3;
    hessXYZ[22](2,2)= 12.0*z2;
    XYZ[23]=x3*y;
    gradXYZ[23][0]= 3.0*x2*y;
    gradXYZ[23][1]= x3;
    hessXYZ[23](0,0)= 6.0*x*y;
    hessXYZ[23](0,1)= hessXYZ[23](1,0) = 3.0*x2;
    XYZ[24]=x3*z;
    gradXYZ[24][0]= 3.0*x2*z;
    gradXYZ[24][2]= x3;
    hessXYZ[24](0,0)= 6.0*x*z;
    hessXYZ[24](0,2)= hessXYZ[24](2,0) = 3.0*x2;
    XYZ[25]=y3*x;
    gradXYZ[25][1]= 3.0*y2*x;
    gradXYZ[25][0]= y3;
    hessXYZ[25](1,1)= 6.0*x*y;
    hessXYZ[25](0,1)= hessXYZ[25](1,0) = 3.0*y2;
    XYZ[26]=y3*z;
    gradXYZ[26][1]= 3.0*y2*z;
    gradXYZ[26][2]= y3;
    hessXYZ[26](1,1)= 6.0*z*y;
    hessXYZ[26](2,1)= hessXYZ[26](1,2) = 3.0*y2;
    XYZ[27]=z3*x;
    gradXYZ[27][2]= 3.0*z2*x;
    gradXYZ[27][0]= z3;
    hessXYZ[27](2,2)= 6.0*z*x;
    hessXYZ[27](2,0)= hessXYZ[27](0,2) = 3.0*z2;
    XYZ[28]=z3*y;
    gradXYZ[28][2]= 3.0*z2*y;
    gradXYZ[28][1]= z3;
    hessXYZ[28](2,2)= 6.0*z*y;
    hessXYZ[28](2,1)= hessXYZ[28](1,2) = 3.0*z2;
    XYZ[29]=x2*y2;
    gradXYZ[29][0]= 2.0*x*y2;
    gradXYZ[29][1]= 2.0*x2*y;
    hessXYZ[29](0,0)= 2.0*y2;
    hessXYZ[29](1,1)= 2.0*x2;
    hessXYZ[29](0,1)= hessXYZ[29](1,0) = 4.0*x*y;
    XYZ[30]=x2*z2;
    gradXYZ[30][0]= 2.0*x*z2;
    gradXYZ[30][2]= 2.0*x2*z;
    hessXYZ[30](0,0)= 2.0*z2;
    hessXYZ[30](2,2)= 2.0*x2;
    hessXYZ[30](0,2)= hessXYZ[30](2,0) = 4.0*x*z;
    XYZ[31]=y2*z2;
    gradXYZ[31][1]= 2.0*y*z2;
    gradXYZ[31][2]= 2.0*y2*z;
    hessXYZ[31](2,2)= 2.0*y2;
    hessXYZ[31](1,1)= 2.0*z2;
    hessXYZ[31](1,2)= hessXYZ[31](2,1) = 4.0*z*y;
    XYZ[32]=x2*y*z;
    gradXYZ[32][0]= 2.0*x*y*z;
    gradXYZ[32][1]= x2*z;
    gradXYZ[32][2]= x2*y;
    hessXYZ[32](0,0)= 2.0*y*z;
    hessXYZ[32](0,1)= hessXYZ[32](1,0) = 2.0*x*z;
    hessXYZ[32](0,2)= hessXYZ[32](2,0) = 2.0*x*y;
    hessXYZ[32](1,2)= hessXYZ[32](2,1) = x*x;
    XYZ[33]=x*y2*z;
    gradXYZ[33][1]= 2.0*x*y*z;
    gradXYZ[33][0]= y2*z;
    gradXYZ[33][2]= y2*x;
    hessXYZ[33](1,1)= 2.0*x*z;
    hessXYZ[33](0,1)= hessXYZ[33](1,0) = 2.0*y*z;
    hessXYZ[33](0,2)= hessXYZ[33](2,0) = y*y;
    hessXYZ[33](1,2)= hessXYZ[33](2,1) = 2.0*x*y;
    XYZ[34]=x*y*z2;
    gradXYZ[34][2]= 2.0*x*y*z;
    gradXYZ[34][1]= z2*x;
    gradXYZ[34][0]= z2*y;
    hessXYZ[34](2,2)= 2.0*y*x;
    hessXYZ[34](0,1)= hessXYZ[34](1,0) = z*z;
    hessXYZ[34](0,2)= hessXYZ[34](2,0) = 2.0*z*y;
    hessXYZ[34](1,2)= hessXYZ[34](2,1) = 2.0*z*x;
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
    XYZ[10]=x3;
    gradXYZ[10][0]= 3.0*x2;
    hessXYZ[10](0,0)= 6.0*x;
    XYZ[11]=y3;
    gradXYZ[11][1]= 3.0*y2;
    hessXYZ[11](1,1)= 6.0*y;
    XYZ[12]=z3;
    gradXYZ[12][2]= 3.0*z2;
    hessXYZ[12](2,2)= 6.0*z;
    XYZ[13]=x2*y;
    gradXYZ[13][0]= 2.0*x*y;
    gradXYZ[13][1]= x2;
    hessXYZ[13](0,0)= 2.0*y;
    hessXYZ[13](0,1)=hessXYZ[13](1,0)= 2.0*x;
    XYZ[14]=x2*z;
    gradXYZ[14][0]= 2.0*x*z;
    gradXYZ[14][2]= x2;
    hessXYZ[14](0,0)= 2.0*z;
    hessXYZ[14](0,2)=hessXYZ[14](2,0)= 2.0*x;
    XYZ[15]=y2*x;
    gradXYZ[15][1]= 2.0*x*y;
    gradXYZ[15][0]= y2;
    hessXYZ[15](1,1)= 2.0*x;
    hessXYZ[15](0,1)=hessXYZ[15](1,0)= 2.0*y;
    XYZ[16]=y2*z;
    gradXYZ[16][1]= 2.0*z*y;
    gradXYZ[16][2]= y2;
    hessXYZ[16](1,1)= 2.0*z;
    hessXYZ[16](2,1)=hessXYZ[16](1,2)= 2.0*y;
    XYZ[17]=z2*x;
    gradXYZ[17][2]= 2.0*x*z;
    gradXYZ[17][0]= z2;
    hessXYZ[17](2,2)= 2.0*x;
    hessXYZ[17](0,2)=hessXYZ[17](2,0)= 2.0*z;
    XYZ[18]=z2*y;
    gradXYZ[18][2]= 2.0*z*y;
    gradXYZ[18][1]= z2;
    hessXYZ[18](2,2)= 2.0*y;
    hessXYZ[18](2,1)=hessXYZ[18](1,2)= 2.0*z;
    XYZ[19]=x*y*z;
    gradXYZ[19][0]= y*z;
    gradXYZ[19][1]= x*z;
    gradXYZ[19][2]= x*y;
    hessXYZ[19](0,1)=hessXYZ[19](1,0)=z;
    hessXYZ[19](0,2)=hessXYZ[19](2,0)=y;
    hessXYZ[19](2,1)=hessXYZ[19](1,2)=x;
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
  case 2:
    XYZ[4]=x*x;
    gradXYZ[4][0]= 2.0*x;
    hessXYZ[4](0,0)= 2.0;
    XYZ[5]=y*y;
    gradXYZ[5][1]= 2.0*y;
    hessXYZ[5](1,1)= 2.0;
    XYZ[6]=z*z;
    gradXYZ[6][2]= 2.0*z;
    hessXYZ[6](2,2)= 2.0;
    XYZ[7]=x*y;
    gradXYZ[7][0]= y;
    gradXYZ[7][1]= x;
    hessXYZ[7](0,1)=hessXYZ[7](1,0)= 1.0;
    XYZ[8]=x*z;
    gradXYZ[8][0]= z;
    gradXYZ[8][2]= x;
    hessXYZ[8](0,2)=hessXYZ[8](2,0)= 1.0;
    XYZ[9]=y*z;
    gradXYZ[9][1]= z;
    gradXYZ[9][2]= y;
    hessXYZ[9](2,1)=hessXYZ[9](1,2)= 1.0;
  case 1:
    XYZ[1]=x;
    gradXYZ[1][0]= 1.0;
    XYZ[2]=y;
    gradXYZ[2][1]= 1.0;
    XYZ[3]=z;
    gradXYZ[3][2]= 1.0;
  case 0:
    XYZ[0]=1.0;
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

template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::getABC(int n, int& a, int& b, int& c)
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

#endif
