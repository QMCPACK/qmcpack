//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "LongRange/LPQHISRCoulombBasis.h"
#include <cassert>

namespace qmcplusplus
{

using std::cos;
using std::sin;
void
LPQHISRCoulombBasis::set_NumKnots(int n)
{
  assert(n>1);
  NumKnots = n;
  if(m_rc != 0.0)
  {
    delta = m_rc / (NumKnots - 1);
    deltainv = 1.0/delta;
  }
  //set the BasisSize to 3*NumKnots
  BasisSize=3*NumKnots;
}


void
LPQHISRCoulombBasis::set_rc(mRealType rc)
{
  m_rc = rc;
  if(NumKnots != 0)
  {
    delta = m_rc / (NumKnots - 1);
    deltainv = 1.0/delta;
  }
}


//LPQHISRCoulombBasis::mRealType
//LPQHISRCoulombBasis::h(int n, mRealType r) {
//  int i=n/3;
//  int alpha = n-3*i;
//  mRealType ra = delta*(i-1);
//  mRealType rb = delta*i;
//  mRealType rc = delta*(i+1);
//  rc = std::min(m_rc, rc);
//  if ((r > ra) && (r <= rb)) {
//    mRealType sum = 0.0;
//    mRealType prod = 1.0;
//    for (int j=0; j<=5; j++) {
//      sum += (S(alpha,j) * prod);
//      prod *= ((rb - r) * deltainv);
//    }
//    for (int j=0; j<alpha; j++)
//      sum *= -1.0;
//    return (sum);
//  }
//  else if ((r > rb) && (r <= rc)) {
//    mRealType sum = 0.0;
//    mRealType prod = 1.0;
//    for (int j=0; j<=5; j++) {
//      sum += S(alpha,j) * prod;
//      prod *= ((r-rb) * deltainv);
//    }
//    return sum;
//  }
//  return 0.0;
//}


LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::hintr2(int n)
{
  int j=n/3;
  int alpha = n-3*j;
  mRealType deltaa2 = std::pow(delta, alpha+2.0);
  int mysign=1-2*(alpha%2);
  mRealType sum=0.0; 
  
  for (int i=0; i<=5; i++)
  {
     if ( j<NumKnots-1 ) sum+=S(alpha,i)*deltaa2*(1.0/mRealType(i+2)+mRealType(j)/mRealType(i+1));
     if ( j>0 ) sum+=mysign*S(alpha,i)*deltaa2*(-1.0/mRealType(i+2)+mRealType(j)/mRealType(i+1));
  }

  return(sum);
}


LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::c(int m, mRealType k)
{
  int i=m/3;
  int alpha = m-3*i;
  mRealType sum = 0.0;
  if (i == 0)
  {
    for (int n=0; n<=5; n++)
    {
      mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n));
    }
  }
  else
    if (i == (NumKnots-1))
    {
      for (int n=0; n<=5; n++)
      {
        mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dminus(i,k,n)*sign);
      }
    }
    else
    {
      for (int n=0; n<=5; n++)
      {
        mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dplus(i,k,n) + Dminus(i,k,n)*sign);
      }
    }
  return std::pow(delta, alpha)*(sum);
}


inline std::complex<LPQHISRCoulombBasis::mRealType>
LPQHISRCoulombBasis::Eplus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  if (n == 0)
  {
    std::complex<mRealType> e1(cos(k*delta)-1.0, sin(k*delta));
    std::complex<mRealType> e2(cos(k*delta*i),   sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else
  {
    std::complex<mRealType> t1, t2;
    mRealType sign = 1.0;
    t1 = std::complex<mRealType>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    t2 =-(mRealType)n/delta*Eplus(i,k,n-1);;
    return (-(eye/k)*(t1+t2));
  }
}


inline std::complex<LPQHISRCoulombBasis::mRealType>
LPQHISRCoulombBasis::Eminus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  if (n == 0)
  {
    std::complex<mRealType> e1(cos(k*delta)-1.0, -sin(k*delta));
    std::complex<mRealType> e2(cos(k*delta*i),    sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else
  {
    std::complex<mRealType> t1, t2;
    mRealType sign = (n & 1) ? -1.0 : 1.0;
    t1 = sign*
         std::complex<mRealType> (cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    t2 =-(mRealType)n/delta*Eminus(i,k,n-1);
    return (-(eye/k)*(t1+t2));
  }
}


inline LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::Dplus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  std::complex<mRealType> Z1 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Lattice.Volume)*(Z1.imag());
}


inline LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::Dminus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  std::complex<mRealType> Z1 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*(Z1.imag());
}

LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::dc_dk(int m, mRealType k)
{
  int i=m/3;
  int alpha = m-3*i;
  mRealType sum = 0.0;
  if (i == 0)
  {
    for (int n=0; n<=5; n++)
    {
      mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus_dG(i,k,n));

    }
  }
  else
    if (i == (NumKnots-1))
    {
      for (int n=0; n<=5; n++)
      {
        mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dminus_dG(i,k,n)*sign);
       
      }
    }
    else
    {
      for (int n=0; n<=5; n++)
      {
        mRealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dplus_dG(i,k,n) + Dminus_dG(i,k,n)*sign);

      }
    }
  return std::pow(delta, alpha)*(sum);
}

inline std::complex<LPQHISRCoulombBasis::mRealType>
LPQHISRCoulombBasis::Eplus_dG(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  mRealType ri = i*delta;
  mRealType kinv=1/mRealType(k);
  std::complex<mRealType> eigd(cos(k*delta), sin(k*delta));
  std::complex<mRealType> eigr(cos(k*ri),    sin(k*ri));
  
  if (n == 0)
  {

    return Eplus(i,k,n)*(eye*ri - kinv) + delta*kinv*eigr*eigd;
  }
  else
  {
    return -kinv*Eplus(i,k,n) - eye*kinv*(eye*(ri+delta)*eigd*eigr - (n/mRealType(delta))*Eplus_dG(i,k,n-1));
  }
}

inline std::complex<LPQHISRCoulombBasis::mRealType>
LPQHISRCoulombBasis::Eminus_dG(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  mRealType ri = i*delta;
  mRealType kinv=1.0/mRealType(k);
  std::complex<mRealType> eigd(cos(k*delta), -sin(k*delta));
  std::complex<mRealType> eigr(cos(k*ri),    sin(k*ri));

  if (n == 0)
  {
    std::complex<mRealType> eigd(cos(k*delta), -sin(k*delta));
    std::complex<mRealType> eigr(cos(k*ri),    sin(k*ri));
    return Eminus(i,k,n)*(eye*ri - kinv) - delta*kinv*eigr*eigd;
  }
  else
  {
    return -kinv*Eminus(i,k,n) - eye*kinv*(mRealType(pow(-1.0,n))*eye*(ri-delta)*eigd*eigr - (n/mRealType(delta))*Eminus_dG(i,k,n-1));
  }
}


inline LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::Dplus_dG(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  mRealType kinv=1.0/mRealType(k);
  std::complex<mRealType> Z1 = Eplus_dG(i,k,n);

  return 4.0*M_PI/(k*Lattice.Volume)*Z1.imag()- kinv*Dplus(i,k,n);
}


inline LPQHISRCoulombBasis::mRealType
LPQHISRCoulombBasis::Dminus_dG(int i, mRealType k, int n)
{
  mRealType kinv=1.0/mRealType(k);
  std::complex<mRealType> eye(0.0, 1.0);
  std::complex<mRealType> Z1 = Eminus_dG(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*Z1.imag()- kinv*Dminus(i,k,n);
}

}

