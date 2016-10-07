//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "LongRange/LPQHIBasis.h"
#include <cassert>

namespace qmcplusplus
{

using std::cos;
using std::sin;
void
LPQHIBasis::set_NumKnots(int n)
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
LPQHIBasis::set_rc(mRealType rc)
{
  m_rc = rc;
  if(NumKnots != 0)
  {
    delta = m_rc / (NumKnots - 1);
    deltainv = 1.0/delta;
  }
}


//LPQHIBasis::RealType
//LPQHIBasis::h(int n, mRealType r) {
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


LPQHIBasis::mRealType
LPQHIBasis::hintr2(int n)
{
  int j=n/3;
  int alpha = n-3*j;
  mRealType deltap3 = delta*delta*delta;
  mRealType min1toalpha=1.0;
  bool alphaeven=true;
  //Constants above correct for alpha==0
  if(alpha == 1)
  {
    min1toalpha = -1.0;
    alphaeven = false;
  }
  mRealType sum=0.0;
  if(j==0)
  {
    for(int i=0; i<=5; i++)
      sum += S(alpha,i)/(i+3);
    sum *= deltap3;
  }
  else
    if(j==(NumKnots-1))
    {
      mRealType prod1 = 1.0/3.0;
      mRealType prod2 = -j;
      mRealType prod3 = j*j;
      for(int i=0; i<=5; i++)
      {
        sum += S(alpha,i)*(prod1+prod2+prod3);
        prod1 *= (i+3.0)/(i+4.0);
        prod2 *= (i+2.0)/(i+3.0);
        prod3 *= (i+1.0)/(i+2.0);
      }
      sum *= deltap3*min1toalpha;
    }
    else
      // expression for 0<j<M
    {
      if(alphaeven)
      {
        mRealType prod1 = 1.0/3.0;
        mRealType prod2 = j*j;
        for(int i=0; i<=5; i++)
        {
          sum += S(alpha,i)*(prod1+prod2);
          prod1 *= (i+3.0)/(i+4.0); //Prepare for next cycle.
          prod2 *= (i+1.0)/(i+2.0); //Prepare for next cycle.
        }
        sum *= 2.*deltap3;
      }
      else
      {
        for(int i=0; i<=5; i++)
          sum += S(alpha,i)/(i+2.0);
        sum *= deltap3*4.*j;
      }
    }
  return(sum);
}


LPQHIBasis::mRealType
LPQHIBasis::c(int m, mRealType k)
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
  return (sum);
}


inline std::complex<LPQHIBasis::mRealType>
LPQHIBasis::Eplus(int i, mRealType k, int n)
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


inline std::complex<LPQHIBasis::mRealType>
LPQHIBasis::Eminus(int i, mRealType k, int n)
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


inline LPQHIBasis::mRealType
LPQHIBasis::Dplus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  std::complex<mRealType> Z1 = Eplus(i,k,n+1);
  std::complex<mRealType> Z2 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Lattice.Volume)*(delta* Z1.imag() + i*delta*Z2.imag());
}


inline LPQHIBasis::mRealType
LPQHIBasis::Dminus(int i, mRealType k, int n)
{
  std::complex<mRealType> eye(0.0, 1.0);
  std::complex<mRealType> Z1 = Eminus(i,k,n+1);
  std::complex<mRealType> Z2 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*(delta* Z1.imag() + i*delta*Z2.imag());
}
}

