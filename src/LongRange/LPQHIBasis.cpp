//Based on code by Ken Esler from PIMC++.

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
LPQHIBasis::set_rc(RealType rc)
{
  m_rc = rc;
  if(NumKnots != 0)
  {
    delta = m_rc / (NumKnots - 1);
    deltainv = 1.0/delta;
  }
}


//LPQHIBasis::RealType
//LPQHIBasis::h(int n, RealType r) {
//  int i=n/3;
//  int alpha = n-3*i;
//  RealType ra = delta*(i-1);
//  RealType rb = delta*i;
//  RealType rc = delta*(i+1);
//  rc = std::min(m_rc, rc);
//  if ((r > ra) && (r <= rb)) {
//    RealType sum = 0.0;
//    RealType prod = 1.0;
//    for (int j=0; j<=5; j++) {
//      sum += (S(alpha,j) * prod);
//      prod *= ((rb - r) * deltainv);
//    }
//    for (int j=0; j<alpha; j++)
//      sum *= -1.0;
//    return (sum);
//  }
//  else if ((r > rb) && (r <= rc)) {
//    RealType sum = 0.0;
//    RealType prod = 1.0;
//    for (int j=0; j<=5; j++) {
//      sum += S(alpha,j) * prod;
//      prod *= ((r-rb) * deltainv);
//    }
//    return sum;
//  }
//  return 0.0;
//}


LPQHIBasis::RealType
LPQHIBasis::hintr2(int n)
{
  int j=n/3;
  int alpha = n-3*j;
  RealType deltap3 = delta*delta*delta;
  RealType min1toalpha=1.0;
  bool alphaeven=true;
  //Constants above correct for alpha==0
  if(alpha == 1)
  {
    min1toalpha = -1.0;
    alphaeven = false;
  }
  RealType sum=0.0;
  if(j==0)
  {
    for(int i=0; i<=5; i++)
      sum += S(alpha,i)/(i+3);
    sum *= deltap3;
  }
  else
    if(j==(NumKnots-1))
    {
      RealType prod1 = 1.0/3.0;
      RealType prod2 = -j;
      RealType prod3 = j*j;
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
        RealType prod1 = 1.0/3.0;
        RealType prod2 = j*j;
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


LPQHIBasis::RealType
LPQHIBasis::c(int m, RealType k)
{
  int i=m/3;
  int alpha = m-3*i;
  RealType sum = 0.0;
  if (i == 0)
  {
    for (int n=0; n<=5; n++)
    {
      RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n));
    }
  }
  else
    if (i == (NumKnots-1))
    {
      for (int n=0; n<=5; n++)
      {
        RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dminus(i,k,n)*sign);
      }
    }
    else
    {
      for (int n=0; n<=5; n++)
      {
        RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dplus(i,k,n) + Dminus(i,k,n)*sign);
      }
    }
  return (sum);
}


inline complex<LPQHIBasis::RealType>
LPQHIBasis::Eplus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  if (n == 0)
  {
    complex<RealType> e1(cos(k*delta)-1.0, sin(k*delta));
    complex<RealType> e2(cos(k*delta*i),   sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else
  {
    complex<RealType> t1, t2;
    RealType sign = 1.0;
    t1 = complex<RealType>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    t2 =-(RealType)n/delta*Eplus(i,k,n-1);;
    return (-(eye/k)*(t1+t2));
  }
}


inline complex<LPQHIBasis::RealType>
LPQHIBasis::Eminus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  if (n == 0)
  {
    complex<RealType> e1(cos(k*delta)-1.0, -sin(k*delta));
    complex<RealType> e2(cos(k*delta*i),    sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else
  {
    complex<RealType> t1, t2;
    RealType sign = (n & 1) ? -1.0 : 1.0;
    t1 = sign*
         complex<RealType> (cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    t2 =-(RealType)n/delta*Eminus(i,k,n-1);
    return (-(eye/k)*(t1+t2));
  }
}


inline LPQHIBasis::RealType
LPQHIBasis::Dplus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  complex<RealType> Z1 = Eplus(i,k,n+1);
  complex<RealType> Z2 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Lattice.Volume)*(delta* Z1.imag() + i*delta*Z2.imag());
}


inline LPQHIBasis::RealType
LPQHIBasis::Dminus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  complex<RealType> Z1 = Eminus(i,k,n+1);
  complex<RealType> Z2 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*(delta* Z1.imag() + i*delta*Z2.imag());
}
}

