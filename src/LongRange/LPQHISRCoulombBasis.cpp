//Based on code by Ken Esler from PIMC++.

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
LPQHISRCoulombBasis::set_rc(RealType rc)
{
  m_rc = rc;
  if(NumKnots != 0)
  {
    delta = m_rc / (NumKnots - 1);
    deltainv = 1.0/delta;
  }
}


//LPQHISRCoulombBasis::RealType
//LPQHISRCoulombBasis::h(int n, RealType r) {
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


LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::hintr2(int n)
{
  int j=n/3;
  int alpha = n-3*j;
  RealType deltaa2 = std::pow(delta, alpha+2.0);
  int mysign=1-2*(alpha%2);
  RealType sum=0.0; 
  
  for (int i=0; i<=5; i++)
  {
     if ( j<NumKnots-1 ) sum+=S(alpha,i)*deltaa2*(1.0/RealType(i+2)+RealType(j)/RealType(i+1));
     if ( j>0 ) sum+=mysign*S(alpha,i)*deltaa2*(-1.0/RealType(i+2)+RealType(j)/RealType(i+1));
  }

  return(sum);
}


LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::c(int m, RealType k)
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
  return std::pow(delta, alpha)*(sum);
}


inline complex<LPQHISRCoulombBasis::RealType>
LPQHISRCoulombBasis::Eplus(int i, RealType k, int n)
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


inline complex<LPQHISRCoulombBasis::RealType>
LPQHISRCoulombBasis::Eminus(int i, RealType k, int n)
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


inline LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::Dplus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  complex<RealType> Z1 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Lattice.Volume)*(Z1.imag());
}


inline LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::Dminus(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  complex<RealType> Z1 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*(Z1.imag());
}

LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::dc_dk(int m, RealType k)
{
  int i=m/3;
  int alpha = m-3*i;
  RealType sum = 0.0;
  if (i == 0)
  {
    for (int n=0; n<=5; n++)
    {
      RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus_dG(i,k,n));

    }
  }
  else
    if (i == (NumKnots-1))
    {
      for (int n=0; n<=5; n++)
      {
        RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dminus_dG(i,k,n)*sign);
       
      }
    }
    else
    {
      for (int n=0; n<=5; n++)
      {
        RealType sign = ((alpha+n)&1) ? -1.0 : 1.0;
        sum += S(alpha, n) * (Dplus_dG(i,k,n) + Dminus_dG(i,k,n)*sign);

      }
    }
  return std::pow(delta, alpha)*(sum);
}

inline complex<LPQHISRCoulombBasis::RealType>
LPQHISRCoulombBasis::Eplus_dG(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  RealType ri = i*delta;
  RealType kinv=1/RealType(k);
  complex<RealType> eigd(cos(k*delta), sin(k*delta));
  complex<RealType> eigr(cos(k*ri),    sin(k*ri));
  
  if (n == 0)
  {

    return Eplus(i,k,n)*(eye*ri - kinv) + delta*kinv*eigr*eigd;
  }
  else
  {
    return -kinv*Eplus(i,k,n) - eye*kinv*(eye*(ri+delta)*eigd*eigr - (n/RealType(delta))*Eplus_dG(i,k,n-1));
  }
}

inline complex<LPQHISRCoulombBasis::RealType>
LPQHISRCoulombBasis::Eminus_dG(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  RealType ri = i*delta;
  RealType kinv=1.0/RealType(k);
  complex<RealType> eigd(cos(k*delta), -sin(k*delta));
  complex<RealType> eigr(cos(k*ri),    sin(k*ri));

  if (n == 0)
  {
    complex<RealType> eigd(cos(k*delta), -sin(k*delta));
    complex<RealType> eigr(cos(k*ri),    sin(k*ri));
    return Eminus(i,k,n)*(eye*ri - kinv) - delta*kinv*eigr*eigd;
  }
  else
  {
    return -kinv*Eminus(i,k,n) - eye*kinv*(RealType(pow(-1.0,n))*eye*(ri-delta)*eigd*eigr - (n/RealType(delta))*Eminus_dG(i,k,n-1));
  }
}


inline LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::Dplus_dG(int i, RealType k, int n)
{
  complex<RealType> eye(0.0, 1.0);
  RealType kinv=1.0/RealType(k);
  complex<RealType> Z1 = Eplus_dG(i,k,n);

  return 4.0*M_PI/(k*Lattice.Volume)*Z1.imag()- kinv*Dplus(i,k,n);
}


inline LPQHISRCoulombBasis::RealType
LPQHISRCoulombBasis::Dminus_dG(int i, RealType k, int n)
{
  RealType kinv=1.0/RealType(k);
  complex<RealType> eye(0.0, 1.0);
  complex<RealType> Z1 = Eminus_dG(i,k,n);
  return -4.0*M_PI/(k*Lattice.Volume)*Z1.imag()- kinv*Dminus(i,k,n);
}

}

