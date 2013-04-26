#ifndef QMCPLUSPLUS_LPQHIBASIS_H
#define QMCPLUSPLUS_LPQHIBASIS_H

#include "LongRange/LRBasis.h"
#include <complex>

namespace qmcplusplus
{

/** @ingroup longrange
 *\brief A derivative of LRBasis class to provide the functionality
 * of the LPQHI basis. Based on code by Ken Esler from PIMC++.
 */

class LPQHIBasis: public LRBasis
{
private:
  int NumKnots; //Number of knots for basis.
  RealType delta, deltainv;
  Matrix<RealType> S; //Coefficients for LPQHI
  Matrix<RealType> S1; //First derivatives
  Matrix<RealType> S2; //Second derivatives
  RealType Mfactor[3];
  vector<RealType> tvec; //Coefficients

  //Helper functions for computing FT of basis functions (used in c(n,k))
  inline complex<RealType> Eplus(int i, RealType k, int n);
  inline complex<RealType> Eminus(int i, RealType k, int n);
  inline RealType Dplus(int i, RealType k, int n);
  inline RealType Dminus(int i, RealType k, int n);

public:

  LPQHIBasis(const LPQHIBasis& b,ParticleLayout_t& ref):
    LRBasis(ref),
    NumKnots(b.NumKnots), delta(b.delta), deltainv(b.deltainv),
    S(b.S),S1(b.S1), S2(b.S2), tvec(b.tvec)
  {
    Mfactor[0]=1.0;
    Mfactor[1]=-1.0;
    Mfactor[2]=1.0;
    BasisSize=b.BasisSize;
    m_rc=b.m_rc;
  }

  inline RealType get_delta() const
  {
    return delta;
  }
  //inline int NumBasisElem() const {return 3*NumKnots;}
  void set_NumKnots(int n); // n >= 2 required
  void set_rc(RealType rc);
  inline RealType h(int n, RealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;
    RealType ra = delta*(i-1);
    RealType rb = delta*i;
    RealType rc = delta*(i+1);
    rc = std::min(m_rc, rc);
    const RealType* restrict Sa(S[alpha]);
    if(r<ra || r>rc)
      return 0.0;
    if (r <= rb)
    {
      RealType x=(rb-r)*deltainv;
      return Mfactor[alpha]*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5])))));
    }
    else
    {
      RealType x=(r-rb)*deltainv;
      return Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5]))));
    }
  }
  inline RealType df(int n, RealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;
    RealType ra = delta*(i-1);
    RealType rb = delta*i;
    RealType rc = delta*(i+1);
    rc = std::min(m_rc, rc);
    const RealType* restrict Sa(S1[alpha]);
    if(r<ra || r>rc)
      return 0.0;
    if (r <= rb)
    {
      RealType x=(rb-r)*deltainv;
      return Mfactor[alpha]*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*Sa[4]))));
    }
    else
    {
      RealType x=(r-rb)*deltainv;
      return Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*Sa[4])));
    }
  }
//    inline TinyVector<RealType,3> getTriplet(int n, RealType r) const {
//      typedef TinyVector<RealType,3> Return_t;
//      int i=n/3;
//      int alpha = n-3*i;
//      RealType ra = delta*(i-1);
//      RealType rb = delta*i;
//      RealType rc = delta*(i+1);
//      rc = std::min(m_rc, rc);
//      if(r<ra || r>rc) return Return_t;
//      const RealType* restrict Sa(S[alpha]);
//      const RealType* restrict S1a(S1[alpha]);
//      const RealType* restrict S2a(S2[alpha]);
//      if (r <= rb) {
//        RealType x=(rb-r)*deltainv;
//        return Return_t(
//            Mfactor[alpha]*(Sa[0] +x*(Sa[1] +x*(Sa[2] +x*(Sa[3] +x*(Sa[4]+x*Sa[5]))))),
//            Mfactor[alpha]*(S1a[0]+x*(S1a[1]+x*(S1a[2]+x*(S1a[3]+x*S1a[4])))),
//            Mfactor[alpha]*(S2a[0]+x*(S2a[1]+x*(S2a[2]+x*S2a[3]))));
//      } else {
//        RealType x=(r-rb)*deltainv;
//        return Return_t(
//            Sa[0] +x*(Sa[1] +x*(Sa[2] +x*(Sa[3] +x*(Sa[4]+x*Sa[5])))),
//            S1a[0]+x*(S1a[1]+x*(S1a[2]+x*(S1a[3]+x*S1a[4]))),
//            S2a[0]+x*(S2a[1]+x*(S2a[2]+x*S2a[3])));
//      }
//    }

  RealType hintr2(int n);
  RealType c(int n, RealType k);
  //Constructor...fill S matrix...call correct base-class constructor
  LPQHIBasis(ParticleLayout_t& ref) : LRBasis(ref), NumKnots(0), delta(0.0)
  {
    S.resize(3,6);
    S(0,0)=1.0;
    S(0,1)=0.0;
    S(0,2)=0.0;
    S(0,3)=-10.0;
    S(0,4)=15.0;
    S(0,5)=-6.0;
    S(1,0)=0.0;
    S(1,1)=1.0;
    S(1,2)=0.0;
    S(1,3)=-6.0;
    S(1,4)=8.0;
    S(1,5)=-3.0;
    S(2,0)=0.0;
    S(2,1)=0.0;
    S(2,2)=0.5;
    S(2,3)=-1.5;
    S(2,4)=1.5;
    S(2,5)=-0.5;
    S1.resize(3,5);
    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        S1(j,i) = static_cast<double>(i+1.0)*S(j,i+1);
      }
    }
    Mfactor[0]=1.0;
    Mfactor[1]=-1.0;
    Mfactor[2]=1.0;
  }
};

}

#endif
