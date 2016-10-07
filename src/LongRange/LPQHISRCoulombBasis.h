//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_LPQHISRCOULOMBBASIS_H
#define QMCPLUSPLUS_LPQHISRCOULOMBBASIS_H

#include "LongRange/LRBasis.h"
#include <complex>

namespace qmcplusplus
{

/** @ingroup longrange
 *\brief A derivative of LRBasis class to provide the functionality
 * of the LPQHI basis.  A 1/r factor is incorporated into the basis functions
 * to faciliate real space representation of the short-ranged coulomb potential, 
 * following Natoli & Ceperley's 1995 paper on Optimized Breakup.  http://dx.doi.org/10.1006/jcph.1995.1054.
 */

class LPQHISRCoulombBasis: public LRBasis
{
private:
  int NumKnots; //Number of knots for basis.
  mRealType delta, deltainv;
  Matrix<mRealType> S; //Coefficients for LPQHI
  Matrix<mRealType> S1; //First derivatives
  Matrix<mRealType> S2; //Second derivatives
  mRealType Mfactor[3];
  std::vector<mRealType> tvec; //Coefficients

  //Helper functions for computing FT of basis functions (used in c(n,k))
  inline std::complex<mRealType> Eplus(int i, mRealType k, int n);
  inline std::complex<mRealType> Eminus(int i, mRealType k, int n);
  inline std::complex<mRealType> Eplus_dG(int i, mRealType k, int n);
  inline std::complex<mRealType> Eminus_dG(int i, mRealType k, int n);
  inline mRealType Dplus(int i, mRealType k, int n);
  inline mRealType Dminus(int i, mRealType k, int n);
  inline mRealType Dplus_dG(int i, mRealType k, int n);
  inline mRealType Dminus_dG(int i, mRealType k, int n);
  

public:

  LPQHISRCoulombBasis(const LPQHISRCoulombBasis& b,ParticleLayout_t& ref):
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

  inline mRealType get_delta() const
  {
    return delta;
  }
  //inline int NumBasisElem() const {return 3*NumKnots;}
  void set_NumKnots(int n); // n >= 2 required
  void set_rc(mRealType rc);
  inline mRealType h(int n, mRealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;
    mRealType ra = delta*(i-1);
    mRealType rb = delta*i;
    mRealType rc = delta*(i+1);
    mRealType rinv = 1.0/r;
    rc = std::min(m_rc, rc);
    const mRealType* restrict Sa(S[alpha]);
    if(r<ra || r>rc)
      return 0.0;
    if (r <= rb)
    {
      mRealType x=(rb-r)*deltainv;
      return std::pow(delta,alpha)*Mfactor[alpha]*rinv*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5])))));
    }
    else
    {
      mRealType x=(r-rb)*deltainv;
      return std::pow(delta,alpha)*rinv*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5])))));
    }
  }
  inline mRealType rh(int n, mRealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;
    mRealType ra = delta*(i-1);
    mRealType rb = delta*i;
    mRealType rc = delta*(i+1);
    mRealType rinv = 1.0/r;
    rc = std::min(m_rc, rc);
    const mRealType* restrict Sa(S[alpha]);
    if(r<ra || r>rc)
      return 0.0;
    if (r <= rb)
    {
      mRealType x=(rb-r)*deltainv;
      return std::pow(delta,alpha)*Mfactor[alpha]*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5])))));
    }
    else
    {
      mRealType x=(r-rb)*deltainv;
      return std::pow(delta,alpha)*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*(Sa[4]+x*Sa[5])))));
    }
  }
  
  inline mRealType dh_dr(int n, mRealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;
    mRealType ra = delta*(i-1);
    mRealType rb = delta*i;
    mRealType rc = delta*(i+1);
    rc = std::min(m_rc, rc);

    mRealType polyderiv=0.0;
    mRealType rinv=1.0/r;
    mRealType hval=h(n,r);
  
    const mRealType* restrict Sa(S1[alpha]);
    if(r<ra || r>rc)
      return 0.0;
    if (r <= rb)
    {
      mRealType x=(rb-r)*deltainv;
      polyderiv=-std::pow(delta,alpha-1)*Mfactor[alpha]*( Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*Sa[4]))));
    }
    else
    {
      mRealType x=(r-rb)*deltainv;
      polyderiv = std::pow(delta,alpha-1)*(Sa[0]+x*(Sa[1]+x*(Sa[2]+x*(Sa[3]+x*Sa[4]))));
    }

    return rinv*(polyderiv-hval);
  }
  
  inline mRealType dh_ddelta(int n, mRealType r) const
  {
    int i=n/3;
    int alpha = n-3*i;

	return h(n,r)*alpha/double(delta) - dh_dr(n,r)*r/double(delta);
  }
  
  
//    inline TinyVector<mRealType,3> getTriplet(int n, mRealType r) const {
//      typedef TinyVector<mRealType,3> Return_t;
//      int i=n/3;
//      int alpha = n-3*i;
//      mRealType ra = delta*(i-1);
//      mRealType rb = delta*i;
//      mRealType rc = delta*(i+1);
//      rc = std::min(m_rc, rc);
//      if(r<ra || r>rc) return Return_t;
//      const mRealType* restrict Sa(S[alpha]);
//      const mRealType* restrict S1a(S1[alpha]);
//      const mRealType* restrict S2a(S2[alpha]);
//      if (r <= rb) {
//        mRealType x=(rb-r)*deltainv;
//        return Return_t(
//            Mfactor[alpha]*(Sa[0] +x*(Sa[1] +x*(Sa[2] +x*(Sa[3] +x*(Sa[4]+x*Sa[5]))))),
//            Mfactor[alpha]*(S1a[0]+x*(S1a[1]+x*(S1a[2]+x*(S1a[3]+x*S1a[4])))),
//            Mfactor[alpha]*(S2a[0]+x*(S2a[1]+x*(S2a[2]+x*S2a[3]))));
//      } else {
//        mRealType x=(r-rb)*deltainv;
//        return Return_t(
//            Sa[0] +x*(Sa[1] +x*(Sa[2] +x*(Sa[3] +x*(Sa[4]+x*Sa[5])))),
//            S1a[0]+x*(S1a[1]+x*(S1a[2]+x*(S1a[3]+x*S1a[4]))),
//            S2a[0]+x*(S2a[1]+x*(S2a[2]+x*S2a[3])));
//      }
//    }

  mRealType hintr2(int n);
  mRealType c(int n, mRealType k);
  mRealType dc_dk(int n, mRealType k);
  //Constructor...fill S matrix...call correct base-class constructor
  LPQHISRCoulombBasis(ParticleLayout_t& ref) : LRBasis(ref), NumKnots(0), delta(0.0)
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
