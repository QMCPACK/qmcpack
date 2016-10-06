//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef _CUSPCORR_
#define _CUSPCORR_

#include<cmath>
#include <iostream>
#include "Configuration.h"
//#include "QMCWaveFunctions/BasisSetFactory.h"
//#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/ParticleSet.h"
#include "OhmmsData/OhmmsElementBase.h"
//#include "OhmmsData/Libxml2Doc.h"
//#include "OhmmsData/AttributeSet.h"
//#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"
//#include "QMCWaveFunctions/LocalizedBasisSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{

template<class BS>
class CuspCorr : public QMCTraits
{

  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
  typedef SPOSetBase*        SPOSetBasePtr;

public:

  CuspCorr(RealType r, int nIntPnts, ParticleSet* targetP, ParticleSet* sourceP, bool print=true): Rc_init(r),Rc(r),Rc_max(r),printOrbs(print)
  {
    sourcePtcl=sourceP;
    targetPtcl=targetP;
    coeff.resize(10);
    alpha.resize(10);
    coeff[0] = 3.25819;
    coeff[1] = -15.0126;
    coeff[2] = 33.7308;
    coeff[3] = -42.8705;
    coeff[4] = 31.2276;
    coeff[5] = -12.1316;
    coeff[6] = 1.94692;
    for(int i=0; i<10; i++)
      alpha[i]=0.0;
    C=0.0;
    sg=1.0;
    nElms=nIntPnts;
  }

  ~CuspCorr() {}

/* 
  RealType execute(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, std::string file, RealType cutoff,RealType* data);

  void executeWithRCLoop(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, std::string file, RealType cutoff,RealType* data);

  void fillRadFunWithPhiBar(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, RealType* data);

  void fillRadFunWithPhi(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb);

*/

  RealType loop(RealType phi0, ValueVector_t X)
  {
    valAtZero = phi0;
    sg = phi0>0.0?1.0:-1.0;
    C = (valAtRc*phi0<0.0)?1.5*valAtRc:0.0;
    evalX(X);
    X2alpha(X);
    getELcurr();
    return getchi2();
  }


  void evaluate(SPOSetBasePtr Psi,TinyVector<RealType,3> r, ValueVector_t& val, GradVector_t  &grad, ValueVector_t &lapl)
  {
    targetPtcl->R[0] = sourcePtcl->R[curCenter];
    TinyVector<RealType,3> ddr2=targetPtcl->makeMove(0,r);
    Psi->evaluate(*targetPtcl,0,val,grad,lapl);
  }

  RealType phi( RealType r)
  {
    TinyVector<RealType,3> dr=0;
    dr[0]=r;
    evaluate(Psi1,dr,val1,grad1,lapl1);
    return val1[curOrb];
  }


  RealType idealEL(RealType r)
  {
    RealType r1=r*r,res=coeff[0]*r1;
    r1*=r;
    res+=coeff[1]*r1;
    r1*=r;
    res+=coeff[2]*r1;
    r1*=r;
    res+=coeff[3]*r1;
    r1*=r;
    res+=coeff[4]*r1;
    r1*=r;
    res+=coeff[5]*r1;
    r1*=r;
    res+=coeff[6]*r1;
    return (res+beta0)*Z*Z;
  }

  RealType getELorig()
  {
    RealType dx = Rc*1.2/nElms;
    TinyVector<RealType,3> ddr=0;
    evaluate(Psi2,ddr,val2,grad2,lapl2);  // eta(0)
    evaluate(Psi1,ddr,val1,grad1,lapl1);  // phi(0)
    RealType Zeff = Z*(1.0+ val2[curOrb]/val1[curOrb]);
    for( int i=0; i<nElms; i++)
    {
      ddr[0]=(i+1.0)*dx;
      pos[i]=(i+1.0)*dx;
      evaluate(Psi1,ddr,val1,grad1,lapl1);
      ELorig[i] = -0.5*lapl1[curOrb]/val1[curOrb]-Zeff/pos[i];
    }
    ddr[0]=Rc;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    return -0.5*lapl1[curOrb]/val1[curOrb]-Zeff/Rc;
  }

  void getELideal(RealType el)
  {
    RealType dx = Rc*1.2/nElms;
    beta0=0.0;
    RealType tmp=idealEL(Rc);
    beta0=(el-tmp)/Z/Z;
    for( int i=0; i<nElms; i++)
    {
      pos[i]=(i+1.0)*dx;
      ELideal[i] = idealEL(pos[i]);
    }
  }

  void getELcurr()
  {
    RealType dx = Rc*1.2/nElms;
    RealType Zeff = Z*(1.0+ eta0/phiBar(0.0));
    RealType dp;
    TinyVector<RealType,3> ddr=0;
    ddr[0]=Rc;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    RealType dE = ELorigAtRc - (-0.5*lapl1[curOrb]/val1[curOrb]-Zeff/Rc);
    for( int i=0; i<nElms; i++)
    {
      pos[i]=(i+1.0)*dx;
      if(pos[i] <= Rc)
      {
        dp=dpr(pos[i]);
        ELcurr[i] = -0.5*Rr(pos[i])*(2.0*dp/pos[i]+d2pr(pos[i])+dp*dp)/phiBar(pos[i])-Zeff/pos[i] + dE;
      }
      else
      {
        ddr[0]=pos[i];
        evaluate(Psi1,ddr,val1,grad1,lapl1);
        ELcurr[i] = -0.5*lapl1[curOrb]/val1[curOrb]-Zeff/pos[i] + dE;
      }
    }
  }

  RealType getchi2()
  {
    RealType res=0.0;
    for(int i=0; i<nElms; i++)
      res+=(ELcurr[i]-ELideal[i])*(ELcurr[i]-ELideal[i]);
    return res;
  }

  inline RealType pr(RealType r)
  {
    return alpha[0] + alpha[1]*r + alpha[2]*r*r
           + alpha[3]*r*r*r + alpha[4]*r*r*r*r;
  }

  inline RealType dpr(RealType r)
  {
    return  alpha[1] + 2.0*alpha[2]*r
            + 3.0*alpha[3]*r*r + 4.0*alpha[4]*r*r*r;
  }

  inline RealType d2pr(RealType r)
  {
    return  2.0*alpha[2] + 6.0*alpha[3]*r + 12.0*alpha[4]*r*r;
  }

  inline RealType Rr(RealType r)
  {
    return sg*std::exp(pr(r));
  }

  inline RealType phiBar(RealType r)
  {
    if(r <= Rc)
      return C + Rr(r);
    else
      return phi(r);
  }

  RealType numDeriv(RealType r)
  {
    TinyVector<RealType,3> ddr=0;
    ddr[0]=r+0.0001;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    RealType res = val1[curOrb];
    ddr[0]=r;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    res -= val2[curOrb];
    return res/0.0001;
  }

  RealType num2Deriv(RealType r)
  {
    TinyVector<RealType,3> ddr=0;
    ddr[0]=r+0.0001;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    RealType res = val1[curOrb];
    ddr[0]=r-0.0001;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    res += val1[curOrb];
    ddr[0]=r;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    res -= 2.0*val1[curOrb];
    return res/0.0001/0.0001;
  }


  void evalX(ValueVector_t &X)
  {
    TinyVector<RealType,3> ddr=0;
    ddr[0]=Rc;
    evaluate(Psi1,ddr,val1,grad1,lapl1);
    X[0] = std::log(std::abs(val1[curOrb]-C));
    X[1] = grad1[curOrb][0]/(val1[curOrb]-C);  // since vec{r}=vec{x}
    X[2] = (lapl1[curOrb]-2.0*grad1[curOrb][0]/Rc)/(val1[curOrb]-C);
    X[3] = -Z*(valAtZero+eta0)/(valAtZero-C);
    X[4] = std::log(std::abs(valAtZero-C));
  }

  void X2alpha(const ValueVector_t X)
  {
    RealType RcInv=1.0/Rc, RcInv2=RcInv*RcInv;
    alpha[0] = X[4];
    alpha[1] = X[3];
    alpha[2] = 6.0*X[0]*RcInv2 - 3.0*X[1]*RcInv + X[2]*0.5
               - 3.0*X[3]*RcInv - 6.0*X[4]*RcInv2 - 0.5*X[1]*X[1];
    alpha[3] = -8.0*X[0]*RcInv2*RcInv + 5.0*X[1]*RcInv2 - X[2]*RcInv
               + 3.0*X[3]*RcInv2 + 8.0*X[4]*RcInv2*RcInv + X[1]*X[1]*RcInv;
    alpha[4] = 3.0*X[0]*RcInv2*RcInv2 - 2.0*X[1]*RcInv2*RcInv
               + 0.5*X[2]*RcInv2 - X[3]*RcInv2*RcInv - 3.0*X[4]*RcInv2*RcInv2
               - 0.5*X[1]*X[1]*RcInv2;
  }

#include "QMCWaveFunctions/Experimental/CuspCorr.cpp"

private:
  ///target ParticleSet
  ParticleSet *targetPtcl;
  ///source ParticleSet
  ParticleSet *sourcePtcl;

  SPOSetBasePtr Psi1,Psi2;

  // cutoff
  RealType beta0,DX,eta0, ELorigAtRc;
  RealType Rc_init,Rc,C,sg,Z,valAtZero,valAtRc,Rc_max;
  int nElms,curOrb,curCenter;
  ValueVector_t alpha,coeff;
  bool printOrbs;

  ValueVector_t ELideal;
  ValueVector_t ELorig;
  ValueVector_t ELcurr;

  ValueVector_t val1;
  GradVector_t  grad1;
  ValueVector_t lapl1;

  ValueVector_t val2;
  GradVector_t  grad2;
  ValueVector_t lapl2;

  ValueVector_t pos;

};

}

#endif
