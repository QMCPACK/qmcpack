//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
//#include "Numerics/OhmmsBlas.h"
//#include "Numerics/MatrixOperators.h"
#include "simd/simd.hpp"

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
RNDiracDeterminantBase::RNDiracDeterminantBase(SPOSetBasePtr const &spos, int first):
  DiracDeterminantBase(spos, first), logepsilon(0.0)
{
  OrbitalName="RNDiracDeterminantBase";
}

RNDiracDeterminantBase::RNDiracDeterminantBase(const RNDiracDeterminantBase& s):
  DiracDeterminantBase(s), logepsilon(0.0)
{
  this->resize(s.NumPtcls,s.NumOrbitals);
  setLogEpsilon(s.logepsilon);
//     app_log()<<"setting logepsilon "<<s.logepsilon<<" "<<logepsilon<<endl;
}

///default destructor
RNDiracDeterminantBase::~RNDiracDeterminantBase() {}

RNDiracDeterminantBase& RNDiracDeterminantBase::operator=(const RNDiracDeterminantBase& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  setLogEpsilon(s.logepsilon);
  return *this;
}

void RNDiracDeterminantBase::resize(int nel, int morb)
{
  int norb=morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  dpsiM.resize(nel,norb);
  d2psiM.resize(nel,norb);
  psiM_temp.resize(nel,norb);
  dpsiM_temp.resize(nel,norb);
  d2psiM_temp.resize(nel,norb);
  psiMinv.resize(nel,norb);
  psiV.resize(norb);
  myG_alternate.resize(nel);
  myL_alternate.resize(nel);
  WorkSpace.resize(nel);
  Pivot.resize(nel);
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;
  // For forces
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
}

RNDiracDeterminantBase::RealType
RNDiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  if (NP == 0) //first time, allocate once
  {
    //int norb = cols();
    dpsiV.resize(NumOrbitals);
    d2psiV.resize(NumOrbitals);
    workV1.resize(NumOrbitals);
    workV2.resize(NumOrbitals);
    NP=P.getTotalNum();
    myG.resize(NP);
    myL.resize(NP);
    myG_temp.resize(NP);
    myL_temp.resize(NP);
    myG_alternate.resize(NP);
    myL_alternate.resize(NP);
    FirstAddressOfG = &myG[0][0];
    LastAddressOfG = FirstAddressOfG + NP*DIM;
    FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  }
  myG=0.0;
  myL=0.0;
  myG_alternate=0.0;
  myL_alternate=0.0;
  //ValueType x=evaluate(P,myG,myL);
  evaluateLog(P,myG,myL);
  P.G += myG;
  P.L += myL;
  //add the data: determinant, inverse, gradient and laplacians
  buf.add(psiM.first_address(),psiM.last_address());
  buf.add(FirstAddressOfdV,LastAddressOfdV);
  buf.add(d2psiM.first_address(),d2psiM.last_address());
  buf.add(myL.first_address(), myL.last_address());
  buf.add(FirstAddressOfG,LastAddressOfG);
  buf.add(LogValue);
  buf.add(alternateLogValue);
  buf.add(alternatePhaseValue);
  return LogValue;
}

RNDiracDeterminantBase::RealType RNDiracDeterminantBase::updateBuffer(ParticleSet& P,
    PooledData<RealType>& buf, bool fromscratch)
{
  myG=0.0;
  myL=0.0;
  myG_alternate=0.0;
  myL_alternate=0.0;
  if (fromscratch)
  {
    evaluateLog(P,myG,myL);
    UpdateTimer.start();
  }
  else
  {
    if (UpdateMode == ORB_PBYP_RATIO)
      Phi->evaluate(P, FirstIndex, LastIndex, psiM_temp,dpsiM, d2psiM);
    RealType cp = std::exp(logepsilon -2.0*alternateLogValue);
    RealType bp = 1.0/(1+cp);
    UpdateTimer.start();
    if (NumPtcls==1)
    {
      ValueType y=psiM(0,0);
      GradType rv = y*dpsiM(0,0);
      ValueType rv2=dot(rv,rv);
      myG_alternate(FirstIndex) += rv;
      myL_alternate(FirstIndex) += y*d2psiM(0,0) - rv2;
      myG(FirstIndex) += bp*rv;
      myL(FirstIndex) += bp*(y*d2psiM(0,0) + (1-2*bp)*rv2);
    }
    else
    {
      const ValueType* restrict yptr=psiM.data();
      const ValueType* restrict d2yptr=d2psiM.data();
      const GradType* restrict dyptr=dpsiM.data();
      for (int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
      {
        GradType rv;
        ValueType lap=0.0;
        for (int j=0; j<NumOrbitals; j++,yptr++)
        {
          rv += *yptr * *dyptr++;
          lap += *yptr * *d2yptr++;
        }
        ValueType rv2=dot(rv,rv);
        myG_alternate(iat) += rv;
        myL_alternate(iat) += lap - rv2;
        myG(iat) += bp*rv;
        myL(iat) += bp*(lap + (1-2*bp)*rv2);
      }
    }
  }
  P.G += myG;
  P.L += myL;
  //copy psiM to psiM_temp
  psiM_temp=psiM;
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(LogValue);
  buf.put(alternateLogValue);
  buf.put(alternatePhaseValue);
  UpdateTimer.stop();
  return LogValue;
}

void RNDiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.get(psiM.first_address(),psiM.last_address());
  buf.get(FirstAddressOfdV,LastAddressOfdV);
  buf.get(d2psiM.first_address(),d2psiM.last_address());
  buf.get(myL.first_address(), myL.last_address());
  buf.get(FirstAddressOfG,LastAddressOfG);
  buf.get(LogValue);
  buf.get(alternateLogValue);
  buf.get(alternatePhaseValue);
  //re-evaluate it for testing
  //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
  //need extra copy for gradient/laplacian calculations without updating it
  psiM_temp = psiM;
  dpsiM_temp = dpsiM;
  d2psiM_temp = d2psiM;
}


/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
RNDiracDeterminantBase::ValueType RNDiracDeterminantBase::ratio(ParticleSet& P, int iat)
{
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  Phi->evaluate(P, iat, psiV);
  RatioTimer.start();
  alternateCurRatio = DetRatioByRow(psiM, psiV,WorkingIndex);
  if (abs(alternateCurRatio)< numeric_limits<RealType>::epsilon())
  {
    app_log()<<"stepped on node: ratioGrad"<<endl;
    RatioTimer.stop();
    return 0.0;
  }
  RealType R = std::abs(alternateCurRatio);
  RealType logR = std::log(R);
  RealType bp = 1+std::exp(logepsilon-2.0*alternateLogValue-2.0*logR);
  curRatio=R*std::sqrt(bp/(1+std::exp(logepsilon-2.0*alternateLogValue)));
  RatioTimer.stop();
  return curRatio;
}

RNDiracDeterminantBase::ValueType RNDiracDeterminantBase::alternateRatio(ParticleSet& P)
{
  //returns psi_T/psi_G
  for (int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
  {
    P.G(iat) += myG_alternate(iat) - myG(iat);
    P.L(iat) += myL_alternate(iat) - myL(iat);
  }
  RealType sgn = std::cos(alternatePhaseValue);
  return sgn*std::exp(alternateLogValue-LogValue);
}

void RNDiracDeterminantBase::alternateGrad(ParticleSet::ParticleGradient_t& G)
{
  //returns psi_T/psi_G
  for (int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    G(iat) += myG_alternate(iat);
}


RNDiracDeterminantBase::GradType
RNDiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
{
  WorkingIndex = iat-FirstIndex;
  RatioTimer.start();
  RNDiracDeterminantBase::GradType g = simd::dot(psiM[WorkingIndex],dpsiM[WorkingIndex],NumOrbitals);
  g *= 1.0/(1.0+std::exp(logepsilon-2.0*alternateLogValue));
  RatioTimer.stop();
  return g;
}

RNDiracDeterminantBase::GradType
RNDiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
                                       int iat)
{
  APP_ABORT("What?! released node and forces?");
  return GradType();
}

RNDiracDeterminantBase::GradType
RNDiracDeterminantBase::evalGradSourcep
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  APP_ABORT("What?! released node and forces?");
  return GradType();
}


RNDiracDeterminantBase::GradType
RNDiracDeterminantBase::evalGradSource
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  APP_ABORT("What?! released node and forces?");
  return GradType();
}

void RNDiracDeterminantBase::restore(int iat)
{
  if(UpdateMode == ORB_PBYP_ALL)
  {
    psiM_temp = psiM;
    std::copy(dpsiM[WorkingIndex],dpsiM[WorkingIndex+1],dpsiM_temp[WorkingIndex]);
    std::copy(d2psiM[WorkingIndex],d2psiM[WorkingIndex+1],d2psiM_temp[WorkingIndex]);
  }
  curRatio=1.0;
  alternateCurRatio=1.0;
}


RNDiracDeterminantBase::ValueType
RNDiracDeterminantBase::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  RatioTimer.start();
  WorkingIndex = iat-FirstIndex;
  UpdateMode=ORB_PBYP_PARTIAL;
  alternateCurRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
  if (abs(alternateCurRatio)< numeric_limits<RealType>::epsilon())
  {
    app_log()<<"stepped on node: ratioGrad"<<endl;
    RatioTimer.stop();
    return 0.0;
  }
  RealType R = std::abs(alternateCurRatio);
  RealType logR = std::log(R);
  RealType bp = 1+std::exp(logepsilon-2.0*alternateLogValue-2.0*logR) ;
  curRatio=R*std::sqrt(bp/(1.0+std::exp(logepsilon-2.0*alternateLogValue)));
//     bp=std::sqrt(1.0/bp);
  bp = 1.0/bp;
  GradType rv=simd::dot(psiM[WorkingIndex],dpsiV.data(),NumOrbitals);
  grad_iat += (1.0/alternateCurRatio) *bp* rv;
  RatioTimer.stop();
  return curRatio;
  ////////////////////////////////////////
  ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
  //{
  //  int kat=FirstIndex;
  //  for(int j=0; j<NumOrbitals; j++) {
  //    dpsiM_temp(WorkingIndex,j)=dpsiV[j];
  //    d2psiM_temp(WorkingIndex,j)=d2psiV[j];
  //  }
  //  const ValueType* restrict yptr=psiM_temp.data();
  //  const ValueType* restrict d2yptr=d2psiM_temp.data();
  //  const GradType* restrict dyptr=dpsiM_temp.data();
  //  for(int i=0; i<NumPtcls; i++,kat++) {
  //    //This mimics gemm with loop optimization
  //    GradType rv;
  //    ValueType lap=0.0;
  //    for(int j=0; j<NumOrbitals; j++,yptr++) {
  //      rv += *yptr * *dyptr++;
  //      lap += *yptr * *d2yptr++;
  //    }
  //    lap -= dot(rv,rv);
  //    myG_temp[kat]=rv;
  //    myL_temp[kat]=lap;
  //  }
  //}
  ////////////////////////////////////////
}

/** return the ratio
 * @param P current configuration
 * @param iat particle whose position is moved
 * @param dG differential Gradients
 * @param dL differential Laplacians
 *
 * Data member *_temp contain the data assuming that the move is accepted
 * and are used to evaluate differential Gradients and Laplacians.
 */
RNDiracDeterminantBase::ValueType RNDiracDeterminantBase::ratio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  UpdateMode=ORB_PBYP_ALL;
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  RatioTimer.start();
  WorkingIndex = iat-FirstIndex;
  //psiM_temp = psiM;
  alternateCurRatio= DetRatioByRow(psiM_temp, psiV, WorkingIndex);
  RatioTimer.stop();
  if (abs(alternateCurRatio)<numeric_limits<RealType>::epsilon())
  {
    app_log()<<"stepped on node"<<endl;
    UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse
    return 0.0;
  }
  RealType R = std::abs(alternateCurRatio);
  RealType logR = std::log(R);
  RealType bp = 1+std::exp(logepsilon-2.0*alternateLogValue-2.0*logR);
  curRatio = R*std::sqrt(bp/(1+std::exp(logepsilon-2.0*alternateLogValue)));
  bp = 1.0/bp;
  UpdateTimer.start();
  //update psiM_temp with the row substituted
  InverseUpdateByRow(psiM_temp,psiV,workV1,workV2,WorkingIndex,alternateCurRatio);
  //update dpsiM_temp and d2psiM_temp
  std::copy(dpsiV.begin(),dpsiV.end(),dpsiM_temp[WorkingIndex]);
  std::copy(d2psiV.begin(),d2psiV.end(),d2psiM_temp[WorkingIndex]);
  UpdateTimer.stop();
  RatioTimer.start();
  for (int i=0,kat=FirstIndex; i<NumPtcls; i++,kat++)
  {
    //using inline dot functions
    GradType rv=simd::dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
    ValueType lap=simd::dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);
    ValueType rv2 = dot(rv,rv);
    dG[kat] += bp*rv - myG[kat];
    myG_temp[kat]= bp*rv;
    dL[kat] += bp*(lap +(1.0-2.0*bp)*rv2) -myL[kat];
    myL_temp[kat]= bp*(lap +(1.0-2.0*bp)*rv2) ;
  }
  RatioTimer.stop();
  return curRatio;
}

/** move was accepted, update the real container
*/
void RNDiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
{
  alternatePhaseValue += evaluatePhase(alternateCurRatio);
  alternateLogValue +=std::log(std::abs(alternateCurRatio));
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  switch (UpdateMode)
  {
  case ORB_PBYP_RATIO:
    InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,alternateCurRatio);
    break;
  case ORB_PBYP_PARTIAL:
    //psiM = psiM_temp;
    InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,alternateCurRatio);
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    //////////////////////////////////////
    ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
    //myG = myG_temp;
    //myL = myL_temp;
    ///////////////////////
    break;
  default:
    myG = myG_temp;
    myL = myL_temp;
    psiM = psiM_temp;
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    break;
  }
  UpdateTimer.stop();
  alternateCurRatio=1.0;
  curRatio=1.0;
}


void RNDiracDeterminantBase::update(ParticleSet& P,
                                    ParticleSet::ParticleGradient_t& dG,
                                    ParticleSet::ParticleLaplacian_t& dL,
                                    int iat)
{
  UpdateTimer.start();
  InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,alternateCurRatio);
  for (int j=0; j<NumOrbitals; j++)
  {
    dpsiM(WorkingIndex,j)=dpsiV[j];
    d2psiM(WorkingIndex,j)=d2psiV[j];
  }
  UpdateTimer.stop();
  RatioTimer.start();
  int kat=FirstIndex;
  RealType R = std::abs(alternateCurRatio);
  RealType logR = std::log(R);
  RealType bp = 1+std::exp(logepsilon-2.0*alternateLogValue-2.0*logR);
  curRatio = R*std::sqrt(bp/(1+std::exp(logepsilon-2.0*alternateLogValue)));
  bp = 1.0/bp;
  for (int i=0; i<NumPtcls; i++,kat++)
  {
    GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
    ValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
    ValueType rv2 = dot(rv,rv);
    dG[kat] += bp*rv - myG[kat];
    myG_alternate[kat] =rv;
    dL[kat] += bp*(lap +(1-2.0*bp)*rv2) -myL[kat];
    myL_alternate[kat] =lap-rv2;
  }
  RatioTimer.stop();
  alternatePhaseValue += evaluatePhase(alternateCurRatio);
  alternateLogValue +=std::log(std::abs(alternateCurRatio));
  LogValue +=std::log(std::abs(curRatio));
  curRatio=1.0;
  alternateCurRatio=1.0;
}

RNDiracDeterminantBase::RealType
RNDiracDeterminantBase::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(LogValue);
  buf.put(alternateLogValue);
  buf.put(alternatePhaseValue);
  return LogValue;
}



RNDiracDeterminantBase::RealType
RNDiracDeterminantBase::evaluateLog(ParticleSet& P,
                                    ParticleSet::ParticleGradient_t& G,
                                    ParticleSet::ParticleLaplacian_t& L)
{
  //      cerr<<"I'm calling evaluate log"<<endl;
  Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  myG_alternate=0.0;
  myL_alternate=0.0;
//     myG=0.0;
//     myL=0.0;
  if (NumPtcls==1)
  {
    //CurrentDet=psiM(0,0);
    ValueType det=psiM(0,0);
    alternateLogValue = evaluateLogAndPhase(det,alternatePhaseValue);
//         alternatePhaseValue=0.0;
    LogValue = alternateLogValue + 0.5*std::log(1.0+std::exp(logepsilon-2.0*alternateLogValue));
    RealType cp = std::exp(logepsilon -2.0*alternateLogValue);
    RealType bp = 1.0/(1+cp);
    ValueType y=1.0/det;
    psiM(0,0)=y;
    GradType rv = y*dpsiM(0,0);
    ValueType rv2=dot(rv,rv);
    myG_alternate(FirstIndex) += rv;
    myL_alternate(FirstIndex) += y*d2psiM(0,0) - rv2;
    G(FirstIndex) += bp*rv;
    L(FirstIndex) += bp*(y*d2psiM(0,0) + (1-2*bp)*rv2);
//         myG(FirstIndex) += bp*rv;
//         myL(FirstIndex) += bp*(y*d2psiM(0,0) + (1-2*bp)*dot(rv,rv));
  }
  else
  {
    InverseTimer.start();
    alternateLogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),alternatePhaseValue);
//         alternatePhaseValue=0.0;
    LogValue = alternateLogValue + 0.5*std::log(1.0+std::exp(logepsilon-2.0*alternateLogValue));
    RealType cp = std::exp(logepsilon -2.0*alternateLogValue);
    RealType bp = 1.0/(1+cp);
    InverseTimer.stop();
    RatioTimer.start();
    for (int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      ValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
      G(iat) += bp*rv;
      myG_alternate(iat) += rv;
      ValueType rv2=dot(rv,rv);
      myL_alternate(iat) += lap - rv2;
      L(iat) += bp*(lap + (1-2*bp)*rv2);
    }
    RatioTimer.stop();
  }
  psiM_temp = psiM;
  return LogValue;
}

RNDiracDeterminantBase::DiracDeterminantBase* RNDiracDeterminantBase::makeCopy(SPOSetBase* spo) const
{
  RNDiracDeterminantBase* dclone= new RNDiracDeterminantBase(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  dclone->setLogEpsilon(logepsilon);
  return dclone;
}



}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4433 $   $Date: 2009-12-02 15:58:06 -0600 (Wed, 02 Dec 2009) $
 * $Id: RNDiracDeterminantBase.cpp 4433 2009-12-02 21:58:06Z jmcminis $
 ***************************************************************************/
