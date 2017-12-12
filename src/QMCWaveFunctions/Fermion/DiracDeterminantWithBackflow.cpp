//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsPETE/Tensor.h"
#include <simd/simd.hpp>

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminantWithBackflow::DiracDeterminantWithBackflow(ParticleSet &ptcl, SPOSetBasePtr const &spos, BackflowTransformation * BF, int first): DiracDeterminantBase(spos,first)
{
  Optimizable=true;
  OrbitalName="DiracDeterminantWithBackflow";
  registerTimers();
  BFTrans=BF;
  NumParticles = ptcl.getTotalNum();
  NP=0;
}

///default destructor
DiracDeterminantWithBackflow::~DiracDeterminantWithBackflow() {}

DiracDeterminantWithBackflow& DiracDeterminantWithBackflow::operator=(const DiracDeterminantWithBackflow& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}


///reset the size: with the number of particles and number of orbtials
void DiracDeterminantWithBackflow::resize(int nel, int morb)
{
  int norb=morb;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  dpsiM.resize(nel,norb);
  psiMinv.resize(nel,norb);
  WorkSpace.resize(nel);
  Pivot.resize(nel);
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;
  grad_grad_psiM.resize(nel,norb);
  grad_grad_grad_psiM.resize(nel,norb);
  Gtemp.resize(NumParticles); // not correct for spin polarized...
  myG.resize(NumParticles); // not correct for spin polarized...
  myL.resize(NumParticles); // not correct for spin polarized...
  dFa.resize(nel,norb);
  Ajk_sum.resize(nel,norb);
  Qmat.resize(nel,norb);
  Fmat.resize(nel,norb);
  Fmatdiag.resize(norb);
  psiMinv_temp.resize(NumPtcls,norb);
  psiV.resize(norb);
  psiM_temp.resize(NumPtcls,norb);
  // For forces
  /*  not used
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
  */
}

/** replace of SPOSetBase::evaluate function with the removal of t_logpsi */
void DiracDeterminantWithBackflow::evaluate_SPO(ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  Phi->evaluate_notranspose(BFTrans->QP, FirstIndex, LastIndex, psiM_temp, dlogdet, grad_grad_logdet);
  simd::transpose(psiM_temp.data(), NumOrbitals, psiM_temp.cols(), logdet.data(), NumOrbitals, logdet.cols());
}

void DiracDeterminantWithBackflow::evaluate_SPO(ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
{
  Phi->evaluate_notranspose(BFTrans->QP, FirstIndex, LastIndex, psiM_temp, dlogdet, grad_grad_logdet,grad_grad_grad_logdet);
  simd::transpose(psiM_temp.data(), NumOrbitals, psiM_temp.cols(), logdet.data(), NumOrbitals, logdet.cols());
}

void
DiracDeterminantWithBackflow::registerData(ParticleSet& P, WFBufferType& buf)
{
  if(NP == 0)
    //first time, allocate once
  {
    int norb = NumOrbitals;
    NP=P.getTotalNum();
    NumParticles=P.getTotalNum();
    dpsiM_temp.resize(NumPtcls,norb);
    grad_grad_psiM_temp.resize(NumPtcls,norb);
    dpsiV.resize(norb);
    d2psiV.resize(norb);
    grad_gradV.resize(norb);
    Fmatdiag_temp.resize(norb);
    myG_temp.resize(NP);
    myL_temp.resize(NP);
    resize(NumPtcls,NumOrbitals);
    FirstAddressOfG = &myG[0][0];
    LastAddressOfG = FirstAddressOfG + NP*DIM;
    FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
    FirstAddressOfGGG = grad_grad_psiM(0,0).begin(); //[0];
    LastAddressOfGGG = FirstAddressOfGGG + NumPtcls*norb*DIM*DIM;
    FirstAddressOfFm = &(Fmatdiag[0][0]); //[0];
    LastAddressOfFm = FirstAddressOfFm + NumOrbitals*DIM;
  }
  myG_temp=0.0;
  myL_temp=0.0;
  //ValueType x=evaluate(P,myG,myL);
  LogValue=evaluateLog(P,myG_temp,myL_temp);
  P.G += myG_temp;
  P.L += myL_temp;
  //add the data: determinant, inverse, gradient and laplacians
  buf.add(psiM.first_address(),psiM.last_address());
  buf.add(FirstAddressOfdV,LastAddressOfdV);
  buf.add(FirstAddressOfGGG,LastAddressOfGGG);
  buf.add(myL.first_address(), myL.last_address());
  buf.add(FirstAddressOfG,LastAddressOfG);
  buf.add(FirstAddressOfFm,LastAddressOfFm);
  buf.add(psiMinv.first_address(),psiMinv.last_address());
  buf.add(LogValue);
  buf.add(PhaseValue);
}

DiracDeterminantWithBackflow::RealType DiracDeterminantWithBackflow::updateBuffer(ParticleSet& P,
    WFBufferType& buf, bool fromscratch)
{
  // for now, always recalculate from scratch
  // enable from_scratch = true later
  UpdateTimer.start();
  myG_temp=0.0;
  myL_temp=0.0;
  LogValue=evaluateLog(P,myG_temp,myL_temp);
  P.G += myG_temp;
  P.L += myL_temp;
  //copy psiM to psiM_temp
  psiM_temp=psiM;
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(FirstAddressOfGGG,LastAddressOfGGG);
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(FirstAddressOfFm,LastAddressOfFm);
  buf.put(psiMinv.first_address(),psiMinv.last_address());
  buf.put(LogValue);
  buf.put(PhaseValue);
  UpdateTimer.stop();
  return LogValue;
}

void DiracDeterminantWithBackflow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.get(psiM.first_address(),psiM.last_address());
  buf.get(FirstAddressOfdV,LastAddressOfdV);
  buf.get(FirstAddressOfGGG,LastAddressOfGGG);
  buf.get(myL.first_address(), myL.last_address());
  buf.get(FirstAddressOfG,LastAddressOfG);
  buf.get(FirstAddressOfFm,LastAddressOfFm);
  buf.get(psiMinv.first_address(),psiMinv.last_address());
  buf.get(LogValue);
  buf.get(PhaseValue);
  //re-evaluate it for testing
  //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
  //need extra copy for gradient/laplacian calculations without updating it
  //psiM_temp = psiM;
  //dpsiM_temp = dpsiM;
  //d2psiM_temp = d2psiM;
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
DiracDeterminantWithBackflow::ValueType DiracDeterminantWithBackflow::ratio(ParticleSet& P, int iat)
{
  // FIX FIX FIX : code Woodbury formula
  psiM_temp=psiM;
  // either code Woodbury or do multiple single particle updates
  UpdateMode=ORB_PBYP_RATIO;
  std::vector<int>::iterator it = BFTrans->indexQP.begin();
  std::vector<int>::iterator it_end = BFTrans->indexQP.end();
  while(it != it_end)
  {
    if(*it<FirstIndex || *it>=LastIndex )
    {
      ++it;
      continue;
    }
    int jat = *it-FirstIndex;
    PosType dr = BFTrans->newQP[*it] - BFTrans->QP.R[*it];
    BFTrans->QP.makeMoveAndCheck(*it,dr);
    Phi->evaluate(BFTrans->QP, *it, psiV);
    for(int orb=0; orb<psiV.size(); orb++)
      psiM_temp(orb,jat) = psiV[orb];
    BFTrans->QP.rejectMove(*it);
    it++;
  }
  // FIX FIX FIX : code Woodbury formula
  psiMinv_temp = psiM_temp;
  // FIX FIX FIX : code Woodbury formula
  InverseTimer.start();
  RealType NewPhase;
  RealType NewLog=InvertWithLog(psiMinv_temp.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),NewPhase);
  InverseTimer.stop();
#if defined(QMC_COMPLEX)
  RealType ratioMag = std::exp(NewLog-LogValue);
  return curRatio = std::complex<OHMMS_PRECISION>(std::cos(NewPhase-PhaseValue)*ratioMag,std::sin(NewPhase-PhaseValue)*ratioMag);
#else
  return curRatio = std::cos(NewPhase-PhaseValue)*std::exp(NewLog-LogValue);
#endif
}

void DiracDeterminantWithBackflow::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evaluateRatiosAlltoOne. \n");
}

DiracDeterminantWithBackflow::GradType
DiracDeterminantWithBackflow::evalGrad(ParticleSet& P, int iat)
{
  GradType g;
  g=0.0;
  for(int j=0; j<NumPtcls; j++)
  {
    g += dot(BFTrans->Amat(iat,FirstIndex+j),Fmatdiag[j]);
  }
  return g;
}

DiracDeterminantWithBackflow::GradType
DiracDeterminantWithBackflow::evalGradSource(ParticleSet& P, ParticleSet& source,
    int iat)
{
  APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSource() \n");
  return GradType();
}

DiracDeterminantWithBackflow::GradType
DiracDeterminantWithBackflow::evalGradSourcep
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSourcep() \n");
  return GradType();
}


DiracDeterminantWithBackflow::GradType
DiracDeterminantWithBackflow::evalGradSource
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSource() \n");
  return GradType();
}

DiracDeterminantWithBackflow::ValueType
DiracDeterminantWithBackflow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  // FIX FIX FIX : code Woodbury formula
  psiM_temp=psiM;
  dpsiM_temp=dpsiM;
  UpdateMode=ORB_PBYP_PARTIAL;
  std::vector<int>::iterator it = BFTrans->indexQP.begin();
  std::vector<int>::iterator it_end = BFTrans->indexQP.end();
  ParticleSet::ParticlePos_t dr;
  while(it != it_end)
  {
    if(*it<FirstIndex || *it>=LastIndex )
    {
      ++it;
      continue;
    }
    int jat = *it-FirstIndex;
    PosType dr = BFTrans->newQP[*it] - BFTrans->QP.R[*it];
    BFTrans->QP.makeMoveAndCheck(*it,dr);
    Phi->evaluate(BFTrans->QP, *it, psiV, dpsiV, d2psiV);
    for(int orb=0; orb<psiV.size(); orb++)
      psiM_temp(orb,jat) = psiV[orb];
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM_temp.begin(jat));
    std::copy(grad_gradV.begin(),grad_gradV.end(),grad_grad_psiM_temp.begin(jat));
    BFTrans->QP.rejectMove(*it);
    it++;
  }
  // FIX FIX FIX : code Woodbury formula
  psiMinv_temp = psiM_temp;
  // FIX FIX FIX : code Woodbury formula
  InverseTimer.start();
  RealType NewPhase;
  RealType NewLog=InvertWithLog(psiMinv_temp.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),NewPhase);
  InverseTimer.stop();
  // update Fmatdiag_temp
  for(int j=0; j<NumPtcls; j++)
  {
    Fmatdiag_temp[j]=simd::dot(psiMinv_temp[j],dpsiM_temp[j],NumOrbitals);
    grad_iat += dot(BFTrans->Amat_temp(iat,FirstIndex+j),Fmatdiag_temp[j]);
  }
#if defined(QMC_COMPLEX)
  RealType ratioMag = std::exp(NewLog-LogValue);
  return curRatio = std::complex<OHMMS_PRECISION>(std::cos(NewPhase-PhaseValue)*ratioMag,std::sin(NewPhase-PhaseValue)*ratioMag);
#else
  return curRatio = std::cos(NewPhase-PhaseValue)*std::exp(NewLog-LogValue);
#endif
}

void DiracDeterminantWithBackflow::testL(ParticleSet& P)
{
  GradMatrix_t Fmat_p,Fmat_m;
  GradVector_t Fdiag_p,Fdiag_m;
  HessMatrix_t Kij,Qij; // finite difference and analytic derivative of Fmat
  typedef Tensor<RealType,OHMMS_DIM>     HessType_0;
  typedef TinyVector<RealType,DIM>       GradType_0;
  Matrix<GradType_0> Bij,dAij;
  Matrix<HessType_0> Aij_p,Aij_m;
  Fdiag_p.resize(NumOrbitals);
  Fdiag_m.resize(NumOrbitals);
  Fmat_p.resize(NumPtcls,NumOrbitals);
  Fmat_m.resize(NumPtcls,NumOrbitals);
  Kij.resize(NumParticles,NumPtcls);
  Qij.resize(NumParticles,NumPtcls);
  Bij.resize(NumParticles,NumParticles);
  dAij.resize(NumParticles,NumParticles);
  Aij_p.resize(NumParticles,NumParticles);
  Aij_m.resize(NumParticles,NumParticles);
  Fmat_p = 0.0;
  Fmat_m = 0.0;
  Kij = 0.0;
  Qij = 0.0;
  Bij=0.0;
  dAij = 0.0;
  Aij_p = 0.0;
  Aij_m = 0.0;
  P.update();
  BFTrans->evaluate(P);
  // calculate backflow matrix, 1st and 2nd derivatives
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
  {
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
    Fmatdiag[i] = Fmat(i,i);
  }
  // calculate gradients and first piece of laplacians
  GradType temp;
  ValueType temp2;
  int num = P.getTotalNum();
  for(int i=0; i<num; i++)
    for(int j=0; j<NumPtcls; j++)
      for(int a=0; a<3; a++)
        (Bij(i,j))[a] = (BFTrans->Bmat_full(i,FirstIndex+j))[a];
  for(int i=0; i<num; i++)
  {
    for(int j=0; j<NumPtcls; j++)
    {
      HessType q_j;
      q_j=0.0;
      for(int k=0; k<NumPtcls; k++)
        q_j += psiMinv(j,k)*grad_grad_psiM(j,k);
      for(int a=0; a<3; a++)
      {
        for(int b=0; b<3; b++)
        {
          ValueType& qijab = (Qij(i,j))(a,b);
          for(int k=0; k<NumPtcls; k++)
          {
            if(j==k)
            {
              for(int c=0; c<3; c++)
                qijab -= ((BFTrans->Amat(i,FirstIndex+k))(a,c)) * ( ((Fmat(j,k))[c]) * ((Fmat(k,j))[b]) - q_j(b,c) );
            }
            else
            {
              for(int c=0; c<3; c++)
                qijab -= ((BFTrans->Amat(i,FirstIndex+k))(a,c)) * ( ((Fmat(j,k))[c]) * ((Fmat(k,j))[b])  );
            }
          }
        }
      }
    }
  }
  RealType dr = 0.000001;
  for(int i=0; i<num; i++)
  {
    for(int a=0; a<3; a++)
    {
      (P.R[i])[a] += dr;
      P.update();
      BFTrans->evaluate(P);
      evaluate_SPO(psiM, dpsiM, grad_grad_psiM);
      psiMinv=psiM;
      InverseTimer.start();
      LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
      InverseTimer.stop();
      for(int j=0; j<NumPtcls; j++)
      {
        Fdiag_p[j]=simd::dot(psiMinv[j],dpsiM[j],NumOrbitals);
      }
      for(int j=0; j<NumPtcls; j++)
        for(int b=0; b<3; b++)
        {
          (Aij_p(i,j))(a,b) = (BFTrans->Amat(i,FirstIndex+j))(a,b);
        }
      (P.R[i])[a] -= 2.0*dr;
      P.update();
      BFTrans->evaluate(P);
      evaluate_SPO(psiM, dpsiM, grad_grad_psiM);
      psiMinv=psiM;
      InverseTimer.start();
      LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
      InverseTimer.stop();
      for(int j=0; j<NumPtcls; j++)
      {
        Fdiag_m[j]=simd::dot(psiMinv[j],dpsiM[j],NumOrbitals);
      }
      for(int j=0; j<NumPtcls; j++)
        for(int b=0; b<3; b++)
        {
          (Aij_m(i,j))(a,b) = (BFTrans->Amat(i,FirstIndex+j))(a,b);
        }
      (P.R[i])[a] += dr;
      P.update();
      BFTrans->evaluate(P);
      for(int j=0; j<NumPtcls; j++)
        for(int b=0; b<3; b++)
        {
          (Kij(i,j))(a,b) = ( (Fdiag_p[j])[b] - (Fdiag_m[j])[b]   )/((RealType)2.0*dr);
        }
      for(int j=0; j<NumPtcls; j++)
        for(int b=0; b<3; b++)
        {
          (dAij(i,j))[b] += ( (Aij_p(i,j))(a,b) - (Aij_m(i,j))(a,b) ) / ((RealType)2.0*dr);
        }
    }
  }
  std::cout <<"Testing derivative of Fjj in complex case: \n";
  for(int i=0; i<num; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      std::cout <<"i,j: " <<i <<" " <<j << std::endl;
      for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
        {
          std::cout <<a <<" " <<b <<" " <<((Kij(i,j))(a,b)) - ((Qij(i,j))(a,b)) <<" -- " <<((Kij(i,j))(a,b)) <<" -- "  <<((Qij(i,j))(a,b)) << std::endl;
        }
    }
  std::cout <<"Testing derivative of Aij in complex case: \n";
  for(int i=0; i<num; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      std::cout <<"i,j: " <<i <<" " <<j << std::endl;
      for(int a=0; a<3; a++)
      {
        std::cout <<a <<" " <<((dAij(i,j))[a]) - ((Bij(i,j))[a]) <<" -- " <<((dAij(i,j))[a]) <<" -- "  <<((Bij(i,j))[a]) << std::endl;
      }
    }
  std::cout.flush();
  APP_ABORT("Finished testL: Aborting \n");
}

DiracDeterminantWithBackflow::RealType
DiracDeterminantWithBackflow::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L)
{
  //testGG(P);
  //testL(P);
  // calculate backflow matrix, 1st and 2nd derivatives
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
  {
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
    Fmatdiag[i] = Fmat(i,i);
  }
  // calculate gradients and first piece of laplacians
  GradType temp;
  ValueType temp2;
  myG=0.0;
  myL=0.0;
  int num = P.getTotalNum();
  for(int i=0; i<num; i++)
  {
    temp=0.0;
    temp2=0.0;
    for(int j=0; j<NumPtcls; j++)
    {
      for(int k=0; k<OHMMS_DIM; k++)
        temp2 += BFTrans->Bmat_full(i,FirstIndex+j)[k]*Fmat(j,j)[k];
      temp  += dot(BFTrans->Amat(i,FirstIndex+j),Fmat(j,j));
      //temp2 += rcdot(BFTrans->Bmat_full(i,FirstIndex+j),Fmat(j,j));
    }
    myG[i] += temp;
    myL[i] += temp2;
  }
// NOTE: check derivatives of Fjj and Amat numerically here, the problem has to come from somewhere
  for(int j=0; j<NumPtcls; j++)
  {
    HessType q_j;
    q_j=0.0;
    for(int k=0; k<NumPtcls; k++)
      q_j += psiMinv(j,k)*grad_grad_psiM(j,k);
    for(int i=0; i<num; i++)
    {
      Tensor<RealType,OHMMS_DIM> AA = dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j));
      myL[i] += traceAtB(AA,q_j);
      //myL[i] += traceAtB(dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j)),q_j);
    }
    for(int k=0; k<NumPtcls; k++)
    {
      for(int i=0; i<num; i++)
      {
        Tensor<RealType,OHMMS_DIM> AA = dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k));
        HessType FF = outerProduct(Fmat(k,j),Fmat(j,k));
        myL[i] -= traceAtB(AA,FF);
        //myL[i] -= traceAtB(dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k)), outerProduct(Fmat(k,j),Fmat(j,k)));
      }
    }
  }
  for(int i=0; i<num; i++)
  {
    L[i] += myL[i];
    G[i] += myG[i];
  }
  return LogValue;
}


/** move was accepted, update the real container
*/
void DiracDeterminantWithBackflow::acceptMove(ParticleSet& P, int iat)
{
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  switch(UpdateMode)
  {
  case ORB_PBYP_RATIO:
    psiMinv = psiMinv_temp;
    psiM = psiM_temp;
    break;
  case ORB_PBYP_PARTIAL:
    psiMinv = psiMinv_temp;
    psiM = psiM_temp;
    dpsiM = dpsiM_temp;
    grad_grad_psiM = grad_grad_psiM_temp;
    Fmatdiag = Fmatdiag_temp;
    break;
  default:
    myG = myG_temp;
    myL = myL_temp;
    psiMinv = psiMinv_temp;
    psiM = psiM_temp;
    dpsiM = dpsiM_temp;
    grad_grad_psiM = grad_grad_psiM_temp;
    Fmatdiag = Fmatdiag_temp;
    break;
  }
  UpdateTimer.stop();
  curRatio=1.0;
}

/** move was rejected. Nothing to restore for now.
*/
void DiracDeterminantWithBackflow::restore(int iat)
{
  curRatio=1.0;
}


/** Calculate the value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
DiracDeterminantWithBackflow::ValueType
DiracDeterminantWithBackflow::evaluate(ParticleSet& P,
                                       ParticleSet::ParticleGradient_t& G,
                                       ParticleSet::ParticleLaplacian_t& L)
{
  RealType logval = evaluateLog(P, G, L);
#if defined(QMC_COMPLEX)
  RealType ratioMag = std::exp(logval);
  return std::complex<OHMMS_PRECISION>(std::cos(PhaseValue)*ratioMag,std::sin(PhaseValue)*ratioMag);
#else
  return std::cos(PhaseValue)*std::exp(logval);
#endif
}

void
DiracDeterminantWithBackflow::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
  /*  Note:
   *    Since evaluateDerivatives seems to always be called after
   *    evaluateDeltaLog, which in turn calls evaluateLog, many
   *    of the structures calculated here do not need to be calculated
   *    again. The only one that need to be called is a routine
   *    to calculate grad_grad_grad_psiM. The following structures
   *    should already be known here:
   *       -psiM
   *       -dpsiM
   *       -grad_grad_psiM
   *       -psiM_inv
   *       -Fmat
   */
  {
//       must compute if didn;t call earlier function
    evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
//       copy(psiM.begin(),psiM.end(),psiMinv.begin());
    psiMinv=psiM;
//       invert backflow matrix
    InverseTimer.start();
    LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    InverseTimer.stop();
//       calculate F matrix (gradients wrt bf coordinates)
//       could use dgemv with increments of 3*nCols
    for(int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumPtcls; j++)
      {
        Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
      }
  }
  int num = P.getTotalNum();
//mmorales: cheap trick for now
  for(int j=0; j<NumPtcls; j++)
    for(int k=0; k<NumPtcls; k++)
    {
      HessType& q_jk = Qmat(j,k);
      q_jk=0;
      for(int n=0; n<NumPtcls; n++)
        q_jk += psiMinv(j,n)*grad_grad_psiM(k,n);
      HessType& a_jk = Ajk_sum(j,k);
      a_jk=0;
      for(int n=0; n<num; n++)
        a_jk += dot(transpose(BFTrans->Amat(n,FirstIndex+j)),BFTrans->Amat(n,FirstIndex+k));
    }
  // this is a mess, there should be a better way
  // to rearrange this
  for (int pa=0; pa<BFTrans->optIndexMap.size(); ++pa)
    //for (int pa=0; pa<BFTrans->numParams; ++pa)
  {
    ValueType dpsia=0;
    Gtemp=0;
    ValueType dLa=0;
    GradType temp;
    temp=0;
    ValueType temp2(0);
    for(int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumPtcls; j++)
      {
        GradType f_a;
        f_a=0;
        PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
        for(int k=0; k<NumPtcls; k++)
        {
          f_a += (psiMinv(i,k)*dot(grad_grad_psiM(j,k),cj)
                  -  Fmat(k,j)*rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(i,k)));
        }
        dFa(i,j)=f_a;
      }
    for(int i=0; i<num; i++)
    {
      temp=0;
      for(int j=0; j<NumPtcls; j++)
        temp += (dot(BFTrans->Xmat(pa,i,FirstIndex+j),Fmat(j,j))
                 + dot(BFTrans->Amat(i,FirstIndex+j),dFa(j,j)));
      Gtemp[i] += temp;
    }
    for(int j=0; j<NumPtcls; j++)
    {
      GradType B_j;
      B_j=0;
      for(int i=0; i<num; i++)
        B_j += BFTrans->Bmat_full(i,FirstIndex+j);
      dLa += (rcdot(Fmat(j,j),BFTrans->Ymat(pa,FirstIndex+j)) +
              dot(B_j,dFa(j,j)));
      dpsia += rcdot(Fmat(j,j),BFTrans->Cmat(pa,FirstIndex+j));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      HessType a_j_prime;
      a_j_prime=0;
      for(int i=0; i<num; i++)
        a_j_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+j)) );
      HessType q_j_prime;
      q_j_prime=0;
      PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
      for(int k=0; k<NumPtcls; k++)
      {
        //tmp=0.0;
        //for(int n=0; n<NumPtcls; n++) tmp+=psiMinv(k,n)*grad_grad_psiM(j,n);
        //tmp *= dot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(j,k));
        q_j_prime += ( psiMinv(j,k)*(cj[0]*grad_grad_grad_psiM(j,k)[0]
                                     + cj[1]*grad_grad_grad_psiM(j,k)[1]
#if OHMMS_DIM==3
                                     + cj[2]*grad_grad_grad_psiM(j,k)[2]
#endif
                                    ) - rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(j,k))
                       *Qmat(k,j) );
      }
      dLa += (traceAtB(a_j_prime,Qmat(j,j)) + traceAtB(Ajk_sum(j,j),q_j_prime));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      for(int k=0; k<NumPtcls; k++)
      {
        HessType a_jk_prime;
        a_jk_prime=0;
        for(int i=0; i<num; i++)
          a_jk_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+k)) );
        dLa -= (traceAtB(a_jk_prime, outerProduct(Fmat(k,j),Fmat(j,k)))
                + traceAtB(Ajk_sum(j,k), outerProduct(dFa(k,j),Fmat(j,k))
                           + outerProduct(Fmat(k,j),dFa(j,k)) ));
      }  // k
    }   // j
    //int kk = pa; //BFTrans->optIndexMap[pa];
    int kk = BFTrans->optIndexMap[pa];
#if defined(QMC_COMPLEX)
    dlogpsi[kk]+=real(dpsia);
    dhpsioverpsi[kk] -= real(0.5*static_cast<ParticleSet::ParticleValue_t>(dLa)+Dot(P.G,Gtemp));
#else
    dlogpsi[kk]+=dpsia;
    dhpsioverpsi[kk] -= (0.5*static_cast<ParticleSet::ParticleValue_t>(dLa)+Dot(P.G,Gtemp));
#endif
  }
}

/* Used in MultiSlaterDeterminantWithBackflow::evaluateDerivatives */
void
DiracDeterminantWithBackflow::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    int offset,
    Matrix<RealType>& dlogpsi,
    Array<GradType,OHMMS_DIM>& dG,
    Matrix<RealType>& dL)
{
  /*  Note:
   *    Since evaluateDerivatives seems to always be called after
   *    evaluateDeltaLog, which in turn calls evaluateLog, many
   *    of the structures calculated here do not need to be calculated
   *    again. The only one that need to be called is a routine
   *    to calculate grad_grad_grad_psiM. The following structures
   *    should already be known here:
   *       -psiM
   *       -dpsiM
   *       -grad_grad_psiM
   *       -psiM_inv
   *       -Fmat
   */
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  //for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
  // G(iat) += Fmat(i,i);
  const ValueType ConstZero(0.0);
  int num = P.getTotalNum();
  for(int j=0; j<NumPtcls; j++)
    for(int k=0; k<NumPtcls; k++)
    {
      HessType& q_jk = Qmat(j,k);
      q_jk=ConstZero;
      for(int n=0; n<NumPtcls; n++)
        q_jk += psiMinv(j,n)*grad_grad_psiM(k,n);
      HessType& a_jk = Ajk_sum(j,k);
      a_jk=ConstZero;
      for(int n=0; n<num; n++)
        a_jk += dot(transpose(BFTrans->Amat(n,FirstIndex+j)),BFTrans->Amat(n,FirstIndex+k));
    }
  ValueType sumL = Sum(myL);
  ValueType dotG = Dot(myG,myG);
  // this is a mess, there should be a better way
  // to rearrange this
  for (int pa=0; pa<BFTrans->optIndexMap.size(); ++pa)
    //for (int pa=0; pa<BFTrans->numParams; ++pa)
  {
    ValueType dpsia=ConstZero;
    Gtemp=ConstZero;
    ValueType dLa=ConstZero;
    GradType temp;
    ValueType temp2;
    for(int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumPtcls; j++)
      {
        GradType f_a;
        PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
        for(int k=0; k<NumPtcls; k++)
        {
          f_a += (psiMinv(i,k)*dot(grad_grad_psiM(j,k),cj)
                  -  Fmat(k,j)*rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(i,k)));
        }
        dFa(i,j)=f_a;
      }
    for(int i=0; i<num; i++)
    {
      temp=ConstZero;
      for(int j=0; j<NumPtcls; j++)
        temp += (dot(BFTrans->Xmat(pa,i,FirstIndex+j),Fmat(j,j))
                 + dot(BFTrans->Amat(i,FirstIndex+j),dFa(j,j)));
      Gtemp[i] += temp;
    }
    for(int j=0; j<NumPtcls; j++)
    {
      GradType B_j;
      for(int i=0; i<num; i++)
        B_j += BFTrans->Bmat_full(i,FirstIndex+j);
      dLa += (rcdot(Fmat(j,j),BFTrans->Ymat(pa,FirstIndex+j)) +
              dot(B_j,dFa(j,j)));
      dpsia += rcdot(Fmat(j,j),BFTrans->Cmat(pa,FirstIndex+j));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      HessType a_j_prime;
      for(int i=0; i<num; i++)
        a_j_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+j)) );
      HessType q_j_prime;
      PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
      for(int k=0; k<NumPtcls; k++)
      {
        //tmp=0.0;
        //for(int n=0; n<NumPtcls; n++) tmp+=psiMinv(k,n)*grad_grad_psiM(j,n);
        //tmp *= dot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(j,k));
        q_j_prime += ( psiMinv(j,k)*(cj[0]*grad_grad_grad_psiM(j,k)[0]
                                     + cj[1]*grad_grad_grad_psiM(j,k)[1]
                                     + cj[2]*grad_grad_grad_psiM(j,k)[2])
                       - rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(j,k))
                       *Qmat(k,j) );
      }
      dLa += (traceAtB(a_j_prime,Qmat(j,j)) + traceAtB(Ajk_sum(j,j),q_j_prime));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      for(int k=0; k<NumPtcls; k++)
      {
        HessType a_jk_prime;
        for(int i=0; i<num; i++)
          a_jk_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+k)) );
        dLa -= (traceAtB(a_jk_prime, outerProduct(Fmat(k,j),Fmat(j,k)))
                + traceAtB(Ajk_sum(j,k), outerProduct(dFa(k,j),Fmat(j,k))
                           + outerProduct(Fmat(k,j),dFa(j,k)) ));
      }  // k
    }   // j
#if defined(QMC_COMPLEX)
    convert(dpsia,dlogpsi(offset,pa));
    convert(dLa + sumL*dpsia + dotG*dpsia + static_cast<ValueType>(2.0*Dot(myG,Gtemp)), dL(offset,pa));
#else
    dlogpsi(offset,pa) = dpsia; // \nabla_pa ln(D)
    dL(offset,pa) = dLa + sumL*dpsia + dotG*dpsia + static_cast<ValueType>(2.0*Dot(myG,Gtemp));
#endif
    // \sum_i (\nabla_pa  \nabla2_i D) / D
    for(int k=0; k<num; k++)
      dG(offset,pa,k) = Gtemp[k] + myG[k]*static_cast<ParticleSet::ParticleValue_t>(dpsia);  // (\nabla_pa \nabla_i D) / D
  }
}

void DiracDeterminantWithBackflow::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi,
    ParticleSet::ParticleGradient_t* G0,
    ParticleSet::ParticleLaplacian_t* L0,
    int pa )
{
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  //for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
  // G(iat) += Fmat(i,i);
  // this is a mess, there should be a better way
  // to rearrange this
  //for (int pa=0; pa<BFTrans->optIndexMap.size(); ++pa)
  //for (int pa=0; pa<BFTrans->numParams; ++pa)
  const ValueType ConstZero(0.0);
  {
    ValueType dpsia=ConstZero;
    Gtemp=ConstZero;
    //ValueType dLa=0.0;
    La1=La2=La3=ConstZero;
    GradType temp;
    ValueType temp2;
    int num = P.getTotalNum();
    for(int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumPtcls; j++)
      {
        GradType f_a;
        PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
        for(int k=0; k<NumPtcls; k++)
        {
          f_a += (psiMinv(i,k)*dot(grad_grad_psiM(j,k),cj)
                  -  Fmat(k,j)*rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(i,k)));
        }
        dFa(i,j)=f_a;
      }
    for(int i=0; i<num; i++)
    {
      temp=ConstZero;
      for(int j=0; j<NumPtcls; j++)
        temp += (dot(BFTrans->Xmat(pa,i,FirstIndex+j),Fmat(j,j))
                 + dot(BFTrans->Amat(i,FirstIndex+j),dFa(j,j)));
      Gtemp[i] += temp;
    }
    for(int j=0; j<NumPtcls; j++)
    {
      GradType B_j;
      for(int i=0; i<num; i++)
        B_j += BFTrans->Bmat_full(i,FirstIndex+j);
      La1 += (rcdot(Fmat(j,j),BFTrans->Ymat(pa,FirstIndex+j)) +
              dot(B_j,dFa(j,j)));
      dpsia += rcdot(Fmat(j,j),BFTrans->Cmat(pa,FirstIndex+j));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      HessType q_j;
      for(int k=0; k<NumPtcls; k++)
        q_j += psiMinv(j,k)*grad_grad_psiM(j,k);
// define later the dot product with a transpose tensor
      HessType a_j;
      for(int i=0; i<num; i++)
        a_j += dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j));
      HessType a_j_prime;
      for(int i=0; i<num; i++)
        a_j_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+j)) );
      HessType q_j_prime, tmp;
      PosType& cj = BFTrans->Cmat(pa,FirstIndex+j);
      for(int k=0; k<NumPtcls; k++)
      {
        tmp=ConstZero;
        for(int n=0; n<NumPtcls; n++)
          tmp+=psiMinv(k,n)*grad_grad_psiM(j,n);
        tmp *= rcdot(BFTrans->Cmat(pa,FirstIndex+k),Fmat(j,k));
        q_j_prime += ( psiMinv(j,k)*(cj[0]*grad_grad_grad_psiM(j,k)[0]
                                     + cj[1]*grad_grad_grad_psiM(j,k)[1]
                                     + cj[2]*grad_grad_grad_psiM(j,k)[2]) - tmp);
      }
      La2 += (traceAtB(a_j_prime,q_j) + traceAtB(a_j,q_j_prime));
    }
    for(int j=0; j<NumPtcls; j++)
    {
      for(int k=0; k<NumPtcls; k++)
      {
        HessType a_jk;
        for(int i=0; i<num; i++)
          a_jk += dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k));
        HessType a_jk_prime;
        for(int i=0; i<num; i++)
          a_jk_prime += ( dot(transpose(BFTrans->Xmat(pa,i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k)) + dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Xmat(pa,i,FirstIndex+k)) );
        La3 -= (traceAtB(a_jk_prime, outerProduct(Fmat(k,j),Fmat(j,k)))
                + traceAtB(a_jk, outerProduct(dFa(k,j),Fmat(j,k))
                           + outerProduct(Fmat(k,j),dFa(j,k)) ));
      }  // k
    }   // j
    int kk = pa; //BFTrans->optIndexMap[pa];
#if defined(QMC_COMPLEX)
    dlogpsi[kk]+=real(dpsia);
    dhpsioverpsi[kk] -= real(0.5*static_cast<ParticleSet::ParticleValue_t>(La1+La2+La3)+Dot(P.G,Gtemp));
#else
    dlogpsi[kk]+=dpsia;
    dhpsioverpsi[kk] -= (0.5*static_cast<ParticleSet::ParticleValue_t>(La1+La2+La3)+Dot(P.G,Gtemp));
#endif
    *G0 += Gtemp;
    (*L0)[0] += La1+La2+La3;
  }
}

OrbitalBasePtr DiracDeterminantWithBackflow::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use DiracDeterminantWithBackflow::makeClone");
  return 0;
}

DiracDeterminantWithBackflow* DiracDeterminantWithBackflow::makeCopy(SPOSetBasePtr spo) const
{
//    BackflowTransformation *BF = BFTrans->makeClone();
  // mmorales: particle set is only needed to get number of particles, so using QP set here
  DiracDeterminantWithBackflow* dclone= new DiracDeterminantWithBackflow(BFTrans->QP,spo,BFTrans);
  dclone->resize(NumPtcls, NumOrbitals);
  dclone->Optimizable=Optimizable;
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}

DiracDeterminantWithBackflow::DiracDeterminantWithBackflow(const DiracDeterminantWithBackflow& s):
  DiracDeterminantBase(s),BFTrans(s.BFTrans)

{
  registerTimers();
  this->resize(s.NumPtcls,s.NumOrbitals);
}

//SPOSetBasePtr  DiracDeterminantWithBackflow::clonePhi() const
//{
//  return Phi->makelone();
//}

void DiracDeterminantWithBackflow::testGG(ParticleSet& P)
{
  ParticleSet::ParticlePos_t qp_0;
  qp_0.resize(BFTrans->QP.getTotalNum());
  ValueMatrix_t psiM_1,psiM_2;
  ValueMatrix_t psiM_3,psiM_4;
  GradMatrix_t dpsiM_1,dpsiM_2;
  HessMatrix_t dgM, ggM, ggM0;
  psiM_1.resize(NumPtcls,NumOrbitals);
  psiM_2.resize(NumPtcls,NumOrbitals);
  psiM_3.resize(NumPtcls,NumOrbitals);
  psiM_4.resize(NumPtcls,NumOrbitals);
  dpsiM_1.resize(NumPtcls,NumOrbitals);
  dgM.resize(NumPtcls,NumOrbitals);
  ggM.resize(NumPtcls,NumOrbitals);
  ggM0.resize(NumPtcls,NumOrbitals);
  const RealType dh = 0.0000000001;//PREC_WARNING
  for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
    qp_0[i] = BFTrans->QP.R[i];
  Phi->evaluate_notranspose(BFTrans->QP, FirstIndex, LastIndex, psiM, dpsiM, ggM);
  app_log() <<"Testing GGType calculation: " << std::endl;
  for(int lx=0; lx<3; lx++)
  {
    for(int ly=0; ly<3; ly++)
    {
      if(lx == ly )
      {
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] + dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_1, dpsiM_1, ggM0);
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] - dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_2, dpsiM_1, ggM0);
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        BFTrans->QP.update();
        evaluate_SPO(psiM_3, dpsiM_1, ggM0);
        for(int i=0; i<NumPtcls; i++)
          for(int j=0; j<NumOrbitals; j++)
            (dgM(i,j))(lx,ly) = (psiM_1(i,j)+psiM_2(i,j)-(RealType)2.0*psiM_3(i,j))/(dh*dh);
      }
      else
      {
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] + dh;
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][ly] = qp_0[i][ly] + dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_1, dpsiM_1, ggM0);
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] - dh;
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][ly] = qp_0[i][ly] - dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_2, dpsiM_1, ggM0);
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] + dh;
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][ly] = qp_0[i][ly] - dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_3,dpsiM_1,ggM0);
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i] = qp_0[i];
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][lx] = qp_0[i][lx] - dh;
        for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
          BFTrans->QP.R[i][ly] = qp_0[i][ly] + dh;
        BFTrans->QP.update();
        evaluate_SPO(psiM_4, dpsiM_1, ggM0);
        for(int i=0; i<NumPtcls; i++)
          for(int j=0; j<NumOrbitals; j++)
            (dgM(i,j))(lx,ly) = (psiM_1(i,j)+psiM_2(i,j)-psiM_3(i,j)-psiM_4(i,j))/((RealType)4.0*dh*dh);
      }
    }
  }
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumOrbitals; j++)
    {
      std::cout <<"i,j: " <<i <<" " <<j << std::endl;
      for(int lx=0; lx<3; lx++)
        for(int ly=0; ly<3; ly++)
        {
          std::cout <<"a,b: " <<lx <<" " <<ly <<(dgM(i,j))(lx,ly)-(ggM(i,j))(lx,ly) <<" -- " <<(dgM(i,j))(lx,ly) <<" -- " <<(ggM(i,j))(lx,ly) << std::endl;
        }
    }
  for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
    BFTrans->QP.R[i] = qp_0[i];
  BFTrans->QP.update();
}

void DiracDeterminantWithBackflow::testGGG(ParticleSet& P)
{
  ParticleSet::ParticlePos_t qp_0;
  qp_0.resize(BFTrans->QP.getTotalNum());
  ValueMatrix_t psiM_1,psiM_2;
  GradMatrix_t dpsiM_1,dpsiM_2;
  HessMatrix_t ggM_1, ggM_2;
  psiM_1.resize(NumPtcls,NumOrbitals);
  psiM_2.resize(NumPtcls,NumOrbitals);
  dpsiM_1.resize(NumPtcls,NumOrbitals);
  dpsiM_2.resize(NumPtcls,NumOrbitals);
  ggM_1.resize(NumPtcls,NumOrbitals);
  ggM_2.resize(NumPtcls,NumOrbitals);
  GGGMatrix_t  ggg_psiM1,ggg_psiM2;
  ggg_psiM1.resize(NumPtcls,NumOrbitals);
  ggg_psiM2.resize(NumPtcls,NumOrbitals);
  const RealType dh = 0.000001; //PREC_WARNING
  for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
    qp_0[i] = BFTrans->QP.R[i];
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  app_log() <<"Testing GGGType calculation: " << std::endl;
  for(int lc=0; lc<3; lc++)
  {
    for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
      BFTrans->QP.R[i] = qp_0[i];
    for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
      BFTrans->QP.R[i][lc] = qp_0[i][lc] + dh;
    BFTrans->QP.update();
    evaluate_SPO(psiM_1, dpsiM_1, ggM_1, ggg_psiM1);
    for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
      BFTrans->QP.R[i][lc] = qp_0[i][lc] - dh;
    BFTrans->QP.update();
    evaluate_SPO(psiM_2, dpsiM_2, ggM_2, ggg_psiM2);
    const RealType dh2=RealType(0.5/dh);
    RealType maxD(0);
    ValueType av(0);
    RealType cnt(0);
    for(int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumOrbitals; j++)
      {
        //HessType dG = (ggM_1(i,j)-ggM_2(i,j))/((RealType)2.0*dh)-(grad_grad_grad_psiM(i,j))[lc];
        HessType dG = (ggM_1(i,j)-ggM_2(i,j))*dh2-grad_grad_grad_psiM(i,j)[lc];
        for(int la=0; la<9; la++)
        {
          cnt++;
          av+=dG[la];
#if defined(QMC_COMPLEX)
          if( std::abs(dG[la].real()) > maxD)
            maxD = std::abs(dG[la].real());
#else
          if( std::abs(dG[la]) > maxD)
            maxD = std::abs(dG[la]);
#endif
        }
        app_log() <<i <<"  " <<j <<"\n"
                  <<"dG : \n" <<dG << std::endl
                  <<"GGG: \n" <<(grad_grad_grad_psiM(i,j))[lc] << std::endl;
      }
    app_log() <<"lc, av, max: " <<lc <<"  " <<av/cnt <<"  "
              <<maxD << std::endl;
  }
  for(int i=0; i<BFTrans->QP.getTotalNum(); i++)
    BFTrans->QP.R[i] = qp_0[i];
  BFTrans->QP.update();
}

void DiracDeterminantWithBackflow::testDerivFjj(ParticleSet& P, int pa)
{
  app_log() <<" Testing derivatives of the F matrix, prm: " <<pa << std::endl;
  opt_variables_type wfVars,wfvar_prime;
  BFTrans->checkInVariables(wfVars);
  BFTrans->checkOutVariables(wfVars);
  int Nvars= wfVars.size();
  wfvar_prime= wfVars;
  GradMatrix_t dpsiM_1, dpsiM_2, dpsiM_0;
  dpsiM_0.resize(NumPtcls,NumOrbitals);
  dpsiM_1.resize(NumPtcls,NumOrbitals);
  dpsiM_2.resize(NumPtcls,NumOrbitals);
  const RealType dh(0.00001); //PREC_WARN
  int pr = pa;
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  wfvar_prime[pr] = wfVars[pr]+ dh;
  BFTrans->resetParameters(wfvar_prime);
  BFTrans->evaluateDerivatives(P);
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      dpsiM_1(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  wfvar_prime[pr] = wfVars[pr]- dh;
  BFTrans->resetParameters(wfvar_prime);
  BFTrans->evaluateDerivatives(P);
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      dpsiM_2(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  BFTrans->resetParameters(wfvar_prime);
  BFTrans->evaluateDerivatives(P);
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM, grad_grad_grad_psiM);
  //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
  psiMinv=psiM;
  // invert backflow matrix
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  // calculate F matrix (gradients wrt bf coordinates)
  // could use dgemv with increments of 3*nCols
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  RealType cnt=0,av=0;
  RealType maxD(0);
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
      for(int lb=0; lb<3; lb++)
      {
        dpsiM_0(i,j)[lb]=0.0;
        for(int k=0; k<NumPtcls; k++)
          for(int lc=0; lc<3; lc++)
          {
            dpsiM_0(i,j)[lb] += (BFTrans->Cmat(pr,FirstIndex+j)[lc]*psiMinv(i,k)*grad_grad_psiM(j,k)(lb,lc) - BFTrans->Cmat(pr,FirstIndex+k)[lc]*Fmat(i,k)[lc]*Fmat(k,j)[lb]);
          }
        cnt++;
#if defined(QMC_COMPLEX)
        RealType t=real(dpsiM_0(i,j)[lb]-( dpsiM_1(i,j)[lb]-dpsiM_2(i,j)[lb] ))/(2*dh);
        av+=t;
        maxD=std::max(maxD,std::abs(t));
#else
        av+=dpsiM_0(i,j)[lb]-( dpsiM_1(i,j)[lb]-dpsiM_2(i,j)[lb] )/(2*dh);
        if( std::abs(dpsiM_0(i,j)[lb]-( dpsiM_1(i,j)[lb]-dpsiM_2(i,j)[lb] )/(2*dh)) > maxD  )
          maxD=dpsiM_0(i,j)[lb]-( dpsiM_1(i,j)[lb]-dpsiM_2(i,j)[lb] )/(2*dh);
#endif
      }
  app_log() <<" av,max : " <<av/cnt <<"  " <<maxD << std::endl;
  //APP_ABORT("testing Fij \n");
}


void DiracDeterminantWithBackflow::testDerivLi(ParticleSet& P, int pa)
{
  //app_log() <<"Testing new L[i]: \n";
  opt_variables_type wfVars,wfvar_prime;
  BFTrans->checkInVariables(wfVars);
  BFTrans->checkOutVariables(wfVars);
  int Nvars= wfVars.size();
  wfvar_prime= wfVars;
  RealType dh=0.00001;
  //BFTrans->evaluate(P);
  ValueType L1a,L2a,L3a,L0a;
  ValueType L1b,L2b,L3b,L0b;
  ValueType L1c,L2c,L3c,L0c;
  //dummyEvalLi(L1,L2,L3);
  myG=0.0;
  myL=0.0;
  //ValueType ps = evaluateLog(P,myG,myL);
  //L0 = Sum(myL);
  //app_log() <<"L old, L new: " <<L0 <<"  " <<L1+L2+L3 << std::endl;
  app_log() << std::endl <<" Testing derivatives of L[i] matrix. " << std::endl;
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  wfvar_prime[pa] = wfVars[pa]+ dh;
  BFTrans->resetParameters(wfvar_prime);
  BFTrans->evaluate(P);
  dummyEvalLi(L1a,L2a,L3a);
  for (int j=0; j<Nvars; j++)
    wfvar_prime[j]=wfVars[j];
  wfvar_prime[pa] = wfVars[pa]- dh;
  BFTrans->resetParameters(wfvar_prime);
  BFTrans->evaluate(P);
  dummyEvalLi(L1b,L2b,L3b);
  BFTrans->resetParameters(wfVars);
  BFTrans->evaluateDerivatives(P);
  std::vector<RealType> dlogpsi;
  std::vector<RealType> dhpsi;
  dlogpsi.resize(Nvars);
  dhpsi.resize(Nvars);
  evaluateDerivatives(P,wfVars,dlogpsi,dhpsi,&myG,&myL,pa);
  app_log() <<"pa: " <<pa << std::endl
            <<"dL Numrical: "
            <<(L1a-L1b)/(2*dh) <<"  "
            <<(L2a-L2b)/(2*dh) <<"  "
            <<(L3a-L3b)/(2*dh) <<"\n"
            <<"dL Analitival: "
            <<La1 <<"  "
            <<La2 <<"  "
            <<La3 << std::endl
            <<" dLDiff: "
            <<(L1a-L1b)/(2*dh)-La1 <<"  "
            <<(L2a-L2b)/(2*dh)-La2 <<"  "
            <<(L3a-L3b)/(2*dh)-La3 << std::endl;
}


// evaluate \sum_i L[i] splitted into three pieces
void DiracDeterminantWithBackflow::dummyEvalLi(ValueType& L1, ValueType& L2, ValueType& L3)
{
  L1=L2=L3=0.0;
  evaluate_SPO(psiM, dpsiM, grad_grad_psiM);
  psiMinv=psiM;
  InverseTimer.start();
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
      Fmat(i,j)=simd::dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
  GradType temp;
  ValueType temp2;
  int num = BFTrans->QP.getTotalNum();
  for(int i=0; i<num; i++)
  {
    for(int j=0; j<NumPtcls; j++)
      L1 += rcdot(Fmat(j,j),BFTrans->Bmat_full(i,FirstIndex+j));
  }
  for(int j=0; j<NumPtcls; j++)
  {
    HessType q_j;
    for(int k=0; k<NumPtcls; k++)
      q_j += psiMinv(j,k)*grad_grad_psiM(j,k);
    HessType a_j;
    for(int i=0; i<num; i++)
      a_j += dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+j));
    L2 += traceAtB(a_j,q_j);
    for(int k=0; k<NumPtcls; k++)
    {
      HessType a_jk;
      for(int i=0; i<num; i++)
        a_jk += dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k));
      L3 -= traceAtB(a_jk, outerProduct(Fmat(k,j),Fmat(j,k)));
    }
  }
}

}
