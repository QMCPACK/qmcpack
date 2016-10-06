//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminantAFM.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{

DiracDeterminantAFM::DiracDeterminantAFM
(ParticleSet &ptcl, SPOSetBasePtr const &gs_spos, int first) :
  DiracDeterminantBase(gs_spos, first)
{
  targetPtcl = &ptcl;
  NumOrbitals = gs_spos->OrbitalSetSize;
  resize(NumOrbitals,NumOrbitals);
  BasisVals.resize(NumOrbitals,NumOrbitals);
  BasisGrad.resize(NumOrbitals,NumOrbitals);
  BasisLapl.resize(NumOrbitals,NumOrbitals);
  dlogdet_dC.resize(NumOrbitals, NumOrbitals);
  G_gamma.resize(NumOrbitals, NumOrbitals);
  L_gamma.resize(NumOrbitals, NumOrbitals);
  dgrad_dC.resize(NumOrbitals,NumOrbitals);
  dlapl_dC.resize(NumOrbitals,NumOrbitals);
  Gamma.resize(NumOrbitals,NumOrbitals);
  MyG.resize(NumOrbitals);
  Optimizable = gs_spos->Optimizable;
//    if (first==0) Phi->setpm(+1);
//    else Phi->setpm(-1);
}

DiracDeterminantBase*
DiracDeterminantAFM::makeCopy(SPOSetBasePtr spo) const
{
  DiracDeterminantBase* dclone= new DiracDeterminantAFM(*targetPtcl, spo, FirstIndex);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}



//   // Note:  Currently, this calls Phi-evaluate.  This should not be
//   // necessary if the GS orbitals and basis orbitals are cacheed.
//   void
//   DiracDeterminantAFM::resetParameters(const opt_variables_type& optvars)
//   {
//     if(Optimizable) Phi->resetParameters(optvars);
//     // Update the direct matrices
//
//     Phi->evaluate(*targetPtcl, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
//
//     // Invert PsiM
//     if(NumPtcls==1)
//       psiM(0,0) = 1.0/psiM(0,0);
//     else {
//       InverseTimer.start();
//       LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
// 			     WorkSpace.data(),Pivot.data(),PhaseValue);
//       InverseTimer.stop();
//     }
//     psiM_temp = psiM;
//   }

void
DiracDeterminantAFM::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
  if(!Optimizable)
    return;
  resetParameters(active);
  int loc=Phi->myVars.where(0);
  Phi->evaluateForDeriv(P, FirstIndex, LastIndex, BasisVals, BasisGrad, BasisLapl);
  BLAS::gemm ('N', 'T', NumOrbitals, NumOrbitals, NumOrbitals, 1.0,
              BasisVals.data(), NumOrbitals, psiM.data(), NumOrbitals,
              0.0, dlogdet_dC.data(), NumOrbitals);
#ifdef QMC_COMPLEX
  for (int i=0; i<NumOrbitals; i++)
    dlogpsi[loc]+=dlogdet_dC(i,i).real();
#else
  for (int i=0; i<NumOrbitals; i++)
    dlogpsi[loc]+=dlogdet_dC(i,i);
#endif
  L_gamma = BasisLapl;
  BLAS::gemm ('N', 'N', NumOrbitals, NumOrbitals, NumOrbitals, -1.0,
              dlogdet_dC.data(), NumOrbitals, d2psiM.data(), NumOrbitals,
              1.0, L_gamma.data(), NumOrbitals);
  for (int l=0; l<NumOrbitals; l++)
    for (int j=0; j<NumOrbitals; j++)
    {
      G_gamma(l,j) = BasisGrad(l,j);
      for (int n=0; n<NumOrbitals; n++)
        G_gamma(l,j) -=  dpsiM(l,n)*dlogdet_dC(n,j);
    }
  dlapl_dC=0;
  for (int i=0; i<NumOrbitals; i++)
    for (int j=0; j<NumOrbitals; j++)
      dlapl_dC(i,i) += -(RealType)0.5*L_gamma(i,j)*psiM(i,j);
  for (int ptcl=0; ptcl<NumOrbitals; ptcl++)
  {
    MyG[ptcl] = PosType();
    for (int orb=0; orb<NumOrbitals; orb++)
      MyG[ptcl] += dpsiM(ptcl,orb)*psiM(ptcl,orb);
  }
  for (int i=0; i<NumOrbitals; i++)
    for (int l=0; l<NumOrbitals; l++)
    {
      GradType g = P.G[FirstIndex+l]-MyG[l];
      GradType dg = psiM(l,i)*G_gamma(l,i);
      dlapl_dC(i,i) -= dot(g, dg);
    }
#ifdef QMC_COMPLEX
  for (int i=0; i<NumOrbitals; i++)
    dhpsioverpsi[loc] += dlapl_dC(i,i).real();
#else
  for (int i=0; i<NumOrbitals; i++)
    dhpsioverpsi[loc] += dlapl_dC(i,i);
#endif
}

void
DiracDeterminantAFM::checkInVariables(opt_variables_type& active)
{
  if(Optimizable)
    Phi->checkInVariables(active);
}


void
DiracDeterminantAFM::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable)
    Phi->checkOutVariables(active);
}

}
