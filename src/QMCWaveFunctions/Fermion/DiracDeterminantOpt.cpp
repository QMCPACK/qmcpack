//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminantOpt.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{

DiracDeterminantOpt::DiracDeterminantOpt
(ParticleSet &ptcl, SPOSetBasePtr const &gs_spos, int first) :
  DiracDeterminantBase(gs_spos, first)
{
  targetPtcl = &ptcl;
  NumOrbitals = gs_spos->OrbitalSetSize;
  NumBasis    = gs_spos->BasisSetSize;
  BasisVals.resize(NumOrbitals,NumBasis);
  BasisGrad.resize(NumOrbitals,NumBasis);
  BasisLapl.resize(NumOrbitals,NumBasis);
  dlogdet_dC.resize(NumOrbitals, NumBasis);
  G_gamma.resize(NumOrbitals, NumBasis);
  L_gamma.resize(NumOrbitals, NumBasis);
  dgrad_dC.resize(NumOrbitals,NumBasis);
  dlapl_dC.resize(NumOrbitals,NumBasis);
  Gamma.resize(NumOrbitals,NumBasis);
  MyG.resize(NumOrbitals);
  Optimizable = true;
}

DiracDeterminantBase*
DiracDeterminantOpt::makeCopy(SPOSetBasePtr spo) const
{
  DiracDeterminantBase* dclone= new DiracDeterminantOpt(*targetPtcl, spo, FirstIndex);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}



// Note:  Currently, this calls Phi-evaluate.  This should not be
// necessary if the GS orbitals and basis orbitals are cacheed.
void
DiracDeterminantOpt::resetParameters(const opt_variables_type& optvars)
{
  Phi->resetParameters(optvars);
  // Update the direct matrices
  Phi->evaluate_notranspose(*targetPtcl, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
  // Invert PsiM
  if(NumPtcls==1)
    psiM(0,0) = (RealType)1.0/psiM_temp(0,0);
  else
  {
    InverseTimer.start();
    invertPsiM(psiM_temp,psiM);
    InverseTimer.stop();
  }
  psiM_temp = psiM;
}

void
DiracDeterminantOpt::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
  // The dlogpsi part is simple -- just ratios
  // First, evaluate the basis functions
  // std::cerr << "GEMM 1:\n";
  // fprintf (stderr, "FirstIndex = %d  LastIndex=%d\n", FirstIndex,
  // LastIndex);
  resetParameters(active);
  Phi->evaluateBasis (P, FirstIndex, LastIndex, BasisVals, BasisGrad, BasisLapl);
  BLAS::gemm ('N', 'T', NumBasis, NumOrbitals, NumOrbitals, 1.0,
              BasisVals.data(), NumBasis, psiM.data(), NumOrbitals,
              0.0, dlogdet_dC.data(), NumBasis);
//     for (int i=0; i<NumOrbitals; i++)
//       for (int j=0; j<NumBasis; j++) {
// 	dlogdet_dC(i,j) = 0.0;
// 	for (int n=0; n<NumOrbitals; n++) {
// 	  dlogdet_dC(i,j) += psiM(n,i) * BasisVals(n,j);
// 	  // fprintf (stderr, "BasisVals(%d,%d) = %12.6e\n", n, j, BasisVals(n,j));
// 	}
//       }
  Gamma = dlogdet_dC;
//     for (int n=0; n<NumOrbitals; n++)
//       for (int j=0; j<NumBasis; j++) {
// 	Gamma(n,j) = 0.0;
// 	for (int k=0; k<NumOrbitals; k++) {
// 	  Gamma(n,j) += psiM(k,n) * BasisVals(k,j);
// 	  // fprintf (stderr, "BasisVals(%d,%d) = %12.6e\n", n, j, BasisVals(n,j));
// 	}
//       }
  // Now, d_dC should hold d/dC_{ij} log(det).
  // Multiply L matrix by gamma matrix, as shown in eq. 17 of
  // docs/OrbitalOptimization.tex
  L_gamma = BasisLapl;
  BLAS::gemm ('N', 'N', NumBasis, NumOrbitals, NumOrbitals, -1.0,
              Gamma.data(), NumBasis, d2psiM.data(), NumOrbitals,
              1.0, L_gamma.data(), NumBasis);
//     for (int l=0; l<NumOrbitals; l++)
//       for (int j=0; j<NumBasis; j++) {
// 	// L_gamma(l,j) = BasisLapl(l,j);
// 	L_gamma(l,j) = BasisLapl(l,j);
// 	G_gamma(l,j) = BasisGrad(l,j);
// 	for (int n=0; n<NumOrbitals; n++) {
// 	  L_gamma(l,j) -= d2psiM(l,n)*Gamma(n,j);
// 	  G_gamma(l,j) -=  dpsiM(l,n)*Gamma(n,j);
// 	}
//       }
  for (int l=0; l<NumOrbitals; l++)
    for (int j=0; j<NumBasis; j++)
    {
      G_gamma(l,j) = BasisGrad(l,j);
      for (int n=0; n<NumOrbitals; n++)
        G_gamma(l,j) -=  dpsiM(l,n)*Gamma(n,j);
    }
  // Now, compute d/dC_{ij} lapl(det)/det by multiplying by Ainv
  BLAS::gemm('N', 'T', NumBasis, NumOrbitals, NumOrbitals, -0.5,
             L_gamma.data(), NumBasis, psiM.data(), NumOrbitals,
             0.0, dlapl_dC.data(), NumBasis);
  for (int ptcl=0; ptcl<NumOrbitals; ptcl++)
  {
    MyG[ptcl] = PosType();
    for (int orb=0; orb<NumOrbitals; orb++)
      MyG[ptcl] += dpsiM(ptcl,orb)*psiM(ptcl,orb);
//       fprintf (stderr, "myG  = %11.4e %11.4e %11.4e\n",
// 	       myG[ptcl][0],  myG[ptcl][1],  myG[ptcl][2]);
//       fprintf (stderr, "MyG = %11.4e %11.4e %11.4e\n",
// 	       MyG[ptcl][0], MyG[ptcl][1], MyG[ptcl][2]);
//       fprintf (stderr, "P.G  = %11.4e %11.4e %11.4e\n",
// 	       P.G[ptcl+FirstIndex][0],
// 	       P.G[ptcl+FirstIndex][1],
// 	       P.G[ptcl+FirstIndex][2]);
  }
  for (int i=0; i<NumOrbitals; i++)
    for (int j=0; j<NumBasis; j++)
      for (int l=0; l<NumOrbitals; l++)
      {
        GradType g = P.G[FirstIndex+l];
        GradType dg1 = psiM(l,i)*G_gamma(l,j);
        GradType dg2 = g*dlogdet_dC(i,j);
        GradType dg = dg1;// - dg2;
// 	  fprintf (stderr, "dg = %9.4f %9.4f %9.4f\n",
// 		   dg[0], dg[1], dg[2]);
// 	  fprintf (stderr, "g1 = %11.4e %11.4e %11.4e  g2 = %11.4e %11.4e %11.4e\n",
// 		   dg1[0], dg1[1], dg1[2], dg2[0], dg2[1], dg2[2]);
        g -= MyG[l];
        dlapl_dC(i,j) -= dot(g, dg);
      }
//     for (int i=0; i<NumOrbitals; i++)
//       for (int j=0; j<NumBasis; j++) {
// 	dlapl_dC(i,j) = 0.0;
// 	for (int l=0; l<NumOrbitals; l++) {
// 	  // dlapl_dC(i,j) += psiM(n,i)*L_gamma(n,j);
// 	  // dlapl_dC(i,j) += psiM(l,i)*L_gamma(l,j);
// 	  dlapl_dC(i,j) -= 0.5*psiM(l,i)*L_gamma(l,j);
// 	  //dlapl_dC(i,j) += psiM(l,i)*BasisLapl(l,j);
// 	}
//       }
  // Pull elements from dense d_dC matrices and put into parameter
  // derivatives, dlogpsi and dhpsioverpsi
  Phi->copyParamsFromMatrix(active, dlogdet_dC, dlogpsi);
  Phi->copyParamsFromMatrix(active,   dlapl_dC, dhpsioverpsi);
}

void
DiracDeterminantOpt::checkInVariables(opt_variables_type& active)
{
  Phi->checkInVariables(active);
}


void
DiracDeterminantOpt::checkOutVariables(const opt_variables_type& active)
{
  Phi->checkOutVariables(active);
}

}
