//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/Fermion/DiracDeterminantOpt.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{

  DiracDeterminantOpt::DiracDeterminantOpt
  (SPOSetBasePtr const &gs_spos, int first) :
    DiracDeterminantBase(gs_spos, first)
  {
    
  }

  
  // Note:  Currently, this calls Phi-evaluate.  This should not be
  // necessary if the GS orbitals and basis orbitals are cacheed.
  void
  DiracDeterminantOpt::resetParameters(const opt_variables_type& optvars)
  {
    Phi->resetParameters(optvars);
    // Update the direct matrices
    Phi->evaluate(*targetPtcl, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);

    // Invert PsiM
    if(NumPtcls==1) 
      psiM(0,0) = 1.0/psiM(0,0);
    else {
      InverseTimer.start();
      LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
			     WorkSpace.data(),Pivot.data(),PhaseValue);
      InverseTimer.stop();
    }
    psiM_temp = psiM;
  }

  void
  DiracDeterminantOpt::evaluateDerivatives(ParticleSet& P,
					   const opt_variables_type& active,
					   vector<RealType>& dlogpsi,
					   vector<RealType>& dhpsioverpsi)
  {
    // The dlogpsi part is simple -- just ratios
    // First, evaluate the basis functions
    Phi->evaluateBasis (P, FirstIndex, LastIndex, BasisVals, BasisGrad, BasisLapl);
    BLAS::gemm ('N', 'T', NumOrbitals, NumBasis, NumOrbitals, 1.0, psiM.data(), 
		NumOrbitals, BasisVals.data(), NumBasis, 0.0, dlogdet_dC.data(), NumBasis);
    // Now, d_dC should hold d/dC_{ij} log(det).
    
    // Multiply L matrix by gamma matrix, as shown in eq. 17 of
    // docs/OrbitalOptimization.tex 
    BLAS::gemm('T', 'T', NumOrbitals, NumBasis, NumOrbitals, -1.0, d2psiM.data(),
	       NumOrbitals, dlogdet_dC.data(), NumBasis, 0.0, L_gamma.data(), NumBasis);
    
    // Add on BasisLapl matrix
    L_gamma += BasisLapl;

    // Now, compute d/dC_{ij} lapl(det)/det by multiplying by Ainv
    BLAS::gemm ('T', 'T', NumOrbitals, NumBasis, NumOrbitals, 1.0, psiM.data(), 
		NumOrbitals, L_gamma.data(), NumBasis, 0.0, dlapl_dC.data(), NumBasis);

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
