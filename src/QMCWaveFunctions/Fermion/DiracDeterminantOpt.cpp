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
    L_gamma.resize(NumOrbitals, NumBasis);
    dlapl_dC.resize(NumOrbitals,NumBasis);
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
    // cerr << "NumOrbitals = " << NumOrbitals << endl;
    // cerr << "NumBasis    = " << NumBasis << endl;
    
    // The dlogpsi part is simple -- just ratios
    // First, evaluate the basis functions
    // cerr << "GEMM 1:\n";
    // fprintf (stderr, "FirstIndex = %d  LastIndex=%d\n", FirstIndex, LastIndex);
    Phi->evaluateBasis (P, FirstIndex, LastIndex, BasisVals, BasisGrad, BasisLapl);
    // BLAS::gemm ('N', 'T', NumOrbitals, NumBasis, NumOrbitals, 1.0, psiM.data(), 
    // 		NumOrbitals, BasisVals.data(), NumBasis, 0.0, dlogdet_dC.data(), NumOrbitals);

    // BLAS::gemm ('T', 'N', NumBasis, NumOrbitals, NumBasis, 1.0, 
    // 		BasisVals.data(), NumBasis, psiM.data(), NumOrbitals, 0.0, dlogdet_dC.data(), NumOrbitals);
    for (int i=0; i<NumOrbitals; i++)
      for (int j=0; j<NumBasis; j++) {
	dlogdet_dC(i,j) = 0.0;
	for (int n=0; n<NumOrbitals; n++) {
	  dlogdet_dC(i,j) += psiM(n,i) * BasisVals(n,j);
	  //	  fprintf (stderr, "BasisVals(%d,%d) = %12.6e\n", n, j, BasisVals(n,j));
	}
      }

    // Now, d_dC should hold d/dC_{ij} log(det).
    
    // Multiply L matrix by gamma matrix, as shown in eq. 17 of
    // docs/OrbitalOptimization.tex 
    // cerr << "GEMM 2:\n";
    // BLAS::gemm('T', 'T', NumOrbitals, NumBasis, NumOrbitals, -1.0, d2psiM.data(),
    // 	       NumOrbitals, dlogdet_dC.data(), NumBasis, 0.0, L_gamma.data(), NumBasis);

    // BLAS::gemm('T', 'T', NumBasis, NumOrbitals, NumOrbitals, -1.0, dlogdet_dC.data(), NumOrbitals, 
    // 	       d2psiM.data(), NumOrbitals, 0.0, L_gamma.data(), NumBasis);
    for (int n=0; n<NumOrbitals; n++)
      for (int j=0; j<NumBasis; j++) {
	L_gamma(n,j) = 0.0;
	for (int m=0; m<NumOrbitals; m++)
	  L_gamma(n,j) -= d2psiM(n,m)*dlogdet_dC(m,j);
      }
    

    // Add on BasisLapl matrix
    L_gamma += BasisLapl;

    // Now, compute d/dC_{ij} lapl(det)/det by multiplying by Ainv
    // cerr << "GEMM 3:\n";
    // BLAS::gemm ('T', 'N', NumBasis, NumOrbitals, NumOrbitals, 1.0, L_gamma.data(), NumOrbitals, 
    // 		psiM.data(), NumOrbitals, 0.0, dlapl_dC.data(), NumBasis);
    for (int i=0; i<NumOrbitals; i++)
      for (int j=0; j<NumBasis; j++) {
	dlapl_dC(i,j) = 0.0;
	for (int n=0; n<NumOrbitals; n++)
	  dlapl_dC(i,j) += psiM(n,i)*L_gamma(n,j);
      }


    // Pull elements from dense d_dC matrices and put into parameter
    // derivatives, dlogpsi and dhpsioverpsi    
    Phi->copyParamsFromMatrix(active, dlogdet_dC, dlogpsi);
    Phi->copyParamsFromMatrix(active,   dlapl_dC, dhpsioverpsi);
    // fprintf (stderr, "dlogpsi:\n");
    // for (int i=0; i<dlogpsi.size(); i++)
    //   fprintf (stderr, "%12.6e ", dlogpsi[i]);
    // fprintf (stderr, "\n");
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
