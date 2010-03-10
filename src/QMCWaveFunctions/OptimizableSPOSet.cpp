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

#include "QMCWaveFunctions/OptimizableSPOSet.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

  bool
  OptimizableSPOSet::put (xmlNodePtr node)
  {

  }

  void 
  OptimizableSPOSet::checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }
  
  void 
  OptimizableSPOSet::checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }


  void
  OptimizableSPOSet::resetParameters(const opt_variables_type& active)
  {
    for (int i=0; i<ParamPointers.size(); i++) {
      int loc=myVars.where(i);
      if (loc>=0) *(ParamPointers[i])=myVars[i]=active[loc];
    }
  }

  // Obsolete
  void
  OptimizableSPOSet::evaluateDerivatives
  (ParticleSet& P, int iat, const opt_variables_type& active,
   ValueMatrix_t& d_phi, ValueMatrix_t& d_lapl_phi)
  {
    for (int i=0; i<d_phi.size(); i++)
      for (int j=0; j<N; j++)
	d_phi[i][j] = d_lapl_phi[i][j] = ValueType();
    
    // Evaluate basis states
    if (BasisOrbitals) {
      BasisOrbitals->evaluate(P,iat,BasisVal);
      vector<TinyVector<int,2> >::iterator iter;
      vector<TinyVector<int,2> >& act = ActiveBasis[iat];
      
      for (iter=act.begin(); iter != act.end(); iter++) {
	int elem  = (*iter)[0];
	int param = (*iter)[1];
	int loc = myVars.where(param);
	if (loc >= 0);
      }

    }
    else {

    }
  }

  void
  OptimizableSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    GSOrbitals->evaluate(P,iat,GSVal);
    if (BasisOrbitals) {
      BasisOrbitals->evaluate(P,iat,BasisVal);
      BLAS::gemv (N, M, C.data(), &(GSVal[N]), &(psi[0]));
    }
    else 
      BLAS::gemv (N, M, C.data(), &(GSVal[N]), &(psi[0]));

    for (int i=0; i<N; i++)	psi[i] += GSVal[i];
  }

  void
  OptimizableSPOSet::evaluate(const ParticleSet& P, PosType r, 
			      vector<RealType> &psi)
  {

  }
  
  void
  OptimizableSPOSet::evaluate(const ParticleSet& P, int iat, 
			      ValueVector_t& psi, GradVector_t& dpsi, 
			      ValueVector_t& d2psi)
  {
    GSOrbitals->evaluate(P,iat,GSVal,GSGrad,GSLapl);
    if (BasisOrbitals) {
      BasisOrbitals->evaluate(P,iat,BasisVal,BasisGrad,BasisLapl);
      BLAS::gemv (N, M, C.data(), &(GSVal[N]),  &(psi[0]));
      BLAS::gemv (N, M, C.data(), &(GSLapl[N]), &(d2psi[0]));
    }
    else {
      BLAS::gemv (N, M, C.data(), &(GSVal[N]),  &(psi[0]));
      BLAS::gemv (N, M, C.data(), &(GSLapl[N]), &(d2psi[0]));
    }
    
    for (int i=0; i<N; i++) {
      psi[i]  += GSVal[i];
      d2psi[i] += GSLapl[i];
    }
  }


  void 
  OptimizableSPOSet::evaluateBasis (const ParticleSet &P, int first, int last,
				    ValueMatrix_t &basis_val, GradMatrix_t &basis_grad,
				    ValueMatrix_t &basis_lapl)
  {
    if (BasisOrbitals) 
      BasisOrbitals->evaluate_notranspose(P, first, last, basis_val, basis_grad, basis_lapl);
    else {
      for (int iat=first; iat<last; iat++) {
	GSOrbitals->evaluate (P, iat, GSVal, GSGrad, GSLapl);
	for (int i=0; i<M; i++) {
	  basis_val (iat,i) = GSVal[N+i];
      	  basis_grad(iat,i) = GSGrad[N+i];
      	  basis_lapl(iat,i) = GSLapl[N+i];
	}
      }
    }
      
  }

  void 
  OptimizableSPOSet::copyParamsFromMatrix (const opt_variables_type& active,
					   const Matrix<RealType> &mat,
					   vector<RealType> &destVec)
  {
    for (int ip=0; ip<myVars.size(); ip++) {
      int loc = myVars.where(ip);
      if (loc >= 0) {
	TinyVector<int,2> idx = ParamIndex[ip];
	destVec[loc] = mat(idx[0], idx[1]);
      }
    }
  }

  void 
  OptimizableSPOSet::copyParamsFromMatrix (const opt_variables_type& active,
					   const Matrix<ComplexType> &mat,
					   vector<RealType> &destVec)
  {
    for (int ip=0; ip<myVars.size(); ip+=2) {
      int loc = myVars.where(ip);
      if (loc >= 0) {
	TinyVector<int,2> idx = ParamIndex[ip];
	assert (ParamIndex[ip+1][0] == idx[0] &&
		ParamIndex[ip+1][1] == idx[1]);
	destVec[loc] = mat(idx[0], idx[1]).real();
	loc = myVars.where(ip+1);
	assert (loc >= 0);
	destVec[loc] = mat(idx[0], idx[1]).imag();
      }
    }
  }



  void
  OptimizableSPOSet::evaluate_notranspose
  (const ParticleSet& P, int first, int last,
   ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    GSOrbitals->evaluate_notranspose
      (P, first, last, GSValMatrix, GSGradMatrix, GSLaplMatrix);
    if (BasisOrbitals) {
      BasisOrbitals->evaluate_notranspose
	(P, first, last, BasisValMatrix, BasisGradMatrix, BasisLaplMatrix);
      BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
		   M, BasisValMatrix.data(), M, 0.0, logdet.data(), N);
      logdet += GSValMatrix;
      BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
		   M, BasisLaplMatrix.data(), M, 0.0, d2logdet.data(), N);
      d2logdet += GSLaplMatrix;

      // Gradient part.  
      for (int dim=0; dim<OHMMS_DIM; dim++) {
	for (int i=0; i<M; i++)
	  for (int j=0; j<N; j++)
	    GradTmpSrc(i,j) = BasisGradMatrix(i,j)[dim];
	BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(), M, 
		     GradTmpSrc.data(), M, 0.0, GradTmpDest.data(), N);
	for (int i=0; i<N; i++)
	  for (int j=0; j<N; j++)
	    dlogdet(i,j)[dim] = GradTmpDest(i,j) + GSGradMatrix(i,j)[dim];
      }
    }
    else {
      BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
		  M, GSValMatrix.data()+N, M+N, 0.0, logdet.data(), N);
      logdet += GSValMatrix;
      BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
		  M, GSLaplMatrix.data()+N, M+N, 0.0, d2logdet.data(), N);
      d2logdet += GSLaplMatrix;

      // Gradient part.  
      for (int dim=0; dim<OHMMS_DIM; dim++) {
	for (int i=0; i<M; i++)
	  for (int j=0; j<N; j++)
	    GradTmpSrc(i,j) = GSGradMatrix(i,j+N)[dim];
	BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(), M, 
		     GradTmpSrc.data(), M, 0.0, GradTmpDest.data(), N);
	for (int i=0; i<N; i++)
	  for (int j=0; j<N; j++)
	    dlogdet(i,j)[dim] = GradTmpDest(i,j) + GSGradMatrix(i,j)[dim];
      }
    }
  }

  SPOSetBase*
  OptimizableSPOSet::makeClone() const
  {

  }


}
