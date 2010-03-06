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

  void
  OptimizableSPOSet::evaluateDerivatives
  (ParticleSet& P,  const opt_variables_type& active,
   vector<RealType>& d_phi, vector<RealType>& d_lapl_phi)
  {

  }

  void
  OptimizableSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {

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

  }

  void
  OptimizableSPOSet::evaluate(const ParticleSet& P, int first, int last,
			      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
			      ValueMatrix_t& d2logdet)
  {

  }

  void
  OptimizableSPOSet::evaluate_notranspose
  (const ParticleSet& P, int first, int last,
   ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {

  }

  SPOSetBase*
  OptimizableSPOSet::makeClone() const
  {

  }


}
