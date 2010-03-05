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

namespace qmcplusplus
{

  DiracDeterminantOpt::DiracDeterminantOpt
  (SPOSetBasePtr const &gs_spos, 
   SPOSetBasePtr const &basis_spos, int first) :
    DiracDeterminantBase(gs_spos, first),
    Basis(basis_spos)
  {
    
  }

  
  void
  DiracDeterminantOpt::resetParameters(const opt_variables_type& optvars)
  {

  }

  void
  DiracDeterminantOpt::evaluateDerivatives(ParticleSet& P,
					   const opt_variables_type& active,
					   vector<RealType>& dlogpsi,
					   vector<RealType>& dhpsioverpsi)
  {
    // The dlogpsi part is simple -- just ratios
    // First, evaluate the basis functions
    


  }

  void
  DiracDeterminantOpt::checkInVariables(opt_variables_type& active)
  {
  }


  void
  DiracDeterminantOpt::checkOutVariables(const opt_variables_type& active)
  {
  }

}
