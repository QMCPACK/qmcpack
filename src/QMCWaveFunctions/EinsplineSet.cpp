//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)             //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/EinsplineSet.h"

namespace qmcplusplus {
  EinsplineSetBase::UnitCellType
  EinsplineSetBase::GetLattice()
  {
    
  }
  
  void
  EinsplineSetBase::resetParameters(VarRegistry<RealType>& vlist)
  {
  }
  
  void
  EinsplineSetBase::resetTargetParticleSet(ParticleSet& e)
  {
  }
  
  void
  EinsplineSetBase:: setOrbitalSetSize(int norbs)
  {
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int iat, 
				  ValueVector_t& psi)
  {
    
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int iat, 
				  ValueVector_t& psi, GradVector_t& dpsi, 
				  ValueVector_t& d2psi)
  {
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int first, int last,
				  ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
				  ValueMatrix_t& d2logdet)
  {
    
  }
  
  void 
  EinsplineSetLocalized::evaluate (const ParticleSet& P, int iat, 
				   ValueVector_t& psi)
  {
    
  }
  
  void 
  EinsplineSetLocalized::evaluate (const ParticleSet& P, int iat, 
				   ValueVector_t& psi, GradVector_t& dpsi, 
				   ValueVector_t& d2psi)
  {

  }
    
  void 
  EinsplineSetLocalized::evaluate (const ParticleSet& P, int first, int last,
				   ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
				   ValueMatrix_t& d2logdet)
  {
    
  }
}
