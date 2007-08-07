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
    OrbitalSetSize = norbs;
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int iat, 
				  ValueVector_t& psi)
  {
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    for(int j=0; j<OrbitalSetSize; j++) 
      Orbitals[j].evaluate(ru, psi[j]); 
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int iat, 
				  ValueVector_t& psi, GradVector_t& dpsi, 
				  ValueVector_t& d2psi)
  {
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ValueType val;
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
    for(int j=0; j<OrbitalSetSize; j++) {
      Orbitals[j].evaluate(ru, val, gu, hess);
      psi[j]   = val;
      dpsi[j]  = dot(PrimLattice.G, gu);
      d2psi[j] = trace(hess, GGt);
    }
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int first, int last,
				  ValueMatrix_t& vals, GradMatrix_t& grads, 
				  ValueMatrix_t& lapls)
  {
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType ru(PrimLattice.toUnit(P.R[iat]));
      ValueType val;
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
      for(int j=0; j<OrbitalSetSize; j++) {
	Orbitals[j].evaluate(ru, val, gu, hess);
	vals(j,i)  = val;
	grads(i,j) = dot(PrimLattice.G, gu);
	lapls(i,j) = trace(hess, GGt);
      }
    }
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
