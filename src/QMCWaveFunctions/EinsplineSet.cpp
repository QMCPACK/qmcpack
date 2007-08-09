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
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
//     fprintf (stderr, "  1:  ru = (%8.5f, %8.5f, %8.5f)\n",
// 	     ru[0], ru[1], ru[2]);
    for(int j=0; j<OrbitalSetSize; j++) {
      Orbitals[j].evaluate(ru, psi[j]); 
      double s,c;
      double phase = -dot(P.R[iat], Orbitals[j].kVec);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      psi *= e_mikr;
    }
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int iat, 
				  ValueVector_t& psi, GradVector_t& dpsi, 
				  ValueVector_t& d2psi)
  {
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
//     fprintf (stderr, "  2:  ru = (%8.5f, %8.5f, %8.5f)\n",
// 	     ru[0], ru[1], ru[2]);
    ValueType val;
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
    complex<double> eye (0.0, 1.0);
    for(int j=0; j<OrbitalSetSize; j++) {
      Orbitals[j].evaluate(ru, val, gu, hess);
      complex<double> u(val);
      TinyVector<complex<double>,3> gradu;
      complex<double> laplu;
      u     = val;
      gradu = dot(PrimLattice.G, gu);
      laplu = trace(hess, GGt);

      PosType k = Orbitals[j].kVec;
      double s,c;
      double phase = -dot(P.R[iat], k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      psi[j]   = e_mikr * u;
      dpsi[j]  = e_mikr*(-eye * k * u + gradu);
      d2psi[j] = e_mikr*(dot(k,k)*u - 2.0*eye*dot(k,gradu) + laplu);
    }
  }
  
  void 
  EinsplineSetExtended::evaluate (const ParticleSet& P, int first, int last,
				  ValueMatrix_t& vals, GradMatrix_t& grads, 
				  ValueMatrix_t& lapls)
  {
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType ru(PrimLattice.toUnit(P.R[iat]));
//       fprintf (stderr, "  3:  ru = (%8.5f, %8.5f, %8.5f)\n",
// 	       ru[0], ru[1], ru[2]);
      ru[0] -= std::floor (ru[0]);
      ru[1] -= std::floor (ru[1]);
      ru[2] -= std::floor (ru[2]);
      ValueType val;
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
      complex<double> eye (0.0, 1.0);
      for(int j=0; j<OrbitalSetSize; j++) {
// 	Orbitals[j].evaluate(ru, val, gu, hess);
// 	vals(j,i)  = val;
// 	grads(i,j) = dot(PrimLattice.G, gu);
// 	lapls(i,j) = trace(hess, GGt);

	Orbitals[j].evaluate(ru, val, gu, hess);
	complex<double> u(val);
	TinyVector<complex<double>,3> gradu;
	complex<double> laplu;
	u     = val;
	gradu = dot(PrimLattice.G, gu);
	laplu = trace(hess, GGt);
	
	PosType k = Orbitals[j].kVec;
	double s,c;
	double phase = -dot(P.R[iat], k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
	vals(j,i)  = e_mikr * u;
	grads(i,j) = e_mikr*(-eye * k * u + gradu);
	lapls(i,j) = e_mikr*(dot(k,k)*u - 2.0*eye*dot(k,gradu) + laplu);
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
