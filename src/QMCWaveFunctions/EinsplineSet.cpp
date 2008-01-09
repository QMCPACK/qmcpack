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

  EinsplineSet::UnitCellType
  EinsplineSet::GetLattice()
  {
    return SuperLattice;
  }
  
  void
  EinsplineSet::resetParameters(VarRegistry<RealType>& vlist)
  {
  }
  
  void
  EinsplineSet::resetTargetParticleSet(ParticleSet& e)
  {
  }
  
  void
  EinsplineSet::setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
  }
  
  void 
  EinsplineSet::evaluate (const ParticleSet& P, int iat, 
			  ValueVector_t& psi)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    for(int j=0; j<OrbitalSetSize; j++) {
      complex<double> val;
      Orbitals[j]->evaluate(ru, val);

      double phase = -dot(r, Orbitals[j]->kVec);
      double s,c;
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      val *= e_mikr;
#ifdef QMC_COMPLEX
      psi[j] = val;
#else
      psi[j] = real(val);
#endif
    }
  }
  
  void 
  EinsplineSet::evaluate (const ParticleSet& P, int iat, 
			  ValueVector_t& psi, GradVector_t& dpsi, 
			  ValueVector_t& d2psi)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    complex<double> val;
    TinyVector<complex<double>,3> gu;
    Tensor<complex<double>,3> hess;
    complex<double> eye (0.0, 1.0);
    for(int j=0; j<OrbitalSetSize; j++) {
      complex<double> u;
      TinyVector<complex<double>,3> gradu;
      complex<double> laplu;

      Orbitals[j]->evaluate(ru, val, gu, hess);
      u  = val;
      // Compute gradient in cartesian coordinates
      gradu = dot(PrimLattice.G, gu);
      laplu = trace(hess, GGt);      
      
      PosType k = Orbitals[j]->kVec;
      TinyVector<complex<double>,3> ck;
      ck[0]=k[0];  ck[1]=k[1];  ck[2]=k[2];
      double s,c;
      double phase = -dot(P.R[iat], k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
#ifdef QMC_COMPLEX
      psi[j]   = e_mikr * u;
      dpsi[j]  = e_mikr*(-eye * ck * u + gradu);
      d2psi[j] = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
#else
      psi[j]   = real(e_mikr * u);
      dpsi[j]  = real(e_mikr*(-eye * ck * u + gradu));
      d2psi[j] = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
#endif

    }
  }
  
  void 
  EinsplineSet::evaluate (const ParticleSet& P, int first, int last,
			  ValueMatrix_t& vals, GradMatrix_t& grads, 
			  ValueMatrix_t& lapls)
  {
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType r (P.R[iat]);
      PosType ru(PrimLattice.toUnit(r));
      ru[0] -= std::floor (ru[0]);
      ru[1] -= std::floor (ru[1]);
      ru[2] -= std::floor (ru[2]);
      complex<double> val;
      TinyVector<complex<double>,3> gu;
      Tensor<complex<double>,3> hess;
      complex<double> eye (0.0, 1.0);
      for(int j=0; j<OrbitalSetSize; j++) {
	complex<double> u;
	TinyVector<complex<double>,3> gradu;
	complex<double> laplu;

	Orbitals[j]->evaluate(ru, val, gu, hess);
	u  = val;
	gradu = dot(PrimLattice.G, gu);
	laplu = trace(hess, GGt);
	
	PosType k = Orbitals[j]->kVec;
	TinyVector<complex<double>,3> ck;
	ck[0]=k[0];  ck[1]=k[1];  ck[2]=k[2];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
#ifdef QMC_COMPLEX
	vals(j,i)  = e_mikr * u;
	grads(i,j) = e_mikr*(-eye*u*ck + gradu);
	lapls(i,j) = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
#else
	vals(j,i)  = real(e_mikr * u);
	grads(i,j) = real(e_mikr*(-eye*u*ck + gradu));
	lapls(i,j) = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
#endif

      }
    }
  }
}
