//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"
#include "LongRange/StructFact.h"
#include "Numerics/e2iphi.h"
#include <sstream>
#include <algorithm>

namespace qmcplusplus {

  void
  kSpaceJastrow::StructureFactor(PosType G, std::vector<ComplexType> &rho_G)
  {
    for (int i=0; i<NumIonSpecies; i++)
      rho_G[i] = ComplexType();
    for (int iat=0; iat<Ions.getLocalNum(); iat++) {
      PosType r(Ions.R[iat]);
      RealType phase = dot(r,G);
      int id = Ions.GroupID[iat];
      rho_G[id] += ComplexType(std::cos(phase), std::sin(phase));
    }
  }

  inline bool 
  Include (int i, int j, int k) {
    if (i > 0)
      return true;
    else if (i==0) {
      if (j>0)
	return true;
      else if ((j==0) && (k>0))
	return true;
    }
    return false;
  }

  void
  kSpaceJastrow::setupGvecs(RealType kc, std::vector<PosType> &gvecs)
  {
    gvecs.clear();
    int maxIndex[OHMMS_DIM];
    for (int i=0; i<OHMMS_DIM; i++)
      maxIndex[i] = 2 + (int)std::floor
	(std::sqrt(dot(Ions.Lattice.a(i), Ions.Lattice.a(i)))*kc / (2.0*M_PI));
    std::vector<ComplexType> rho_G(NumIonSpecies);
    for (int i=0; i<=maxIndex[0]; i++)
      for (int j=-maxIndex[1]; j<=maxIndex[1]; j++)
	for (int k=-maxIndex[2]; k<=maxIndex[2]; k++) {
	  // Omit half the G-vectors because of time-reversal symmetry
	  if (Include(i,j,k)) {
	    PosType G = 2.0*M_PI*((RealType)i*Ions.Lattice.Gv[0] + 
				  (RealType)j*Ions.Lattice.Gv[1] +
				  (RealType)k*Ions.Lattice.Gv[2]);
	    if (dot(G,G) <= (kc*kc)) {
	      bool notZero(false);
	      StructureFactor(G, rho_G);	      
	      for (int sp=0; sp<NumIonSpecies; sp++) 
		notZero = notZero || (norm(rho_G[sp]) > 1.0e-12);
	      if (notZero)
		gvecs.push_back(G);
	    }
	  }
	}
  }

  struct magLess
  {
    inline bool operator()(kSpaceJastrow::PosType a, kSpaceJastrow::PosType b)
    {
      return dot(a,a) < dot(b,b);
    }
  };

  bool
  kSpaceJastrow::operator()(PosType G1, PosType G2)
  {
    if (std::fabs(dot(G1,G1) - dot(G2,G2)) > 1.0e-8)
      return dot(G1,G1) < dot(G2,G2);
    // if (Equivalent(G1,G2)) return false;
    vector<ComplexType> rho_G1(NumIonSpecies), rho_G2(NumIonSpecies);
    StructureFactor(G1, rho_G1);
    StructureFactor(G2, rho_G2);
    for (int i=0; i<NumIonSpecies; i++ ) 
      for (int j=i+1; j<NumIonSpecies; j++) { 
	ComplexType zG1 = rho_G1[i]*conj(rho_G1[j]);
	ComplexType zG2 = rho_G2[i]*conj(rho_G2[j]);
	double SG1  = std::real(zG1);
	double SG2  = std::real(zG2);
	if (std::fabs(SG1 - SG2) > 1.0e-8)
	  return SG1 < SG2;
      }
    return false;
  }

  bool
  kSpaceJastrow::Equivalent(PosType G1, PosType G2)
  {
    return (!(*this)(G1,G2) && !(*this)(G2,G1));

    if (std::fabs(dot(G1,G1) - dot(G2,G2)) > 1.0e-8)
      return false;
    vector<ComplexType> rho_G1(NumIonSpecies), rho_G2(NumIonSpecies);
    StructureFactor(G1, rho_G1);
    StructureFactor(G2, rho_G2);
    
    for (int j=0; j<NumIonSpecies; j++) 
      for (int i=j; i<NumIonSpecies; i++ ) {
	ComplexType zG1 = rho_G1[i]*conj(rho_G1[j]);
	ComplexType zG2 = rho_G2[i]*conj(rho_G2[j]);
	double SG1  = std::real(zG1);
	double SG2  = std::real(zG2);
	if (std::fabs(SG1 - SG2) > 1.0e-8)
	  return false;
    }
    return true;
  }

  template<typename T> void
  kSpaceJastrow::sortGvecs(std::vector<PosType> &gvecs, 
			   std::vector<kSpaceCoef<T> > &coefs,
			   SymmetryType symm)
  {
    if (!gvecs.size())
      return;
    if (symm == CRYSTAL) {
      // First pass:  sort all the G-vectors into equivalent groups
      std::sort(gvecs.begin(), gvecs.end(), *this);
      
      // Now, look through the sorted G-vectors and group them
      kSpaceCoef<T> coef;
      coef.cG = T();
      coef.firstIndex = coef.lastIndex = 0;
      for (int i=1; i<gvecs.size(); i++)
	if (Equivalent(gvecs[i], gvecs[i-1]))
	  coef.lastIndex=i;
	else {
	  coefs.push_back(coef);
	  coef.firstIndex = coef.lastIndex = i;
	}
      coefs.push_back(coef);
    }
    else if (symm == ISOTROPIC) {
      magLess comparator;
      std::sort(gvecs.begin(), gvecs.end(), comparator);
      double curMag2 = dot(gvecs[0], gvecs[0]);
      kSpaceCoef<T> coef;
      coef.cG = T();
      coef.firstIndex = coef.lastIndex = 0;
      for (int i=1; i<gvecs.size(); i++) {
	double mag2 = dot(gvecs[i],gvecs[i]);
	if (std::fabs(mag2-curMag2) < 1.0e-10) 
	  coef.lastIndex = i;
	else {
	  coefs.push_back(coef);
	  coef.firstIndex = coef.lastIndex = i;
	  curMag2 = mag2;
	}
      }
      coefs.push_back(coef);
    }
    else if (symm == NOSYMM) {
      coefs.resize(gvecs.size());
      for (int i=0; i<gvecs.size(); i++) {
	coefs[i].cG = T();
	coefs[i].firstIndex = coefs[i].lastIndex = i;
      }
    }
    app_log() << "Using a total of " << gvecs.size() << " G-vectors in " << coefs.size()
	      << " symmetry groups.\n";
    app_log() << "kSpace coefficent groups:\n";
    for (int i=0; i<coefs.size(); i++) {
      app_log() << "  Group " << i << ":\n";
      for (int j=coefs[i].firstIndex; j<=coefs[i].lastIndex; j++) 
	app_log() << "    " << gvecs[j] << "    " << std::sqrt(dot(gvecs[j], gvecs[j])) << endl;
    }
  }

  kSpaceJastrow::kSpaceJastrow(ParticleSet& ions, ParticleSet& elecs,
			       SymmetryType oneBodySymm, RealType oneBodyCutoff, string oneBodyID,
			       SymmetryType twoBodySymm, RealType twoBodyCutoff, string twoBodyID)
    : Ions(ions), Elecs(elecs)
  {
    NumIonSpecies = 0;
    NumElecs = elecs.getLocalNum();
    for (int iat=0; iat<ions.getLocalNum(); iat++) 
      NumIonSpecies = max(NumIonSpecies, ions.GroupID[iat]+1);

    if (oneBodyCutoff > 0.0) {
      setupGvecs(oneBodyCutoff, OneBodyGvecs);
      sortGvecs (OneBodyGvecs, OneBodySymmCoefs, oneBodySymm);
      for (int i=0; i<OneBodySymmCoefs.size(); i++) {
	stringstream name_real, name_imag;
	name_real << oneBodyID << "_" << 2*i;
	name_imag << oneBodyID << "_" << 2*i+1;
	VarMap[name_real.str()] = &(OneBodySymmCoefs[i].cG.real());
	VarMap[name_imag.str()] = &(OneBodySymmCoefs[i].cG.imag());
      }
    }
    if (twoBodyCutoff > 0.0) {
      setupGvecs(twoBodyCutoff, TwoBodyGvecs);
      sortGvecs (TwoBodyGvecs, TwoBodySymmCoefs, twoBodySymm);
      for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
	stringstream name;
	name << twoBodyID << "_" << i;
	VarMap[name.str()] = &(TwoBodySymmCoefs[i].cG);
      }

    }
    // Now resize all buffers
    int nOne = OneBodyGvecs.size();
    OneBodyCoefs.resize(nOne);
    OneBody_rhoG.resize(nOne);
    Ion_rhoG.resize(nOne);
    OneBodyPhase.resize(nOne);
    OneBody_e2iGr.resize(nOne);

    int nTwo = TwoBodyGvecs.size();
    int nElecs = Elecs.getTotalNum();
    TwoBodyCoefs.resize(nTwo);
    TwoBody_rhoG.resize(nOne);
    TwoBodyPhase.resize(nTwo);
    TwoBody_e2iGr_new.resize(nTwo);
    TwoBody_e2iGr_old.resize(nTwo);
    Delta_e2iGr.resize(nElecs,nTwo);
    for (int iat=0; iat<nElecs; iat++)
      for (int i=0; i<nTwo; i++)
	Delta_e2iGr(iat,i) = ComplexType();

    // Set Ion_rhoG
    for (int i=0; i<OneBodyGvecs.size(); i++) {
      Ion_rhoG[0] = ComplexType();
      for (int iat=0; iat<ions.getLocalNum(); iat++) {
	double phase = dot(OneBodyGvecs[i],ions.R[iat]);
	Ion_rhoG[i] += ComplexType(std::cos(phase), std::sin(phase));
      }
    }
	
  }
  
  void 
  kSpaceJastrow::setCoefficients(std::vector<RealType> &oneBodyCoefs,
				 std::vector<RealType> &twoBodyCoefs)
  {
    for (int i=0; i<oneBodyCoefs.size(); i++) 
      cerr << "oneBodyCoefs[" << i << "] = " << oneBodyCoefs[i] << endl;
    if (oneBodyCoefs.size() != 2*OneBodySymmCoefs.size()) {
      app_warning() << "Warning!  Wrong number of coefficients specified in "
		    << "kSpaceJastrow's one-body coefficients.\n"
		    << oneBodyCoefs.size() << " were specified.  Should have been "
		    << 2*OneBodySymmCoefs.size() << endl;
      for (int i=0; i<OneBodySymmCoefs.size(); i++) {
	OneBodySymmCoefs[i].cG = ComplexType();
	OneBodySymmCoefs[i].set (OneBodyCoefs);
      }
    }
    else {
      for (int i=0; i<OneBodySymmCoefs.size(); i++) {
	OneBodySymmCoefs[i].cG = ComplexType (oneBodyCoefs[2*i+0],
					      oneBodyCoefs[2*i+1]);
	OneBodySymmCoefs[i].set (OneBodyCoefs);
      }
    }

    if (twoBodyCoefs.size() != TwoBodySymmCoefs.size()) {
      app_warning() << "Warning!  Wrong number of coefficients specified in "
		    << "kSpaceJastrow's two-body coefficients.\n"
		    << twoBodyCoefs.size() << " were specified.  Should have been "
		    << TwoBodySymmCoefs.size() << endl;
      for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
	TwoBodySymmCoefs[i].cG = 0.0;
	TwoBodySymmCoefs[i].set (TwoBodyCoefs);
      }
    }
    else {
      for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
	TwoBodySymmCoefs[i].cG = twoBodyCoefs[i];
	TwoBodySymmCoefs[i].set (TwoBodyCoefs);
      }
    }

  }

  void 
  kSpaceJastrow::resetTargetParticleSet(ParticleSet& P) 
  {
    for (int i=0; i<TwoBodyGvecs.size(); i++) {
      TwoBody_rhoG[i] = ComplexType();
      for (int iat=0; iat<NumElecs; iat++) {
	double phase, s, c;
	phase = dot (TwoBodyGvecs[i], P.R[iat]);
	sincos(phase, &s, &c);
	TwoBody_rhoG[i] += ComplexType(c,s);
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  //                  Evaluation functions                     //
  ///////////////////////////////////////////////////////////////

  kSpaceJastrow::ValueType 
  kSpaceJastrow::evaluateLog(ParticleSet& P, 
			     ParticleSet::ParticleGradient_t& G, 
			     ParticleSet::ParticleLaplacian_t& L) 
  {
    RealType J1(0.0), J2(0.0);
    int N = P.getTotalNum();
    ComplexType eye(0.0, 1.0);
    int nOne = OneBodyGvecs.size();
    
    for (int iat=0; iat<N; iat++) {
      PosType r(P.R[iat]);
      for (int i=0; i<nOne; i++) 
	OneBodyPhase[i] = dot(OneBodyGvecs[i], r);
      eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
      
      for (int i=0; i<nOne; i++) {
	ComplexType z = OneBodyCoefs[i] * conj(OneBody_e2iGr[i]);
	J1 += real(z);
	G[iat] += -real(z*eye)*OneBodyGvecs[i];
	L[iat] += -dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z);
      }
    }

    // Do two-body part
    int nTwo = TwoBodyGvecs.size();
    for (int i=0; i<nTwo; i++)
      TwoBody_rhoG[i] = ComplexType();
    for (int iat=0; iat<N; iat++) {
      PosType r(P.R[iat]);
      for (int i=0; i<nTwo; i++) 
	TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
      eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
      for (int i=0; i<nTwo; i++)
	TwoBody_rhoG[i] += TwoBody_e2iGr_new[i];
    }

    for (int i=0; i<nTwo; i++) 
      J2 += TwoBodyCoefs[i]*norm(TwoBody_rhoG[i]);
    
    for (int iat=0; iat<N; iat++) {
      PosType r(P.R[iat]);
      for (int i=0; i<nTwo; i++) 
	TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
      eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
      for (int i=0; i<nTwo; i++) {
	PosType Gvec(TwoBodyGvecs[i]);
	ComplexType z = TwoBody_e2iGr_new[i];
	G[iat] += -2.0*Gvec*TwoBodyCoefs[i]*imag(conj(TwoBody_rhoG[i])*z);
	L[iat] += 2.0*TwoBodyCoefs[i]*dot(Gvec,Gvec)*(-real(z*conj(TwoBody_rhoG[i])) + 1.0);
      }
    }
    return J1 + J2;
  }
  
  
  /* evaluate the ratio with P.R[iat]
   *
   */
  kSpaceJastrow::ValueType 
  kSpaceJastrow::ratio(ParticleSet& P, int iat) 
  {
    RealType J1new(0.0), J1old(0.0), J2new(0.0), J2old(0.0);
    PosType rnew(P.R[iat]), rold(P.getOldPos());
    // Compute one-body contribution
    int nOne = OneBodyGvecs.size();
    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);

    for (int i=0; i<nOne; i++) 
      J1new += real(OneBodyCoefs[i] * conj(OneBody_e2iGr[i]));

    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);

    for (int i=0; i<nOne; i++) 
      J1old += real(OneBodyCoefs[i] * conj(OneBody_e2iGr[i]));

    // Now, do two-body part
    int nTwo = TwoBodyGvecs.size();
    for (int i=0; i<nTwo; i++)   TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_old);
    for (int i=0; i<nTwo; i++)   TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);

    for (int i=0; i<nTwo; i++) {
      ComplexType rho_G = TwoBody_rhoG[i];
      J2old += TwoBodyCoefs[i]*std::norm(rho_G);
    }
    for (int i=0; i<nTwo; i++) {
      ComplexType rho_G = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
      J2new += TwoBodyCoefs[i]*std::norm(rho_G);
    }

    return std::exp(J1new+J2new - (J1old + J2old));
  }
  
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::logRatio(ParticleSet& P, int iat,
			  ParticleSet::ParticleGradient_t& dG,
			  ParticleSet::ParticleLaplacian_t& dL) 
  {
    RealType J1(0.0), J2(0.0);
    PosType rnew(P.R[iat]), rold(P.getOldPos());
    ComplexType eye(0.0, 1.0);

    // Compute one-body contribution
    int nOne = OneBodyGvecs.size();
    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);

    for (int i=0; i<nOne; i++) {
      ComplexType z = OneBodyCoefs[i] * conj(OneBody_e2iGr[i]);
      J1 += real(z);
      dG[iat] += -real(z*eye) * OneBodyGvecs[i];
      dL[iat] += -dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z);
    }

    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);

    for (int i=0; i<nOne; i++) {
      ComplexType z = OneBodyCoefs[i] * conj(OneBody_e2iGr[i]);
      J1 -= real(z);
      dG[iat] -= -real(z*eye) * OneBodyGvecs[i];
      dL[iat] -= -dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z);
    }

    // Compute two-body contribution
    int nTwo = TwoBodyGvecs.size();
    for (int i=0; i<nTwo; i++)   TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
    eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_new);
    for (int i=0; i<nTwo; i++)   TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
    eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_old);

    for (int jat=0; jat<NumElecs; jat++)
      for (int i=0; i<nTwo; i++)
	Delta_e2iGr(jat,i) = ComplexType();

    for (int i=0; i<nTwo; i++)
      Delta_e2iGr(iat,i) = TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];

    for (int i=0; i<nTwo; i++) {
      ComplexType rhoG_old = TwoBody_rhoG[i];
      ComplexType rhoG_new = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
      J2 += TwoBodyCoefs[i] * (norm(rhoG_new) - norm(rhoG_old));
      for (int jat=0; jat<NumElecs; jat++) {
	PosType rj_new = P.R[jat];
	PosType rj_old = rj_new;
	if (iat == jat) 
	  rj_old = P.getOldPos();
	PosType Gvec(TwoBodyGvecs[i]);
	ComplexType zold, znew;
	sincos (dot(Gvec,rj_old), &(zold.imag()), &(zold.real()));
	sincos (dot(Gvec,rj_new), &(znew.imag()), &(znew.real()));
	dG[jat] += -2.0*Gvec*TwoBodyCoefs[i]*imag(conj(rhoG_new)*znew);
	dG[jat] -= -2.0*Gvec*TwoBodyCoefs[i]*imag(conj(rhoG_old)*zold);

	dL[jat] += 2.0*TwoBodyCoefs[i]*dot(Gvec,Gvec)*(-real(znew*conj(rhoG_new)) + 1.0);
	dL[jat] -= 2.0*TwoBodyCoefs[i]*dot(Gvec,Gvec)*(-real(zold*conj(rhoG_old)) + 1.0);
      }

    }
    return J1 + J2;
  }
  
  void 
  kSpaceJastrow::restore(int iat) 
  {
    //substract the addition in logRatio
    //if(NeedToRestore) Rhok -= delta_eikr;
  }

  void 
  kSpaceJastrow::acceptMove(ParticleSet& P, int iat) 
  {
    for (int i=0; i<TwoBody_e2iGr_new.size(); i++) 
      TwoBody_rhoG[i] += Delta_e2iGr(iat,i);
    //TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];

    // std::copy(eikr_new.data(),eikr_new.data()+MaxK,eikr[iat]);
    // U += offU;
    // dU += offdU;
    // d2U += offd2U;
  }

  void 
  kSpaceJastrow::update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) 
  {
    app_error() << "kSpaceJastrow::update is INCOMPLETE " << endl;
  }
  
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::registerData(ParticleSet& P, PooledData<RealType>& buf) 
  {
    LogValue=evaluateLog(P,P.G,P.L); 
    // eikr.resize(NumPtcls,MaxK);
    // eikr_new.resize(MaxK);
    // delta_eikr.resize(MaxK);
    
    // for(int iat=0; iat<NumPtcls; iat++)
    //   std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
    
    // buf.add(Rhok.first_address(), Rhok.last_address());
    // buf.add(U.first_address(), U.last_address());
    // buf.add(d2U.first_address(), d2U.last_address());
    // buf.add(FirstAddressOfdU,LastAddressOfdU);
    // return LogValue;
    return LogValue;
  }
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
			      bool fromscratch) 
  {
    LogValue=evaluateLog(P,P.G,P.L); 
    
    // for(int iat=0; iat<NumPtcls; iat++)
    //   std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
    
    // buf.put(Rhok.first_address(), Rhok.last_address());
    // buf.put(U.first_address(), U.last_address());
    // buf.put(d2U.first_address(), d2U.last_address());
    // buf.put(FirstAddressOfdU,LastAddressOfdU);
    // return LogValue;
    return LogValue;
  }

  void 
  kSpaceJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for (int i=0; i<TwoBodyCoefs.size(); i++) 
      TwoBody_rhoG[i] = ComplexType();
    for (int iat=0; iat<NumElecs; iat++) {
      for (int i=0; i<TwoBodyCoefs.size(); i++) 
	TwoBodyPhase[i] = dot(TwoBodyGvecs[i],P.R[iat]);
      eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_new);
      for (int i=0; i<TwoBodyCoefs.size(); i++)
	TwoBody_rhoG[i] += TwoBody_e2iGr_new[i];
    }
  }
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::evaluate(ParticleSet& P, PooledData<RealType>& buf) 
  {
    RealType J1(0.0), J2(0.0);
    int N = P.getTotalNum();
    for (int iat=0; iat<N; iat++) {
      PosType r(P.R[iat]);
      int nOne = OneBodyGvecs.size();
      for (int i=0; i<nOne; i++) 
	OneBodyPhase[i] = dot(OneBodyGvecs[i], r);
      eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
      
      for (int i=0; i<nOne; i++) {
	J1 +=  real(OneBodyCoefs[i]  * conj(OneBody_e2iGr[i]));
      }
    }
    
    // Do two-body part
    int nTwo = OneBodyGvecs.size();
    for (int i=0; i<nTwo; i++)
      TwoBody_rhoG[i] = ComplexType();
    for (int iat=0; iat<N; iat++) {
      PosType r(P.R[iat]);
      for (int i=0; i<nTwo; i++) 
	TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
      eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
      for (int i=0; i<nTwo; i++)
	TwoBody_rhoG[i] += TwoBody_e2iGr_new[i];
    }
    for (int i=0; i<nTwo; i++) 
      J2 += TwoBodyCoefs[i]*norm(TwoBody_rhoG[i]);
    return std::exp(J1 + J2);
  }
  
  void
  kSpaceJastrow::addOptimizables(OptimizableSetType& vlist)
  {
    std::map<std::string,RealType*>::iterator iter;
    for (iter=VarMap.begin(); iter!=VarMap.end(); iter++) {
      string name = iter->first;
      RealType val = *(iter->second);
      vlist[name] = val;
    }
  }

  void 
  kSpaceJastrow::resetParameters(OptimizableSetType& optVariables) 
  { 
    OptimizableSetType::iterator var;
    for (var=optVariables.begin(); var!=optVariables.end(); var++) {
      std::map<std::string,RealType*>::iterator myVar
	= VarMap.find(var->first);
      if (myVar != VarMap.end())  {
	*(myVar->second) = var->second;
      }
    }
    for (int i=0; i<OneBodySymmCoefs.size(); i++)
      OneBodySymmCoefs[i].set(OneBodyCoefs);
    for (int i=0; i<TwoBodySymmCoefs.size(); i++)
      TwoBodySymmCoefs[i].set(TwoBodyCoefs);
  }

  bool
  kSpaceJastrow::put(xmlNodePtr cur, VarRegistry<RealType>& vlist) 
  {
    return true;
  }
}

