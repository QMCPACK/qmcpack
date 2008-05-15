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
      //maxIndex[i] = (int)std::ceil(2.0*cut / 2.0*M_PI*std::sqrt(Ions.Lattice.Gv[0],Ions.Lattice.Gv[i]));
      maxIndex[i] = 1 + (int)std::floor
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
			       SymmetryType oneBodySymm, RealType oneBodyCutoff,
			       SymmetryType twoBodySymm, RealType twoBodyCutoff)
    : Ions(ions), Elecs(elecs)
  {
    NumIonSpecies = 0;
    for (int iat=0; iat<ions.getLocalNum(); iat++) 
      NumIonSpecies = max(NumIonSpecies, ions.GroupID[iat]+1);

    if (oneBodyCutoff > 0.0) {
      setupGvecs(oneBodyCutoff, OneBodyGvecs);
      OneBodyCoefs.resize(OneBodyGvecs.size());
      sortGvecs (OneBodyGvecs, OneBodySymmCoefs, oneBodySymm);
      for (int i=0; i<OneBodySymmCoefs.size(); i++) {
	stringstream name_real, name_imag;
	name_real << "cG1_" << 2*i;
	name_imag << "cG1_" << 2*i+1;
	VarMap[name_real.str()] = &(OneBodySymmCoefs[i].cG.real());
	VarMap[name_imag.str()] = &(OneBodySymmCoefs[i].cG.real());
      }
    }
    if (twoBodyCutoff > 0.0) {
      setupGvecs(twoBodyCutoff, TwoBodyGvecs);
      TwoBodyCoefs.resize(TwoBodyGvecs.size());
      sortGvecs (TwoBodyGvecs, TwoBodySymmCoefs, twoBodySymm);
      for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
	stringstream name;
	name << "cG2_" << i;
	VarMap[name.str()] = &(TwoBodySymmCoefs[i].cG);
      }

    }
  }
      
  void kSpaceJastrow::resetTargetParticleSet(ParticleSet& P) 
  {
  }
  
  kSpaceJastrow::ValueType 
    kSpaceJastrow::evaluateLog(ParticleSet& P, 
			       ParticleSet::ParticleGradient_t& G, 
			       ParticleSet::ParticleLaplacian_t& L) 
  {
    return 0.0;
  }
  
  
  /* evaluate the ratio with P.R[iat]
   *
   */
  kSpaceJastrow::ValueType 
  kSpaceJastrow::ratio(ParticleSet& P, int iat) 
  {
    return 1.0;
  }
  
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::logRatio(ParticleSet& P, int iat,
			  ParticleSet::ParticleGradient_t& dG,
			  ParticleSet::ParticleLaplacian_t& dL) 
  {
    return 0.0;
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
    // LogValue=evaluateLog(P,P.G,P.L); 
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
    return 0.0;
  }
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
			      bool fromscratch) 
  {
    // LogValue=evaluateLog(P,P.G,P.L); 
    
    // for(int iat=0; iat<NumPtcls; iat++)
    //   std::copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
    
    // buf.put(Rhok.first_address(), Rhok.last_address());
    // buf.put(U.first_address(), U.last_address());
    // buf.put(d2U.first_address(), d2U.last_address());
    // buf.put(FirstAddressOfdU,LastAddressOfdU);
    // return LogValue;
    return 0.0;
  }

  void 
  kSpaceJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {    
  }
  
  kSpaceJastrow::ValueType 
  kSpaceJastrow::evaluate(ParticleSet& P, PooledData<RealType>& buf) 
  {
    return 0.0;
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
      if (myVar != VarMap.end()) 
	*(myVar->second) = var->second;
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

