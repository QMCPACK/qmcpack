//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file kSpaceJastrow.h
 * @brief Declaration of Long-range TwoBody Jastrow
 *
 * The initial coefficients are evaluated according to RPA.
 * A set of k-dependent coefficients are generated.
 * Optimization should use rpa_k0, ...rpa_kn as the names of the
 * optimizable parameters.
 */
#ifndef QMCPLUSPLUS_LR_RPAJASTROW_H
#define QMCPLUSPLUS_LR_RPAJASTROW_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "Optimize/VarList.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus {
  template<typename T>
  struct kSpaceCoef
  {
    // The coefficient
    T cG;
    // The first and last indices of the G-vectors to which this
    // coefficient belongs
    int firstIndex, lastIndex;
    inline void set (std::vector<T> &unpackedCoefs) {
      for (int i=firstIndex; i<=lastIndex; i++)
	unpackeCoefs[i] = cG;
    }
  }

  class kSpaceJastrow: public OrbitalBase {
  public:   
    typedef enum { CRYSTAL, ISOTROPIC, NOSYMM } SymmetryType;
  private:
    RealType CellVolume, NormConstant;
    int NumElecs, NumSpins;
    int NumIons, NumIonSpecies;

    // Gvectors included in the summation for the 
    // one-body and two-body Jastrows, respectively
    std::vector<PosType> OneBodyGvecs;
    std::vector<PosType> TwoBodyGvecs;
    
    // Vectors of unique coefficients, one for each group of
    // symmetry-equivalent G-vectors
    std::vector<kSpaceCoefs<ComplexType> > OneBodySymmCoefs;
    std::vector<kSpaceCoefs<RealType   > > TwoBodySymmCoefs;

    // Vector of coefficients, one for each included G-vector
    // in OneBodyGvecs/TwoBodyGvecs
    std::vector<ComplexType> OneBodyCoefs;
    std::vector<RealType>    TwoBodyCoefs;
    
    SymmetryType OneBodySymmType, TwoBodySymmType;

    // Enumerate G-vectors with nonzero structure factors
    void setupGvecs (RealType kcut, std::vector<PosType> &gvecs);
    void setupCoefs();
    
    ComplexType StructureFactor(int speciesNum, PosType G);
    
    ParticleSet &Ions;

  public:
    kSpaceJastrow(ParticleSet& ions, ParticleSet& elecs,
		  SymmetryType oneBodySymm, RealType oneBodyCutoff,
		  SymmetryType twoBodySymm, RealType twoBodyCutoff);

    void resetParameters(OptimizableSetType& optVariables);

    //evaluate the distance table with els
    void resetTargetParticleSet(ParticleSet& P);

    ValueType evaluateLog(ParticleSet& P,
		         ParticleSet::ParticleGradient_t& G, 
		         ParticleSet::ParticleLaplacian_t& L);

    inline ValueType evaluate(ParticleSet& P,
			      ParticleSet::ParticleGradient_t& G, 
			      ParticleSet::ParticleLaplacian_t& L) {
      return std::exp(evaluateLog(P,G,L));
    }


    ValueType ratio(ParticleSet& P, int iat);

    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL)  {
      return std::exp(logRatio(P,iat,dG,dL));
    }

    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

    void restore(int iat);
    void acceptMove(ParticleSet& P, int iat);
    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat);


    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);
    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    ///process input file
    bool put(xmlNodePtr cur, VarRegistry<RealType>& vlist);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2595 $   $Date: 2008-04-17 07:52:58 -0500 (Thu, 17 Apr 2008) $
 * $Id: kSpaceJastrow.h 2595 2008-04-17 12:52:58Z jnkim $ 
 ***************************************************************************/
