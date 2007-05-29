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
/** @file LRTwoBodyJastrow.h
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
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {

  class LRTwoBodyJastrow: public OrbitalBase {

    bool NeedToRestore;
    IndexType NumPtcls;
    IndexType NumSpecies;
    IndexType NumKpts;
    ///maximum K-point
    IndexType MaxK;
    ///number of k-shells
    IndexType MaxKshell;
    /// 1/Cell Volume
    RealType OneOverCellVolume;
    ///Omega 
    RealType Omega;
    ///4*pi*Omega 
    RealType FourPiOmega;
    ///1.0/Omega 
    RealType OneOverOmega;
    ///Rs 
    RealType Rs;
    ///Normalization Constant
    RealType NormConstant;
    RealType curVal, curLap;
    PosType curGrad;

    RealType *FirstAddressOfdU; 
    RealType *LastAddressOfdU;
    StructFact* skRef;
    // handler used to do evalFk
    typedef LRJastrowSingleton::LRHandlerType HandlerType;    
    HandlerType* handler;

    Vector<RealType> U,d2U;
    Vector<PosType> dU;
    Vector<RealType> offU, offd2U;
    Vector<PosType> offdU;

    Matrix<ComplexType> rokbyF;
    Vector<ComplexType> Rhok;

    ////Matrix<ComplexType> rhok;
    Matrix<ComplexType> eikr;
    Vector<ComplexType> eikr_new;
    Vector<ComplexType> delta_eikr;

    vector<int> Kshell;
    vector<RealType> FkbyKK;
    //vector<PosType>  Kcart;
    
  public:
    Vector<RealType> Fk_0; 
    Vector<RealType> Fk_1; 
    ///Coefficients = Fk_0 + Fk_1
    Vector<RealType> Fk; 
    ///A unique Fk sorted by |k|
    Vector<RealType> Fk_symm;

    LRTwoBodyJastrow(ParticleSet& p, HandlerType* inhandler);

    void resetParameters(OptimizableSetType& optVariables);

    void resetInternals();

    void resize();

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
    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    ///process input file
    bool put(xmlNodePtr cur, VarRegistry<RealType>& vlist);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
