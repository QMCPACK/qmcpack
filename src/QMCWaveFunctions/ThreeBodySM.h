//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and John Gergely
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
/**@file ThreeBodySM.h
 * @brief Declaration of three-body Schmidt-Moskowitz form using generic functor
 */
#ifndef QMCPLUSPLUS_GENERIC_THREEBODYSM_H
#define QMCPLUSPLUS_GENERIC_THREEBODYSM_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/NumericalJastrowFunctor.h"

namespace qmcplusplus {
  
  class ThreeBodySM: public OrbitalBase {

    typedef NumericalJastrow<RealType> FuncType;
    typedef vector<FuncType*> BasisType;

    ParticleSet& CenterRef;
    DistanceTableData* dist_ee;
    DistanceTableData* dist_ie;
    BasisType eeBasis;
    vector<BasisType*> ieBasis;

    ///constructor
    ThreeBodySM(ParticleSet& ions, ParticleSet& els);

    ~ThreeBodySM();

    //evaluate the distance table with P
    void resetTargetParticleSet(ParticleSet& P);

    void reset() { }

    ValueType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      return exp(evaluateLog(P,G,L));
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
