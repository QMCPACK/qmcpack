// -*- C++ -*-
#ifndef OHMMS_QMC_ORBITALBASE_H
#define OHMMS_QMC_ORBITALBASE_H
#include "Configuration.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/RecordProperty.h"

/**@file OrbitalBase.h
 *@brief Declaration of OrbitalBase
 */
namespace ohmmsqmc {

  class WalkerSetRef;
  class ParticleSet;

  /** An abstract class for many-body orbitals
   *
   * TrialWaveFunction is a product of  OrbitalBase objects.
   */
  struct OrbitalBase: public QMCTraits {

    ///recasting enum of DistanceTableData to maintain consistency
    enum {SourceIndex  = DistanceTableData::SourceIndex, 
	  VisitorIndex = DistanceTableData::VisitorIndex, 
	  WalkerIndex  = DistanceTableData::WalkerIndex};

    typedef ParticleAttrib<ValueType>     ValueVectorType;

    OrbitalBase(){ }
    virtual ~OrbitalBase() { }

    virtual void reset() = 0;

    ///resize the internal storage with number of walkers if necessary
    virtual void resizeByWalkers(int nwalkers) {}

    virtual ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L) = 0;

    virtual void 
    evaluate(WalkerSetRef& W, 
	     ValueVectorType& psi,
	     WalkerSetRef::WalkerGradient_t& G,
	     WalkerSetRef::WalkerLaplacian_t& L) = 0;
  };


}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

