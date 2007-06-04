//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_STRUCTUREFACTOR_ESTIMATOR_H
#define QMCPLUSPLUS_STRUCTUREFACTOR_ESTIMATOR_H
#include "Estimators/CompositeEstimators.h"
#include "Estimators/VectorEstimatorImpl.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus {

  class SkEstimator: public CompositeEstimatorBase 
  {
    typedef VectorEstimatorImpl<RealType> VectorEstimatorType;
    /** number of species */
    int NumSpecies;
    /** number of kpoints */
    int NumK;
    /** number of kshells */
    int MaxKshell;
    /** normalization factor */
    RealType OneOverN;
    /** kshell counters */
    vector<int> Kshell;
    /** instantaneous structure factor  */
    vector<RealType> Kmag;
    /** 1.0/degenracy for a ksell */
    vector<RealType> OneOverDnk;
    /** \f$rho_k = \sum_{\alpha} \rho_k^{\alpha} \f$ for species index \f$\alpha\f$ */
    Vector<ComplexType> RhokTot;
    /** instantaneous structure factor  */
    Vector<RealType> SkInst;
    /** Structrue factor estimator */
    VectorEstimatorType Sk;
    /** hdf5 handler for Sk */
    HDFAttribIO<VectorEstimatorType>* Sk_h;
    public:

    /** constructor
     * @param source particleset
     */
    SkEstimator(ParticleSet& source);

    /** virtal destrctor */
    ~SkEstimator();

    //@{
    ///implement virtual functions
    void resetTargetParticleSet(ParticleSet& p);
    void open(hid_t hroot);
    void close();
    void startAccumulate();
    void accumulate(ParticleSet& p);
    void stopAccumulate(RealType wgtinv);
    void startBlock(int steps);
    void stopBlock(RealType wgtnorm, RealType errnorm);
    //@}

    private:
    SkEstimator(const SkEstimator& pc) {}
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/
