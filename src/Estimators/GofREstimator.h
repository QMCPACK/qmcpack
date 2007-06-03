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
#ifndef QMCPLUSPLUS_PAIRCORRELATION_ESTIMATOR_H
#define QMCPLUSPLUS_PAIRCORRELATION_ESTIMATOR_H
#include "Estimators/CompositeEstimators.h"
#include "Estimators/VectorEstimatorImpl.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus {

  class GofREstimator: public CompositeEstimatorBase 
  {
    typedef VectorEstimatorImpl<RealType> VectorEstimatorType;
    ///true if source == target
    bool Symmetric;
    /** number of centers */
    int Centers;
    /** number of distinct pair types */
    int NumPairTypes;
    /** number bins for gofr */
    int NumBins;
    /** maximum distance */
    RealType Dmax;
    /** bin size */
    RealType Delta;
    /** one of bin size */
    RealType DeltaInv;
    ///save the source particleset
    const ParticleSet& sourcePtcl;
    /** distance table */
    const DistanceTableData*  myTable;
    /** local copy of pair index */
    vector<int> PairID;
    /** normalization factor for each bin*/
    vector<RealType> normFactor;
    /** gofr[i] = gofr of i pair relations */
    vector<VectorEstimatorType*> gofr;
    /** instantaneous gofr */
    Matrix<RealType> gofrInst;
    /** name of the pair */
    vector<string> PairName;
    //vector<ostream*> fout;
    public:

    /** constructor
     * @param source particleset
     */
    GofREstimator(ParticleSet& source);

    /** constructor
     * @param source particleset
     * @param target particleset
     */
    GofREstimator(const ParticleSet& source, ParticleSet& target);

    /** virtal destrctor */
    ~GofREstimator();

    void resetTargetParticleSet(ParticleSet& p);
    void open(hid_t hroot);
    void close();
    /** prepare data collect */
    void startAccumulate();
    /** accumulate the observables */
    void accumulate(ParticleSet& p);
    /** reweight of the current cummulative  values */
    void stopAccumulate(RealType wgtinv);
    void startBlock(int steps);
    void stopBlock(RealType wgtnorm, RealType errnorm);
    void setBound(RealType dr);

    private:
    GofREstimator(const GofREstimator& pc): sourcePtcl(pc.sourcePtcl) {}
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/
