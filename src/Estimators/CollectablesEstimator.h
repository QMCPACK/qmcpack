//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Jeongnim Kim 
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
#ifndef QMCPLUSPLUS_COLLECTABLES_ESTIMATOR_H
#define QMCPLUSPLUS_COLLECTABLES_ESTIMATOR_H
#include <Estimators/ScalarEstimatorBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>

namespace qmcplusplus {

  /** Handle an ensemble average of Hamiltonian components
  */
  class CollectablesEstimator: public ScalarEstimatorBase
  {
    ///save the reference hamiltonian
    const QMCHamiltonian& refH;
    int Nthreads;

    public:
    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    CollectablesEstimator(QMCHamiltonian& h);
    //LocalEnergyEstimatorHDF(const LocalEnergyEstimatorHDF& est);

    /** implement virtual function
    */
    CollectablesEstimator* clone();

    void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
    void add2Record(RecordListType& record);
    void accumulate(const MCWalkerConfiguration& W
        , WalkerIterator first, WalkerIterator last , RealType wgt);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: mcminis2 $
 * $Revision: 3421 $   $Date: 2008-12-09 10:21:11 -0600 (Tue, 09 Dec 2008) $
 * $Id: LocalEnergyEstimatorHDF.h 3421 2008-12-09 16:21:11Z mcminis2 $ 
 ***************************************************************************/
