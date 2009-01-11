//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LOCALENERGYESTIMATOR_HDF_H
#define QMCPLUSPLUS_LOCALENERGYESTIMATOR_HDF_H
#include <Estimators/ScalarEstimatorBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>
#include <QMCHamiltonians/observable_helper.h>

namespace qmcplusplus {

  /** Handle an ensemble average of Hamiltonian components
  */
  class LocalEnergyEstimatorHDF: public ScalarEstimatorBase {

    //typedef PooledData<T>                            BufferType;
    //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};
    enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

    //int LocalPotentialIndex;
    int FirstHamiltonian;
    int SizeOfHamiltonians;

    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;
    ///save the reference hamiltonian
    const QMCHamiltonian& refH;

  public:

    /** constructor
     * @param h QMCHamiltonian to define the components
     */
    LocalEnergyEstimatorHDF(QMCHamiltonian& h);
    //LocalEnergyEstimatorHDF(const LocalEnergyEstimatorHDF& est);

    /** implement virtual function
     */
    ScalarEstimatorBase* clone();

    void add2Record(RecordListType& record);

    void accumulate(const Walker_t& awalker, RealType wgt);

    inline void accumulate(WalkerIterator first, WalkerIterator last, RealType wgt) {
      while(first != last) {
        accumulate(**first,wgt);
        ++first;
      }
    }
    
    void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: mcminis2 $
 * $Revision: 3421 $   $Date: 2008-12-09 10:21:11 -0600 (Tue, 09 Dec 2008) $
 * $Id: LocalEnergyEstimatorHDF.h 3421 2008-12-09 16:21:11Z mcminis2 $ 
 ***************************************************************************/
