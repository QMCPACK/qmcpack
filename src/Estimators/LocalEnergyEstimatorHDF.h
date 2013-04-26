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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LOCALENERGYESTIMATOR_HDF_H
#define QMCPLUSPLUS_LOCALENERGYESTIMATOR_HDF_H
#include <Estimators/ScalarEstimatorBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>
#include <QMCHamiltonians/observable_helper.h>

namespace qmcplusplus
{

/** Handle an ensemble average of Hamiltonian components
*/
class LocalEnergyEstimatorHDF: public ScalarEstimatorBase
{

  //typedef PooledData<T>                            BufferType;
  //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};
  enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

  //int LocalPotentialIndex;
  int FirstHamiltonian;
  int SizeOfHamiltonians;

  ///save the reference hamiltonian
  const QMCHamiltonian& refH;

public:

  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  LocalEnergyEstimatorHDF(QMCHamiltonian& h);

  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
    const RealType* restrict ePtr = awalker.getPropertyBase();
    //weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
    RealType wwght= wgt* awalker.Weight;
    // RealType wwght= wgt;
    scalars[0](ePtr[LOCALENERGY],wwght);
    scalars[1](ePtr[LOCALPOTENTIAL],wwght);
    for(int target=2, source=FirstHamiltonian; target<scalars.size();
        ++target, ++source)
      scalars[target](ePtr[source],wwght);
  }

  /*@{*/
  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first!=last; ++first)
      accumulate(**first,wgt);
  }
  void add2Record(RecordListType& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: mcminis2 $
 * $Revision: 3421 $   $Date: 2008-12-09 10:21:11 -0600 (Tue, 09 Dec 2008) $
 * $Id: LocalEnergyEstimatorHDF.h 3421 2008-12-09 16:21:11Z mcminis2 $
 ***************************************************************************/
