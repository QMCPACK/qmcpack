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
#ifndef QMCPLUSPLUS_RNENERGYESTIMATOR_H
#define QMCPLUSPLUS_RNENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{

/** Class to accumulate the local energy and components
 *
 * Use Walker::Properties to accumulate Hamiltonian-related quantities.
 */
class ReleasedNodeEnergyEstimator: public ScalarEstimatorBase
{

  enum {ENERGY_INDEX, ENERGY2_INDEX, POTENTIAL_INDEX, FE_NOW, LE_MAX};

  int FirstHamiltonian;
  int SizeOfHamiltonians;
  const QMCHamiltonian& refH;
  int N_rn, Smax;

public:

  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  ReleasedNodeEnergyEstimator(QMCHamiltonian& h, int Sm);

  /** accumulation per walker
   * @param awalker current walker
   * @param wgt weight
   *
   * Weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
   */
  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
    const RealType* restrict ePtr = awalker.getPropertyBase();
    RealType wwght= wgt* awalker.Weight;
    scalars[0](ePtr[LOCALENERGY],wwght);
    scalars[1](ePtr[LOCALENERGY]*ePtr[LOCALENERGY],wwght);
    scalars[2](ePtr[LOCALPOTENTIAL],wwght);
    scalars[3](ePtr[ALTERNATEENERGY],wwght*std::abs(awalker.ReleasedNodeWeight));
    int rni=awalker.ReleasedNodeAge;
    if ((rni<N_rn)&&(rni>=0))
    {
      scalars[4+3*rni].add(awalker.ReleasedNodeWeight*ePtr[ALTERNATEENERGY]*wwght);
      scalars[5+3*rni].add(awalker.ReleasedNodeWeight*ePtr[ALTERNATEENERGY]*ePtr[ALTERNATEENERGY]*wwght);
      scalars[6+3*rni].add(awalker.ReleasedNodeWeight*wwght);
    }
    if (rni==0)
      for(int target=4+3*N_rn, source=FirstHamiltonian; target<scalars.size(); ++target, ++source)
        scalars[target](ePtr[source],wwght);
  }

  /*@{*/
  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first != last; ++first)
      accumulate(**first,wgt);
  }
  void add2Record(RecordListType& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid) {}
  ScalarEstimatorBase* clone();
  /*@}*/
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4163 $   $Date: 2009-08-31 05:47:46 -0500 (Mon, 31 Aug 2009) $
 * $Id: ReleasedNodeEnergyEstimator.h 4163 2009-08-31 10:47:46Z jmcminis $
 ***************************************************************************/
