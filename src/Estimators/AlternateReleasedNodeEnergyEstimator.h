//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_ALTRNENERGYESTIMATOR_H
#define QMCPLUSPLUS_ALTRNENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{

/** Class to accumulate the local energy and components
 *
 * Use Walker::Properties to accumulate Hamiltonian-related quantities.
 */
class AlternateReleasedNodeEnergyEstimator: public ScalarEstimatorBase
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
  AlternateReleasedNodeEnergyEstimator(QMCHamiltonian& h, int Sm);

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
    scalars[3](ePtr[ALTERNATEENERGY],wwght*awalker.ReleasedNodeWeight);
    int rni=awalker.ReleasedNodeAge;
    if ((rni<N_rn)&&(rni>=0))
    {
      scalars[4+3*rni].add(awalker.ReleasedNodeWeight*ePtr[ALTERNATEENERGY]*wwght);
      scalars[5+3*rni].add(awalker.ReleasedNodeWeight*ePtr[ALTERNATEENERGY]*ePtr[ALTERNATEENERGY]*wwght);
//         scalars[4+3*rni].add(awalker.ReleasedNodeWeight*ePtr[LOCALENERGY]*wwght);
//         scalars[5+3*rni].add(awalker.ReleasedNodeWeight*ePtr[LOCALENERGY]*ePtr[LOCALENERGY]*wwght);
      scalars[6+3*rni].add(awalker.ReleasedNodeWeight*wwght);
    }
    for(int target=4+3*N_rn, source=FirstHamiltonian; target<scalars.size();
        ++target, ++source)
      scalars[target](ePtr[source],wwght);
  }

  /*@{*/
  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first != last; ++first)
      std::accumulate(**first,wgt);
  }
  void add2Record(RecordListType& record);
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid) {}
  ScalarEstimatorBase* clone();
  /*@}*/
};
}
#endif
