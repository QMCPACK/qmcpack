//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_WFMCONLYESTIMATOR_H
#define QMCPLUSPLUS_WFMCONLYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{

class WFMCOnlyEstimator: public ScalarEstimatorBase
{

  //typedef PooledData<T>                            BufferType;
  //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, POTENTIAL_INDEX, LE_MAX};
  enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

  //int LocalPotentialIndex;
  int FirstHamiltonian;
  int SizeOfHamiltonians;
  int PsiIndex;

  ///vector to contain the names of all the constituents of the local energy
  std::vector<std::string> elocal_name;

public:

  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  WFMCOnlyEstimator(QMCHamiltonian& h);

  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
    const RealType* restrict ePtr = awalker.getPropertyBase();
    ///weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
    RealType wwght= wgt* awalker.Weight;
//       RealType wwght= wgt;
    scalars[0](ePtr[LOCALENERGY],wwght);
    scalars[1](ePtr[LOCALPOTENTIAL],wwght);
    int target=2;
    int source;
    for(source=FirstHamiltonian; source<(FirstHamiltonian+SizeOfHamiltonians);
        ++target, ++source)
    {
      scalars[target](ePtr[source],wwght);
    }
    scalars[ target](wwght*ePtr[PsiIndex+FirstHamiltonian],1);
    scalars[ target+1](awalker.Weight,1);
    scalars[ target+2](ePtr[LOCALENERGY],1.0);
  }

  /*@{*/
  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first!=last; ++first)
      std::accumulate(**first,wgt);
  }
  void add2Record(RecordListType& record);
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/

};
}
#endif
