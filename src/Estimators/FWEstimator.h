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
    
    



#ifndef QMCPLUSPLUS_FORWARDWALKING_ESTIMATOR_H
#define QMCPLUSPLUS_FORWARDWALKING_ESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include <cassert>

namespace qmcplusplus
{

/** Class to accumulate the local energy and components
 *
 * Use Walker::Properties to accumulate Hamiltonian-related quantities.
 */
class ForwardWalkingEstimator: public ScalarEstimatorBase
{

  enum {ENERGY_INDEX, POTENTIAL_INDEX, LE_MAX};

  int FirstHamiltonian;
  int SizeOfHamiltonians;
  const QMCHamiltonian& refH;

public:

  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  ForwardWalkingEstimator(QMCHamiltonian& h);

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
    RealType wwght2= wwght*wwght;//I forget, for variance, should I also square the weight? I don't think so.
    scalars[0](ePtr[LOCALENERGY],wwght);
    scalars[1](ePtr[LOCALENERGY]*ePtr[LOCALENERGY],wwght);
    scalars[2](ePtr[LOCALPOTENTIAL],wwght);
    scalars[3](ePtr[LOCALPOTENTIAL]*ePtr[LOCALPOTENTIAL],wwght);
    for(int target=4, source=FirstHamiltonian; target<scalars.size(); target+=2, ++source)
    {
      scalars[target](ePtr[source],wwght);
      scalars[target+1](ePtr[source]*ePtr[source],wwght);
    }
  }

  inline void accumulate_fw( HDF5_FW_observables& OBS, HDF5_FW_weights& WGTS, std::vector<int>& Dims )
  {
//       clear();
//       std::cout <<"Calling accumulate"<< std::endl;
//       int maxV = OBS.getLength();
//       int maxW = WGTS.getLength(step);
//       assert (maxV==maxW);
    for(int i=Dims[3]; i<(Dims[2]-Dims[0]); i++)
    {
      std::vector<int> weights;
      std::vector<double> values;
      OBS.readStep(i,values);
      WGTS.readStep(Dims[0],i,weights);
      std::vector<double>::iterator vit(values.begin());
      for(std::vector<int>::iterator wit(weights.begin()); ((wit!=weights.end())&(vit!=values.end())); wit++)
        for(int k=0; k<Dims[1]; k++,vit++)
        {
          scalars[2*k]((*vit), (*wit));
          scalars[2*k+1]((*vit)*(*vit), (*wit));
        }
    }
//       const RealType* restrict ePtr = awalker.getPropertyBase();
//       RealType wwght= wgt* awalker.Weight;
//       RealType wwght2= wwght*wwght;//I forget, for variance, should I also square the weight? I don't think so.
//       scalars[0](ePtr[LOCALENERGY],wwght);
//       scalars[1](ePtr[LOCALENERGY]*ePtr[LOCALENERGY],wwght);
//       scalars[2](ePtr[LOCALPOTENTIAL],wwght);
//       scalars[3](ePtr[LOCALPOTENTIAL]*ePtr[LOCALPOTENTIAL],wwght);
//       for(int target=4, source=FirstHamiltonian; target<scalars.size(); target+=2, ++source)
//       {
//         scalars[target](ePtr[source],wwght);
//         scalars[target+1](ePtr[source]*ePtr[source],wwght);
//       }
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
