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
    
    



#ifndef QMCPLUSPLUS_COLLECTABLES_ESTIMATOR_H
#define QMCPLUSPLUS_COLLECTABLES_ESTIMATOR_H
#include <Estimators/ScalarEstimatorBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>

namespace qmcplusplus
{

/** Handle an ensemble average of Hamiltonian components
*/
class CollectablesEstimator: public ScalarEstimatorBase
{
  ///save the reference hamiltonian
  const QMCHamiltonian& refH;

public:
  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  CollectablesEstimator(QMCHamiltonian& h);
  //LocalEnergyEstimatorHDF(const LocalEnergyEstimatorHDF& est);

  /** implement virtual function
  */
  CollectablesEstimator* clone();

  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  void add2Record(RecordListType& record);
  /** do nothing with accumulate */
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last , RealType wgt)
  { }
  /** accumulate the collectables */
  inline void accumulate_all(const MCWalkerConfiguration::Buffer_t& data
                             , RealType wgt)
  {
    for(int i=0; i<data.size(); ++i)
      scalars[i](data[i],wgt);
  }
};
}
#endif
