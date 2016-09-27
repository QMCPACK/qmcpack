//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_CORRELATEDLOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_CORRELATEDLOCALENERGYESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"
#include "QMCDrivers/SpaceWarp.h"

namespace qmcplusplus
{

class QMCHamiltonian;
class TrialWaveFunction;

struct CSEnergyEstimator: public ScalarEstimatorBase
{

  //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, PE_INDEX, KE_INDEX, LE_INDEX};
  enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, LE_INDEX};

  ///number of correlated systems
  int NumCopies;
  ///index of the starting Hamiltonian component
  int FirstHamiltonian;
  ///index of the ending Hamiltonian component
  int LastHamiltonian;
  ///save umbrella weights
  std::vector<RealType> uweights;
  ///temporary dat
  Matrix<RealType> tmp_data;
  ///name of hamiltonian components
  std::vector<std::string> h_components;
  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  CSEnergyEstimator(QMCHamiltonian& h, int hcopy=1);

  inline RealType getUmbrellaWeight(int ipsi)
  {
    return scalars_saved[ipsi*LE_INDEX+WEIGHT_INDEX].result();
    //return d_data[ipsi*LE_INDEX+WEIGHT_INDEX];
  }


  void accumulate(const Walker_t& awalker, RealType wgt);

  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    //accumulate  the number of times accumulation has occurred.
    //d_wgt+=last-first;
    for(; first != last; ++first)
      accumulate(**first,wgt);
  }
  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();

  void evaluateDiff();
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: CSEnergyEstimator.h 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/
