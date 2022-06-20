//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov,  Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CORRELATEDLOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_CORRELATEDLOCALENERGYESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"
#include "ScalarEstimatorInputs.h"

namespace qmcplusplus
{
class QMCHamiltonian;
class TrialWaveFunction;

struct CSEnergyEstimator : public ScalarEstimatorBase
{
  //enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, PE_INDEX, KE_INDEX, LE_INDEX};
  enum
  {
    ENERGY_INDEX,
    ENERGY_SQ_INDEX,
    WEIGHT_INDEX,
    LE_INDEX
  };

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
  const CSLocalEnergyInput input_;
  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  CSEnergyEstimator(const QMCHamiltonian& h, int hcopy = 1);

  CSEnergyEstimator(CSLocalEnergyInput&& input, const QMCHamiltonian& h);

  inline RealType getUmbrellaWeight(int ipsi)
  {
    return scalars_saved[ipsi * LE_INDEX + WEIGHT_INDEX].result();
    //return d_data[ipsi*LE_INDEX+WEIGHT_INDEX];
  }

  void accumulate(const Walker_t& awalker, RealType wgt);

  inline void accumulate(const MCWalkerConfiguration& W,
                         WalkerIterator first,
                         WalkerIterator last,
                         RealType wgt) override
  {
    //accumulate  the number of times accumulation has occurred.
    //d_wgt+=last-first;
    for (; first != last; ++first)
      accumulate(**first, wgt);
  }

  inline void accumulate(const RefVector<MCPWalker>& walkers) override
  {
    for (MCPWalker& walker : walkers)
      accumulate(walker, 1.0);
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void add2Record(RecordNamedProperty<RealType>& record) override;
  void registerObservables(std::vector<ObservableHelper>& h5dec, hid_t gid) override;
  CSEnergyEstimator* clone() override;
  const std::string& getSubTypeStr() const override { return input_.get_type(); }
  void evaluateDiff();
  // CSEnergyEstimator is the main estimator for
  bool isMainEstimator() const override { return true; }
};

} // namespace qmcplusplus
#endif
