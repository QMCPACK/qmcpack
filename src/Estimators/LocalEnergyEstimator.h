//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov,  Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_LOCALENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "ScalarEstimatorInputs.h"

namespace qmcplusplus
{
/** Class to accumulate the local energy and components
 *
 * Use Walker::Properties to accumulate Hamiltonian-related quantities.
 */
class LocalEnergyEstimator : public ScalarEstimatorBase
{
  using WP = WalkerProperties::Indexes;
  enum
  {
    ENERGY_INDEX,
    ENERGY2_INDEX,
    POTENTIAL_INDEX,
    LE_MAX
  };

  int FirstHamiltonian;
  int SizeOfHamiltonians;
  bool UseHDF5;
  const QMCHamiltonian& refH;
  const LocalEnergyInput input_;

public:
  /** constructor
   *  @param[in] h            QMCHamiltonian to define the components
   *  @param[in] use_hdf5     record local energy values in hdf5?
   */
  LocalEnergyEstimator(const QMCHamiltonian& h, bool use_hdf5);
  /** Construct from LocalEnergyInput and const reference to hamiltonian.
   *  @param[in] input     contains input parameters for LocalEnergyEstimator
   *  @param[in] ham       is taken as a local reference and used to size scalars data
   */
  LocalEnergyEstimator(LocalEnergyInput&& input, const QMCHamiltonian& ham);
  /** accumulation per walker
   *  @param awalker current walker
   *  @param wgt weight
   *
   * Weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
   */
  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
    const RealType* restrict ePtr = awalker.getPropertyBase();
    RealType wwght                = wgt * awalker.Weight;
    scalars[0](ePtr[WP::LOCALENERGY], wwght);
    scalars[1](ePtr[WP::LOCALENERGY] * ePtr[WP::LOCALENERGY], wwght);
    scalars[2](ePtr[WP::LOCALPOTENTIAL], wwght);
    for (int target = 3, source = FirstHamiltonian; target < scalars.size(); ++target, ++source)
      scalars[target](ePtr[source], wwght);
  }

  /** legacy accumulation function
   */
  inline void accumulate(const MCWalkerConfiguration& W,
                         WalkerIterator first,
                         WalkerIterator last,
                         RealType wgt) override
  {
    for (; first != last; ++first)
      accumulate(**first, wgt);
  }

  std::string getName() const override { return "LocalEnergyEstimator"; }

  void add2Record(RecordListType& record) override;
  void registerObservables(std::vector<ObservableHelper>& h5dec, hdf_archive& file) override;
  LocalEnergyEstimator* clone() override;

  /** Accumulate the hamiltonian operator values for system
   *  This is the batched version
   */
  inline void accumulate(const RefVector<MCPWalker>& walkers) override
  {
    for (MCPWalker& walker : walkers)
      accumulate(walker, 1.0);
  }

  /// LocalEnergyEstimator is the main estimator for VMC and DMC
  bool isMainEstimator() const override { return true; }
  const std::string& getSubTypeStr() const override { return input_.get_type(); }
};
} // namespace qmcplusplus
#endif
