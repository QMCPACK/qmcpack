//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PER_PARTICLE_HAMILTONIAN_LOGGER_H
#define QMCPLUSPLUS_PER_PARTICLE_HAMILTONIAN_LOGGER_H
#include "PerParticleHamiltonianLoggerInput.h"
#include <string>
#include <fstream>
#include <mutex>
#include <unordered_map>
#include <optional>
#include <StdRandom.h>
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCHamiltonians/Listener.hpp"
#include "type_traits/OptionalRef.hpp"
namespace qmcplusplus
{

class PerParticleHamiltonianLogger : public OperatorEstBase
{
public:
  using Real                       = QMCTraits::RealType;
  using CrowdLogValues             = std::unordered_map<std::string, std::vector<Vector<Real>>>;
  using LabeledCrowdLogValues = std::unordered_map<std::string, CrowdLogValues&>;

  PerParticleHamiltonianLogger(PerParticleHamiltonianLoggerInput&& input, int rank);
  PerParticleHamiltonianLogger(PerParticleHamiltonianLogger& other, DataLocality data_locality);

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  const RefVector<QMCHamiltonian>& hams,
                  RandomBase<FullPrecRealType>& rng) override;

  UPtr<OperatorEstBase> spawnCrowdClone() const override;
  void startBlock(int steps) override;

  void registerListeners(QMCHamiltonian& ham_leader) override;

  /** return lambda function to register as listener
   *  the purpose of this function is to factor out the production of the lambda for unit testing
   *  param values, the default just uses the local_energy_values_
   */
  ListenerVector<Real>::ReportingFunction getLogger();

  void collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators) override;

  void write(LabeledCrowdLogValues& label_log_values, const std::vector<long>& walkers_ids);

  /** This function is supplied for testing it sums over all values currently stored by the
   *  per particle logger.  This is useful in unit testing some of the other estimators
   *  that make use of Listeners.
   */
  Real sumOverAll();

  Real sumOverSome(std::vector<std::string> which_listeners);

  int get_block() { return block_; }

private:
  ListenerVector<Real>::ReportingFunction getLogger(CrowdLogValues& values);

  bool crowd_clone = false;
  const OptionalRef<PerParticleHamiltonianLogger> rank_estimator_;
  PerParticleHamiltonianLoggerInput input_;
  int rank_;
  CrowdLogValues local_potential_values_;
  CrowdLogValues local_ion_potential_values_;
  CrowdLogValues local_energy_values_;
  CrowdLogValues kinetic_energy_values_;

  /// To allow iteration over different reporting "queues" with labels
  LabeledCrowdLogValues labeled_crowd_log_values_;

  std::vector<long> walker_ids_;
  const std::string name_{"PerParticleHamiltonianLogger"};
  /// rank owned fstream
  std::fstream rank_fstream_;
  /// use it to prevent race when writing per-crowd.
  std::mutex write_lock;
  int block_ = 0;
};

} // namespace qmcplusplus

#endif
