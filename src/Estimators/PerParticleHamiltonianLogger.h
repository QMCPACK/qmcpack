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

namespace qmcplusplus
{

class PerParticleHamiltonianLogger : public OperatorEstBase
{
public:
  using Real = QMCTraits::RealType;
  using CrowdLogValues = std::unordered_map<std::string, std::vector<Vector<Real>>>;

  
  PerParticleHamiltonianLogger(PerParticleHamiltonianLoggerInput&& input, int rank);
  PerParticleHamiltonianLogger(const PerParticleHamiltonianLogger& other, DataLocality data_locality);

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecRealType>& rng) override;

  UPtr<OperatorEstBase> spawnCrowdClone() const override;
  void startBlock(int steps) override;

  void registerListeners(QMCHamiltonian& ham_leader) override;
  /** return lambda function to register as listener
   *  the purpose of this function is to factor out the production of the lambda for unit testing
   *  \param[out] values
   */
  ListenerVector<Real>::ReportingFunction getLogger();

  void collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators) override;

  void write(CrowdLogValues& values, const std::vector<long>& walkers_ids);

  int get_block() { return block_; }
private:
  bool crowd_clone = false;
  PerParticleHamiltonianLogger  * const rank_estimator_;
  PerParticleHamiltonianLoggerInput input_;
  int rank_;
  CrowdLogValues values_;
  std::vector<long> walker_ids_;
  const std::string name_{"PerParticleHamiltonianLogger"};
  std::fstream rank_fstream_;
  std::mutex write_lock;
  int block_ = 0;
};
  
}

#endif
