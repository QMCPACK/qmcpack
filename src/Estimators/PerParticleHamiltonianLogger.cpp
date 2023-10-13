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

#include "PerParticleHamiltonianLogger.h"
#include <algorithm>
#include <functional>
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

namespace qmcplusplus
{
PerParticleHamiltonianLogger::PerParticleHamiltonianLogger(PerParticleHamiltonianLoggerInput&& input, int rank)
    : OperatorEstBase(DataLocality::crowd), rank_estimator_(nullptr), input_(input), rank_(rank)
{
  requires_listener_ = true;
  my_name_           = "PerParticleHamiltonianLogger";

  std::string filename("rank_" + std::to_string(rank_) + "_" + input_.get_name() + ".dat");
  rank_fstream_.open(filename, std::ios::out);
}

PerParticleHamiltonianLogger::PerParticleHamiltonianLogger(const PerParticleHamiltonianLogger& pphl, DataLocality dl)
    : OperatorEstBase(dl), rank_estimator_(const_cast<PerParticleHamiltonianLogger*>(&pphl)), input_(pphl.input_)
{
  requires_listener_ = true;
  my_name_           = pphl.name_;
  data_locality_     = dl;
}

void PerParticleHamiltonianLogger::write(CrowdLogValues& cl_values, const std::vector<long>& walker_ids)
{
  // fstream is not thread safe but it is buffered.  If the buffer isn't too small this
  // should mostly return quickly and the contention for the lock should be manageable.
  const std::lock_guard<std::mutex> lock(write_lock);
  for (auto& [component, values] : cl_values)
  {
    rank_fstream_ << "operator: " << component << '\n'; // " crowd: " << crowd_id <<
    for (int iw = 0; iw < values.size(); ++iw)
      rank_fstream_ << " walker:" << walker_ids[iw] << " " << NativePrint(values[iw]) << '\n';
  }
  if (input_.get_to_stdout())
  {
    for (auto& [component, values] : cl_values)
    {
      std::cout << component << '\n';
      for (int iw = 0; iw < values.size(); ++iw)
        std::cout << " walker: " << walker_ids[iw] << "  " << NativePrint(values[iw]) << '\n';
    }
  }
}

void PerParticleHamiltonianLogger::accumulate(const RefVector<MCPWalker>& walkers,
                                              const RefVector<ParticleSet>& psets,
                                              const RefVector<TrialWaveFunction>& wfns,
                                              RandomBase<FullPrecRealType>& rng)

{
  // The hamiltonian doesn't know the walker ID only its index in the walker elements of the crowd
  // build mapping from that index to global walker id.
  // This could change every call for DMC.
  walker_ids_.clear();
  for (MCPWalker& walker : walkers)
    walker_ids_.push_back(walker.ID);
  rank_estimator_->write(values_, walker_ids_);

  // \todo some per crowd reduction.
  //       clear log values
}

std::unique_ptr<OperatorEstBase> PerParticleHamiltonianLogger::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  auto spawn = std::make_unique<PerParticleHamiltonianLogger>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size, 0.0);
  return spawn;
}

ListenerVector<QMCTraits::RealType>::ReportingFunction PerParticleHamiltonianLogger::getLogger()
{
  auto& local_values = values_;
  return [&local_values](const int walker_index, const std::string& name, const Vector<Real>& inputV) {
    if (walker_index >= local_values[name].size())
      local_values[name].resize(walker_index + 1);
    local_values[name][walker_index] = inputV;
  };
}

void PerParticleHamiltonianLogger::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  int crowd_count = 0;
}

void PerParticleHamiltonianLogger::registerListeners(QMCHamiltonian& ham_leader)
{
  ListenerVector<Real> listener(name_, getLogger());
  QMCHamiltonian::mw_registerLocalEnergyListener(ham_leader, listener);
}

void PerParticleHamiltonianLogger::startBlock(int steps)
{
  ++block_;
  rank_fstream_ << "starting block:  " << block_ << " steps: " << steps << "\n";
}

} // namespace qmcplusplus
