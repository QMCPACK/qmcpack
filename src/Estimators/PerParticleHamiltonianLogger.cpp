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

using namespace std::string_literals;

PerParticleHamiltonianLogger::PerParticleHamiltonianLogger(PerParticleHamiltonianLoggerInput&& input, int rank)
    : OperatorEstBase(DataLocality::crowd), input_(input), rank_(rank)
{
  labeled_crowd_log_values_ = {{"local_potential"s, local_potential_values_},
                               {"local_energy"s, local_energy_values_},
                               {"ion_potential"s, local_ion_potential_values_},
                               {"kinetic_energy"s, kinetic_energy_values_}};

  requires_listener_ = true;
  my_name_           = "PerParticleHamiltonianLogger";

  std::string filename("rank_" + std::to_string(rank_) + "_" + input_.get_name() + ".dat");
  rank_fstream_.open(filename, std::ios::out);
}

PerParticleHamiltonianLogger::PerParticleHamiltonianLogger(PerParticleHamiltonianLogger& pphl, DataLocality dl)
    : OperatorEstBase(dl), rank_estimator_(makeOptionalRef(pphl)), input_(pphl.input_)
{
  labeled_crowd_log_values_ = {{"local_potential"s, local_potential_values_},
                               {"local_energy"s, local_energy_values_},
                               {"ion_potential"s, local_ion_potential_values_},
                               {"kinetic_energy"s, kinetic_energy_values_}};

  requires_listener_ = true;
  my_name_           = pphl.name_;
  data_locality_     = dl;
}

void PerParticleHamiltonianLogger::write(LabeledCrowdLogValues& labeled_cl_values, const std::vector<long>& walker_ids)
{
  // fstream is not thread safe but it is buffered.  If the buffer isn't too small this
  // should mostly return quickly and the contention for the lock should be manageable.
  const std::lock_guard<std::mutex> lock(write_lock);
  for (auto& [label, cl_values] : labeled_cl_values)
  {
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
}

void PerParticleHamiltonianLogger::accumulate(const RefVector<MCPWalker>& walkers,
                                              const RefVector<ParticleSet>& psets,
                                              const RefVector<TrialWaveFunction>& wfns,
                                              const RefVector<QMCHamiltonian>& hams,
                                              RandomBase<FullPrecRealType>& rng)

{
  // The hamiltonian doesn't know the walker ID only its index in the walker elements of the crowd
  // build mapping from that index to global walker id.
  // This could change every call for DMC.
  walker_ids_.clear();
  for (MCPWalker& walker : walkers)
    walker_ids_.push_back(walker.getWalkerID());
  rank_estimator_->get().write(labeled_crowd_log_values_, walker_ids_);

  // \todo some per crowd reduction.
  //       clear log values
}

PerParticleHamiltonianLogger::Real PerParticleHamiltonianLogger::sumOverSome(std::vector<std::string> which_listeners)
{
  Real sum{0};
  for (auto& listener_label : which_listeners)
  {
    const auto& cl_values = labeled_crowd_log_values_.at(listener_label);
    for (auto& [component, values] : cl_values)
      for (auto& a_vector : values)
        sum += std::accumulate(a_vector.begin(), a_vector.end(), 0.0);
  }
  return sum;
}

PerParticleHamiltonianLogger::Real PerParticleHamiltonianLogger::sumOverAll()
{
  Real sum{0};
  for (auto& [label, cl_values] : labeled_crowd_log_values_)
  {
    for (auto& [component, values] : cl_values)
      for (auto& a_vector : values)
        sum += std::accumulate(a_vector.begin(), a_vector.end(), 0.0);
  }
  return sum;
}

std::unique_ptr<OperatorEstBase> PerParticleHamiltonianLogger::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  auto spawn = std::make_unique<PerParticleHamiltonianLogger>(const_cast<PerParticleHamiltonianLogger&>(*this),
                                                              spawn_data_locality);
  spawn->get_data().resize(data_size, 0.0);
  return spawn;
}

ListenerVector<QMCTraits::RealType>::ReportingFunction PerParticleHamiltonianLogger::getLogger()
{
  return getLogger(local_energy_values_);
}

ListenerVector<QMCTraits::RealType>::ReportingFunction PerParticleHamiltonianLogger::getLogger(CrowdLogValues& values)
{
  // \todo Now that I pass in the values this may not be necessary
  auto& local_values = values;
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
  ListenerVector<Real> local_energy_listener(name_, getLogger(local_energy_values_));
  ListenerVector<Real> local_potential_listener(name_, getLogger(local_potential_values_));
  ListenerVector<Real> local_ion_potential_listener(name_, getLogger(local_ion_potential_values_));
  ListenerVector<Real> kinetic_energy_listener(name_, getLogger(kinetic_energy_values_));
  QMCHamiltonian::mw_registerLocalEnergyListener(ham_leader, local_energy_listener);
  QMCHamiltonian::mw_registerLocalPotentialListener(ham_leader, local_potential_listener);
  QMCHamiltonian::mw_registerLocalIonPotentialListener(ham_leader, local_ion_potential_listener);
  QMCHamiltonian::mw_registerKineticListener(ham_leader, kinetic_energy_listener);
}

void PerParticleHamiltonianLogger::startBlock(int steps)
{
  ++block_;
  rank_fstream_ << "starting block:  " << block_ << " steps: " << steps << "\n";
}

} // namespace qmcplusplus
