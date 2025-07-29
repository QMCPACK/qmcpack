//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKERLOGCOLLECTOR_H
#define QMCPLUSPLUS_WALKERLOGCOLLECTOR_H

#include "WalkerLogBuffer.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{

template<class QT, class PT>
class Walker;
class ParticleSet;
class TrialWaveFunction;
class QMCHamiltonian;
using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;


/// Helper struct holding data transferred from WalkerLogManager to WalkerLogCollector following input read
struct WalkerLogState
{
  /// whether logs are active in the current driver
  bool logs_active = false;
  /// period between MC steps for data collection
  int step_period = 1;
  /// controls verbosity of log file writes
  bool verbose = false;
};


/** Crowd-level resource for walker log collection.
 *
 *    Contains data buffers for walker properties and walker particle data.
 *    Data buffers are resized to zero at the start of an MC block.
 *    Data for all walkers is collected into the buffers each MC step in an MC block.
 *    This class is not responsible for I/O.
 */
class WalkerLogCollector
{
public:
  /// MC step information for each walker throughout the MC block
  std::vector<size_t> steps;
  /// LocalEnergy information for each walker throughout the MC block
  std::vector<WLog::Real> energies;
  /// buffer containing integer walker properties
  WalkerLogBuffer<WLog::Int> walker_property_int_buffer;
  /// buffer containing real-valued walker properties
  WalkerLogBuffer<WLog::Real> walker_property_real_buffer;
  /// buffer containing per-particle walker data
  WalkerLogBuffer<WLog::Real> walker_particle_real_buffer;

  /// ParticleSet::PropertyList quantities to include
  std::unordered_set<std::string> properties_include;
  /// indices in ParticleSet::PropertyList for included quantities
  std::vector<size_t> property_indices;
  /// location of LocalEnergy in ParticleSet::PropertyList
  int energy_index;

private:
  static constexpr std::string_view my_name_{"WalkerLogCollector"};
  enum Timer
  {
    START = 0,
    COLLECT,
    CHECK_BUFFERS
  };
  static constexpr std::array<std::string_view, 3> suffixes_{"start", "collect", "check_buffers"};
  static TimerNameList_t<Timer> create_names(const std::string_view& my_name)
  {
    TimerNameList_t<Timer> timer_names;
    using namespace std::string_literals;
    std::string prefix{"WalkerLog:"s + std::string{my_name} + "::"s};
    for (std::size_t i = 0; i < suffixes_.size(); ++i)
      timer_names.push_back({static_cast<Timer>(i), prefix + std::string{suffixes_[i]}});
    return timer_names;
  }
  TimerList_t walker_log_collector_timers_;

  // temporary (contiguous) storage for awful ParticleAttrib<> quantities
  /// tmp storage for walker positions
  Array<WLog::Real, 2> Rtmp;
  /// tmp storage for walker spins
  Array<WLog::Real, 1> Stmp;
  /// tmp storage for walker wavefunction gradients
  Array<WLog::PsiVal, 2> Gtmp;
  /// tmp storage for walker wavefunciton laplacians
  Array<WLog::PsiVal, 1> Ltmp;
  /** Hopefully This is just a bundle of constructor arguments.
   *  It was a reference, so perhaps, the intention was dynamic
   *  manipulation of a group of objects state from afar.
   *  Since it hadn't yet been used this way its just a private member
   *  now. _reviewers_ do not allow it to be made a reference again.
   *
   *  If you want to do a state transform write a function, document
   *  it, call it from a sensible and obvious scope.
   */
  WalkerLogState state_;

public:
  /// constructor. The state should be given by the manager.
  WalkerLogCollector(const WalkerLogState& state);
  WalkerLogCollector(WalkerLogCollector&& other)            = default;
  WalkerLogCollector& operator=(WalkerLogCollector&& other) = default;
  /// resize buffers to zero rows at beginning of each MC block
  void startBlock();

  /// collect all data for one walker into the data buffers
  void collect(const MCPWalker& walker,
               const ParticleSet& pset,
               const TrialWaveFunction& wfn,
               const QMCHamiltonian& ham,
               int step = -1);

  /** Check that all buffers have the same number of rows.
   *    This ensures that the full data for a given walker can be reconstructed due to enforced data alignment in the buffers.
   */
  void checkBuffers();

private:
  /// shared variable setting in constructors
  void init();

  /// resize buffers to zero rows
  void resetBuffers();
};


} // namespace qmcplusplus

#endif
