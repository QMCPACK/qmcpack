//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKERLOGCOLLECTOR_H
#define QMCPLUSPLUS_WALKERLOGCOLLECTOR_H

#include "WalkerLogBuffer.h"


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
  bool logs_active;
  /// period between MC steps for data collection
  int step_period;
  /// controls verbosity of log file writes
  bool verbose;

  WalkerLogState()
  {
    reset();
    step_period = 1;
  }

  inline void reset()
  {
    logs_active = false;
    verbose     = false;
  }
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

  /// state data set by WalkerLogManager
  WalkerLogState state;

private:
  // temporary (contiguous) storage for awful ParticleAttrib<> quantities
  /// tmp storage for walker positions
  Array<WLog::Real, 2> Rtmp;
  /// tmp storage for walker spins
  Array<WLog::Real, 1> Stmp;
  /// tmp storage for walker wavefunction gradients
  Array<WLog::PsiVal, 2> Gtmp;
  /// tmp storage for walker wavefunciton laplacians
  Array<WLog::PsiVal, 1> Ltmp;

public:
  WalkerLogCollector();

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
  /// resize buffers to zero rows
  void resetBuffers();
};


} // namespace qmcplusplus

#endif
