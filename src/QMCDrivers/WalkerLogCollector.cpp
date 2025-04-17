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


#include "WalkerLogCollector.h"

#include "Particle/Walker.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"


namespace qmcplusplus
{

using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

WalkerLogCollector::WalkerLogCollector(const WalkerLogState& state)
    : walker_log_collector_timers_(getGlobalTimerManager(), create_names(my_name_), timer_level_medium), state_(state)
{
  init();
}


void WalkerLogCollector::init()
{
  properties_include = {"R2Accepted", "R2Proposed", "LocalEnergy", "LocalPotential", "Kinetic",
                        "ElecElec",   "ElecIon",    "LocalECP",    "NonLocalECP"};
  energy_index       = -1;
  // empty walker steps and energy vectors for the MC block
  steps.resize(0);
  energies.resize(0);
  // label the buffers for HDF file write
  walker_property_int_buffer.label  = "walker_property_int";
  walker_property_real_buffer.label = "walker_property_real";
  walker_particle_real_buffer.label = "walker_particle_real";
}


void WalkerLogCollector::startBlock()
{
  if (!state_.logs_active)
    return; // no-op for driver if logs are inactive

  ScopedTimer timer(walker_log_collector_timers_[Timer::START]);

  if (state_.verbose)
    app_log() << "WalkerLogCollector::startBlock " << std::endl;
  resetBuffers(); // resize buffers to zero rows
}


void WalkerLogCollector::collect(const MCPWalker& walker,
                                 const ParticleSet& pset,
                                 const TrialWaveFunction& wfn,
                                 const QMCHamiltonian& ham,
                                 int step)
{
  if (!state_.logs_active)
    return; // no-op for driver if logs are inactive

  ScopedTimer timer(walker_log_collector_timers_[Timer::COLLECT]);

  // only collect walker data at steps matching the period (default 1)
  int current_step = (step == -1) ? pset.current_step : step;
  if (current_step % state_.step_period != 0)
    return;

  auto& bsi = walker_property_int_buffer;
  auto& bsr = walker_property_real_buffer;
  auto& bar = walker_particle_real_buffer;

  // collect per-particle walker quantities
  size_t nparticles = walker.R.size();
  size_t ndim       = walker.R[0].size();
  //   per-particle positions (walker.R)
  Rtmp.resize(nparticles, ndim);
  for (size_t p = 0; p < nparticles; ++p)
    for (size_t d = 0; d < ndim; ++d)
      Rtmp(p, d) = (WLog::Real)walker.R[p][d];
  bar.collect("R", Rtmp);
  //   per-particle "spin" (walker.spins)
  if (pset.isSpinor())
  {
    Stmp.resize(nparticles);
    for (size_t p = 0; p < nparticles; ++p)
      Stmp(p) = (WLog::Real)walker.spins[p];
    bar.collect("S", Stmp);
  }
  //   per-particle gradient(log(psi)) (pset.G)
  Gtmp.resize(nparticles, ndim);
  for (size_t p = 0; p < nparticles; ++p)
    for (size_t d = 0; d < ndim; ++d)
      Gtmp(p, d) = (WLog::PsiVal)pset.G[p][d];
  bar.collect("G", Gtmp);
  //   per-particle laplacian(log(psi)) (pset.L)
  Ltmp.resize(nparticles);
  for (size_t p = 0; p < nparticles; ++p)
    Ltmp(p) = (WLog::PsiVal)pset.L[p];
  bar.collect("L", Ltmp);
  bar.resetCollect();

  // collect integer walker properties
  bsi.collect("step", (WLog::Int)current_step);
  bsi.collect("id", (WLog::Int)walker.getWalkerID());
  bsi.collect("parent_id", (WLog::Int)walker.getParentID());
  bsi.collect("age", (WLog::Int)walker.Age);
  bsi.resetCollect();

  // collect real walker properties
  bsr.collect("weight", (WLog::Real)walker.Weight);
  bsr.collect("multiplicity", (WLog::Real)walker.Multiplicity);
  bsr.collect("logpsi", (WLog::Real)wfn.getLogPsi());
  bsr.collect("phase", (WLog::Real)wfn.getPhase());
  //    from PropertyList
  if (bsr.first_collect)
  {
    for (size_t n = 0; n < pset.PropertyList.size(); ++n)
    {
      auto& name  = pset.PropertyList.Names[n];
      auto& value = walker.Properties(0, n);
      if (properties_include.find(name) != properties_include.end())
      {
        bsr.collect(name, (WLog::Real)value);
        property_indices.push_back(n);
      }
      if (name == "LocalEnergy")
        energy_index = n;
    }
    if (energy_index < 0)
      throw std::runtime_error("LogCollector::collect  energy_index must not be negative");
  }
  else
    for (auto n : property_indices)
    {
      auto& name  = pset.PropertyList.Names[n];
      auto& value = walker.Properties(0, n);
      bsr.collect(name, (WLog::Real)value);
    }
  //    nodal proximity measures
  auto& gv      = Gtmp.storage();
  WLog::Real dr = std::numeric_limits<WLog::Real>::max();
  for (size_t n = 0; n < gv.size(); ++n)
    dr = std::min(dr, std::abs(1. / std::real(gv[n])));
  auto dr_node_min    = dr;
  WLog::Real dlogpsi2 = 0.;
  for (size_t n = 0; n < gv.size(); ++n)
    dlogpsi2 += std::real(gv[n]) * std::real(gv[n]);
  WLog::Real dphase2 = 0.;
  for (size_t n = 0; n < gv.size(); ++n)
    dphase2 += std::imag(gv[n]) * std::imag(gv[n]);
  bsr.collect("dr_node_min", dr_node_min);
  bsr.collect("dlogpsi2", dlogpsi2); // dr_node = 1/sqrt(dlogpsi2)
  bsr.collect("dphase2", dphase2);   // dr_node = 1/sqrt(dphase2)
  bsr.resetCollect();

  // save the energy of this walker
  steps.push_back((size_t)current_step);
  energies.push_back((WLog::Real)walker.Properties(0, energy_index));
}


void WalkerLogCollector::resetBuffers()
{
  if (state_.verbose)
    app_log() << "WalkerLogCollector::reset_buffers" << std::endl;
  // resize all buffers to zero rows
  walker_property_int_buffer.resetBuffer();
  walker_property_real_buffer.resetBuffer();
  walker_particle_real_buffer.resetBuffer();
  // similar for step/energy vectors
  steps.resize(0);
  energies.resize(0);
}


void WalkerLogCollector::checkBuffers()
{
  ScopedTimer timer(walker_log_collector_timers_[Timer::CHECK_BUFFERS]);

  if (state_.verbose)
    app_log() << "WalkerLogCollector::checkBuffers" << std::endl;
  size_t nrows       = walker_property_int_buffer.nrows();
  auto prop_real_bad = walker_property_real_buffer.nrows() != nrows;
  auto part_real_bad = walker_particle_real_buffer.nrows() != nrows;
  auto steps_bad     = steps.size() != nrows;
  auto energies_bad  = energies.size() != nrows;
  auto any_bad       = prop_real_bad || part_real_bad || steps_bad || energies_bad;

  if (prop_real_bad)
    app_log() << "WalkerLogCollector::checkBuffers  walker_property_real_buffer row count does not match\n";
  if (part_real_bad)
    app_log() << "WalkerLogCollector::checkBuffers  walker_particle_real_buffer row count does not match\n";
  if (steps_bad)
    app_log() << "WalkerLogCollector::checkBuffers  steps entry count does not match\n";
  if (energies_bad)
    app_log() << "WalkerLogCollector::checkBuffers  energies entry count does not match\n";
  if (any_bad)
    throw std::runtime_error("WalkerLogCollector::checkBuffers  buffer lengths do not match");
}


} // namespace qmcplusplus
