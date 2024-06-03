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


#ifndef QMCPLUSPLUS_WALKERTRACECOLLECTOR_H
#define QMCPLUSPLUS_WALKERTRACECOLLECTOR_H

#include "WalkerTraceBuffer.h"


namespace qmcplusplus
{

template<class QT, class PT> class Walker;
class ParticleSet;
class TrialWaveFunction;
class QMCHamiltonian;
using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;



struct WalkerTraceState 
{
  bool traces_active;
  int step_period;
  bool verbose;

  WalkerTraceState()
  {
    reset();
    step_period = 1;
  }

  inline void reset()
  {
    traces_active = false;
    verbose       = false;
  }
};


class WalkerTraceCollector
{
public:
  std::vector<size_t>     steps;
  std::vector<WTrace::Real> energies;
  WalkerTraceBuffer<WTrace::Int>  walker_property_int_buffer;
  WalkerTraceBuffer<WTrace::Real> walker_property_real_buffer;
  WalkerTraceBuffer<WTrace::Real> walker_particle_real_buffer;

  std::unordered_set<std::string> properties_include;
  std::vector<size_t>             property_indices;
  int energy_index;

  WalkerTraceState state;

private:
  // temporary (contiguous) storage for awful ParticleAttrib<> quantities
  Array<WTrace::Real  , 2> Rtmp;
  Array<WTrace::Real  , 1> Stmp;
  Array<WTrace::PsiVal, 2> Gtmp;
  Array<WTrace::PsiVal, 1> Ltmp;

public:
  WalkerTraceCollector();

  void startBlock();

  void collect(const MCPWalker& walker, const ParticleSet& pset, const TrialWaveFunction& wfn, const QMCHamiltonian& ham, int step=-1);

  void checkBuffers();

private:
  void resetBuffers();
};





} // namespace qmcplusplus

#endif
