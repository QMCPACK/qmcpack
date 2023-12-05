//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HARMONIC_EXTERNAL_POTENTIAL_H
#define QMCPLUSPLUS_HARMONIC_EXTERNAL_POTENTIAL_H

#include "QMCHamiltonians/OperatorBase.h"


namespace qmcplusplus
{
struct HarmonicExternalPotential : public OperatorBase
{
  //data members
  RealType mass;
  RealType energy;
  RealType length;
  PosType center;
  const ParticleSet& Ps;

#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace sample array
  Array<TraceReal, 1>* V_sample;
#endif

  //construction/destruction
  HarmonicExternalPotential(ParticleSet& P) : Ps(P)
  {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
  }

  ~HarmonicExternalPotential() override {}

  std::string getClassName() const override { return "HarmonicExternalPotential"; }
  //unneeded interface functions
  void resetTargetParticleSet(ParticleSet& P) override {}

  //standard interface functions
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) final;

  //functions for physical (hamiltonian component) estimator
  Return_t evaluate(ParticleSet& P) override;
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy) { return evaluate(P); }

#if !defined(REMOVE_TRACEMANAGER)
  //traces interface
  void contributeParticleQuantities() override { request_.contribute_array(name_); }

  void checkoutParticleQuantities(TraceManager& tm) override
  {
    streaming_particles_ = request_.streaming_array(name_);
    if (streaming_particles_)
      V_sample = tm.checkout_real<1>(name_, Ps);
  }

  void deleteParticleQuantities() override
  {
    if (streaming_particles_)
      delete V_sample;
  }

  //  not really for interface, just collects traces
  inline Return_t evaluate_sp(ParticleSet& P);
#endif
};
} // namespace qmcplusplus
#endif
