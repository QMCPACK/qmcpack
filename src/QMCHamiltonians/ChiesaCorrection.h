//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CHIESA_CORRECTION_H
#define QMCPLUSPLUS_CHIESA_CORRECTION_H
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
class ChiesaCorrection : public OperatorBase
{
private:
  const TrialWaveFunction& psi_ref;
  ParticleSet& ptcl_ref;

public:
  ChiesaCorrection(ParticleSet& ptcl, const TrialWaveFunction& psi) : psi_ref(psi), ptcl_ref(ptcl) {}

  std::string getClassName() const override { return "ChiesaCorrection"; }
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

#ifdef QMC_CUDA
  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override;
#endif

  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "Chiesa correction: " << ptcl_ref.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;
};

} // namespace qmcplusplus

#endif
