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
public:
  ChiesaCorrection(ParticleSet& ptcl, const TrialWaveFunction& psi);

  std::string getClassName() const override;

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override;

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

private:
  const TrialWaveFunction& psi_ref_;
  ParticleSet& ptcl_ref_;
};

} // namespace qmcplusplus

#endif
