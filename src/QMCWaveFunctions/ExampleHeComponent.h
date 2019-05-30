//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EXAMPLEHECOMPONENT_H
#define QMCPLUSPLUS_EXAMPLEHECOMPONENT_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"

/**@file ExampleHeComponent.h
 */
namespace qmcplusplus
{
class ExampleHeComponent : public WaveFunctionComponent
{
public:
  const ParticleSet& ions;
  int myTableID;

  ExampleHeComponent(const ParticleSet& ions_, ParticleSet& els) : ions(ions_)
  {
    myTableID = els.addTable(ions, DT_SOA);
  }


  typedef optimize::VariableSet opt_variables_type;

  opt_variables_type myVars;

  void checkInVariables(opt_variables_type& active) override { active.insertFrom(myVars); }
  void checkOutVariables(const opt_variables_type& active) override { myVars.getIndex(active); }
  void resetParameters(const opt_variables_type& active) override;


  void reportStatus(std::ostream& os) override {}

  void resetTargetParticleSet(ParticleSet& P) override {}

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L) override;

  void acceptMove(ParticleSet& P, int iat) override {}

  void restore(int iat) override {}

  ValueType ratio(ParticleSet& P, int iat) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi) override;


  void registerData(ParticleSet& P, WFBufferType& buf) override {}

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  WaveFunctionComponentPtr makeClone(ParticleSet& tpq) const override;

  bool put(xmlNodePtr cur);


  bool Opt_B;
  RealType B;
  std::string ID_B;

  RealType A;
  RealType Z;
};

} // namespace qmcplusplus
#endif
