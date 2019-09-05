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
 *
 * An example wavefunction component for a simple wavefunction for a helium atom.
 * It uses STO orbitals and a Pade Jastrow.  This wavefunction can be created using
 * the input file in examples/molecules/He/he_example_wf.xml.
 * There is more detail in the manual in the Development chapter, in the "Adding a wavefunction" section.
 */
namespace qmcplusplus
{
class ExampleHeComponent : public WaveFunctionComponent
{
public:

  ExampleHeComponent(const ParticleSet& ions, ParticleSet& els)
    : ions_(ions), my_table_ee_idx_(els.addTable(els, DT_SOA)), my_table_ei_idx_(els.addTable(ions, DT_SOA))
  {
    ClassName = "ExampleHeComponent";
  };

  using OptVariablesType = optimize::VariableSet;
  using PtclGrpIndexes = QMCTraits::PtclGrpIndexes;

  void checkInVariables(OptVariablesType& active) override { active.insertFrom(my_vars_); }
  void checkOutVariables(const OptVariablesType& active) override { my_vars_.getIndex(active); }
  void resetParameters(const OptVariablesType& active) override;


  void reportStatus(std::ostream& os) override {}

  void resetTargetParticleSet(ParticleSet& P) override {}

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L) override;

  void acceptMove(ParticleSet& P, int iat) override {}

  void restore(int iat) override {}

  ValueType ratio(ParticleSet& P, int iat) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void evaluateDerivatives(ParticleSet& P,
                           const OptVariablesType& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;


  void registerData(ParticleSet& P, WFBufferType& buf) override {}

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  WaveFunctionComponentPtr makeClone(ParticleSet& tpq) const override;

  bool put(xmlNodePtr cur);

  bool opt_B;
  RealType B;
  std::string ID_B;

  RealType A;
  RealType Z;

private:
  const ParticleSet& ions_;
  const int my_table_ei_idx_;
  const int my_table_ee_idx_;

  OptVariablesType my_vars_;
};

} // namespace qmcplusplus
#endif
