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
class ExampleHeComponent : public WaveFunctionComponent, OptimizableObject
{
public:
  ExampleHeComponent(const ParticleSet& ions, ParticleSet& els)
      : OptimizableObject("example"),
        ions_(ions),
        my_table_ee_idx_(els.addTable(els, DTModes::NEED_TEMP_DATA_ON_HOST | DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
        my_table_ei_idx_(els.addTable(ions, DTModes::NEED_VP_FULL_TABLE_ON_HOST)){};

  using OptVariablesType = optimize::VariableSet;
  using PtclGrpIndexes   = QMCTraits::PtclGrpIndexes;

  std::string getClassName() const override { return "ExampleHeComponent"; }
  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override { opt_obj_refs.push_back(*this); }
  bool isOptimizable() const override { return true; }
  void checkInVariablesExclusive(OptVariablesType& active) override { active.insertFrom(my_vars_); }
  void checkOutVariables(const OptVariablesType& active) override { my_vars_.getIndex(active); }
  void resetParametersExclusive(const OptVariablesType& active) override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override {}

  void restore(int iat) override {}

  PsiValueType ratio(ParticleSet& P, int iat) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void evaluateDerivatives(ParticleSet& P,
                           const OptVariablesType& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override;


  void registerData(ParticleSet& P, WFBufferType& buf) override {}

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tpq) const override;

  bool put(xmlNodePtr cur);

  bool opt_B;
  RealType B;
  std::string ID_B;

  RealType A;
  RealType Z;

private:
  const ParticleSet& ions_;
  const int my_table_ee_idx_;
  const int my_table_ei_idx_;

  OptVariablesType my_vars_;
};

} // namespace qmcplusplus
#endif
