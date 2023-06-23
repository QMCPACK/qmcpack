//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COMPOSITE_SPOSET_H
#define QMCPLUSPLUS_COMPOSITE_SPOSET_H

#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
class CompositeSPOSet : public SPOSet
{
public:
  ///component SPOSets
  std::vector<std::unique_ptr<SPOSet>> components;
  ///temporary storage for values
  std::vector<ValueVector> component_values;
  ///temporary storage for gradients
  std::vector<GradVector> component_gradients;
  ///temporary storage for laplacians
  std::vector<ValueVector> component_laplacians;
  ///store the precomputed offsets
  std::vector<int> component_offsets;

  CompositeSPOSet(const std::string& my_name);
  CompositeSPOSet(const CompositeSPOSet& other);
  ~CompositeSPOSet() override;

  std::string getClassName() const override { return "CompositeSPOSet"; }

  ///add a sposet component to this composite sposet
  void add(std::unique_ptr<SPOSet> component);

  ///print out component info
  void report();

  //SPOSet interface methods
  ///size is determined by component sposets and nothing else
  inline void setOrbitalSetSize(int norbs) override {}

  std::unique_ptr<SPOSet> makeClone() const override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  ///unimplemented functions call this to abort
  inline void not_implemented(const std::string& method)
  {
    APP_ABORT("CompositeSPOSet::" + method + " has not been implemented");
  }

  //methods to be implemented in the future (possibly)
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& ddlogdet) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& ddlogdet,
                            GGGMatrix& dddlogdet) override;
};

struct CompositeSPOSetBuilder : public SPOSetBuilder
{
  CompositeSPOSetBuilder(Communicate* comm, const SPOSetBuilderFactory& factory)
      : SPOSetBuilder("Composite", comm), sposet_builder_factory_(factory)
  {}

  //SPOSetBuilder interface
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

  std::unique_ptr<SPOSet> createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input) override;

  /// reference to the sposet_builder_factory
  const SPOSetBuilderFactory& sposet_builder_factory_;
};
} // namespace qmcplusplus

#endif
