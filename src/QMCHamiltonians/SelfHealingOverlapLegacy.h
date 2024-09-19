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


#ifndef QMCPLUSPLUS_SHOVERLAPLEGACY__H
#define QMCPLUSPLUS_SHOVERLAPLEGACY__H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
/** Class that collects MSD coefficient values via the Self-Healing overlap
 *    (legacy driver version) 
 */
class SelfHealingOverlapLegacy : public OperatorBase
{
public:
  using LatticeType = PtclOnLatticeTraits::ParticleLayout;
  using RealType    = QMCTraits::RealType;
  using ComplexType = QMCTraits::ComplexType;
  using ValueType   = QMCTraits::ValueType;
  using PosType     = QMCTraits::PosType;

  //data members
  size_t ncoef;
  TrialWaveFunction& psi_ref;
  Vector<ValueType> det_ratios;

  //constructor/destructor
  SelfHealingOverlapLegacy(TrialWaveFunction& wfn);
  ~SelfHealingOverlapLegacy() override {}

  //standard interface
  std::string getClassName() const override { return "SelfHealingOverlapLegacy"; }
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) final;
  bool put(xmlNodePtr cur) override;
  Return_t evaluate(ParticleSet& P) override;

  //required for Collectables interface
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P) override {}
  void setObservables(PropertySetType& plist) override {}
  void setParticlePropertyList(PropertySetType& plist, int offset) override {}
#if !defined(REMOVE_TRACEMANAGER)
  void checkout_scalar_arrays(TraceManager& tm) {}
  void collect_scalar_samples() {}
  void delete_scalar_arrays() {}
#endif

  //obsolete?
  bool get(std::ostream& os) const override { return false; }

  //local functions

};
} // namespace qmcplusplus

#endif
