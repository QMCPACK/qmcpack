//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPIN_DENSITY_H
#define QMCPLUSPLUS_SPIN_DENSITY_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class SpinDensity : public OperatorBase
{
public:
  using Lattice_t = ParticleSet::ParticleLayout;
  using dens_t    = std::vector<RealType>;
  using pts_t     = std::vector<PosType>;

  ParticleSet* Ptmp;

  //data members
  int nspecies;
  std::vector<int> species_size;
  std::vector<std::string> species_name;
  Lattice_t cell;
  PosType corner;
  TinyVector<int, DIM> grid;
  TinyVector<int, DIM> gdims;
  int npoints;

  //constructor/destructor
  SpinDensity(ParticleSet& P);
  ~SpinDensity() override {}

  //standard interface
  std::string getClassName() const override { return "SpinDensity"; }
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
  void reset();
  void report(const std::string& pad);
  void test(int moves, ParticleSet& P);
  Return_t test_evaluate(ParticleSet& P, int& pmin, int& pmax);
};

} // namespace qmcplusplus

#endif
