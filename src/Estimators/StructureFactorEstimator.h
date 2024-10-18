//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/{SkEstimator.h, SkAllEstimator.h}
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_STRUCTUREFACTORESTIMATOR_H
#define QMCPLUSPLUS_STRUCTUREFACTORESTIMATOR_H

#include "OperatorEstBase.h"
#include "StructureFactorInput.h"

namespace qmcplusplus
{

class StructureFactorEstimator : public OperatorEstBase
{
public:
  using QMCT         = QMCTraits;
  using Real         = QMCT::RealType;
  using FullPrecReal = QMCT::FullPrecRealType;
  using KPt          = TinyVector<Real, QMCTraits::DIM>;
  
  StructureFactorEstimator(StructureFactorInput& sfi,
                           ParticleSet& pset_ions,
                           ParticleSet& pset_elec,
                           DataLocality data_locality = DataLocality::crowd);

  /** accumulate 1 or more walkers of EnergyDensity samples
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  const RefVector<QMCHamiltonian>& hams,
                  RandomBase<FullPrecReal>& rng) override;

  /** start block entry point
   */
  void startBlock(int steps) override;

  UPtr<OperatorEstBase> spawnCrowdClone() const override;

  void registerOperatorEstimator(hdf_archive& file) override;
  void write(hdf_archive& file);
  void collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators) override;
private:
  int elec_species_;
  ParticleSet& elns_;
  int ion_species_;
  ParticleSet& ions_;

  /// number of k points
  long long num_kpoints_;

  // std::vector<int> kshell_degeneracy_;
  /// kpts which belong to the ith-shell [kshell[i], kshell[i+1])
  std::vector<int> kshell_offsets_;

  /// All the following are indexed by kshell
  /// instantaneous structure factor for a k shell
  std::vector<Real> kmags_;

  std::vector<Real> one_over_degeneracy_kshell_;

  /// Accumulated, its clearer to do it this way that use the base class data_ but means we don't use that base class infrastructure
  Vector<Real> sfk_e_e_;
  Vector<std::complex<Real>> rhok_e_;

  /*@{
   *  work area to reduce over electron species, they should probably just be complex as well.
   *  but particle set stores the real and imaginary parts as two real vectors.
   */  
  Vector<Real> rhok_tot_r_;
  Vector<Real> rhok_tot_i_;
  /*@}*/
};

} // namespace qmcplusplus

#endif
