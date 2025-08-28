//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/PairCorrEstimator.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_PAIRCORRELATIONESTIMATOR_H
#define QMCPLUSPLUS_PAIRCORRELATIONESTIMATOR_H

#include "PairCorrelationInput.h"
#include <ParticleSetPool.h>
#include "OperatorEstBase.h"

namespace qmcplusplus
{
/** gofr estimator
 *
 * Compute pair correlation function for the target particle set and optionally any source particles
 */
class PairCorrelationEstimator : public OperatorEstBase
{
public:
  using QMCT         = QMCTraits;
  using Real         = QMCT::RealType;
  using FullPrecReal = QMCT::FullPrecRealType;

  /** This is the Constructor called by the application
   *  the Pair correlation estimator can specify arbitrary source and
   *  target particle set names.
   */
  PairCorrelationEstimator(const PairCorrelationInput& pci,
                           const PSPool& psp,
                           DataLocality data_locality = DataLocality::crowd);

  PairCorrelationEstimator(const PairCorrelationEstimator& pce, DataLocality dl);

  /** accumulate 1 or more walkers of PairCorrelationEstimator
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  const RefVector<QMCHamiltonian>& hams,
                  RandomBase<FullPrecReal>& rng) override;

  void normalize(Real invToWgt) override;

  /** start block entry point
   */
  void startBlock(int steps) override;

  UPtr<OperatorEstBase> spawnCrowdClone() const override;

  void packData(PooledData<Real>& buffer) const override;
  void unpackData(PooledData<Real>& buffer) override;
  /** Write the one time datasets for Estimator
   *  Its calls to `file` create the root node for the estimator
   *  if it doesn't yet exist.
   *  this estimator does not use ObservableHelper this does not
   *  update any other state. and `file` is returned to whatever node
   *  position it had before calling this function.
   */
  void registerOperatorEstimator(hdf_archive& file) override;

  /** Write node to hdf5 file output
   *  The written node has the following structure
   *
   *  `file` current hdf node
   *    |
   *    + --name(from input)
   *         |
   *         +-- target_name_pair_index_1_pair_index_2
   *         ...
   */
  void write(hdf_archive& file) override;
  void collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators) override;

private:
  void resetTargetParticleSet(ParticleSet& P) override;

  /// generate the unique pair id from the group ids of particle i and j and the number of species
  static int gen_pair_id(const int ig, const int jg, const int ns);

  void report();

  void set_norm_factor();

private:
  PairCorrelationInput& input_;

  // While these seem like they should be the same as in input there
  // is a bunch ways these can get updated.
  Real rmax_;
  Real num_bins_;
  Real delta_;

  Real delta_inv_;
  /// volume of the cell
  Real volume_;
  ///save pair indices
  std::vector<int> pair_ids_;
  /// table indexs for other type
  std::vector<int> other_ids_;
  /// offset of the gofr's associated with others_id
  std::vector<int> other_offsets_;
  /////save source indices
  //vector<int> source_ids;
  ///normalization factor
  Matrix<Real> norm_factor_;
  int num_species_;
  int n_e_;
  std::vector<Real> n_vec_;
  // AA table ID
  int d_aa_id_;
  /////data
  //Matrix<RealType> gof_r;
  ///prefix of each gof_r
  std::vector<std::string> gof_r_prefix_;
  /** resize the internal data
   * @param nbins number of bins for the historgram
   */
  void resize(int nbins);
};

} // namespace qmcplusplus
