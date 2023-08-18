//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerBase.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATORMANAGERCROWD_H
#define QMCPLUSPLUS_ESTIMATORMANAGERCROWD_H

#include <bitset>

#include "Configuration.h"
#include "Utilities/Timer.h"
#include "Pools/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{
class MCWalkerConifugration;
class QMCHamiltonian;

/** Thread local estimator container/accumulator
 *
 *  Stepping away from the CloneManger + clones design which creates EstimatorManagers
 *  Which operate differently based on internal switches.
 *  
 *  see EstimatorManagerNew.h for full description of the new design.
 */
class EstimatorManagerCrowd
{
public:
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType  = EstimatorManagerNew::RealType;
  using FullPrecRealType = EstimatorManagerNew::FullPrecRealType;

  /** EstimatorManagerCrowd are always spawn of an EstimatorManagerNew
   */
  EstimatorManagerCrowd(EstimatorManagerNew& em);

  ///destructor
  ~EstimatorManagerCrowd(){};

  ///return the number of ScalarEstimators
  inline int size() const { return scalar_estimators_.size(); }

  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  void stopBlock();

  /** Accumulate over all scalar estimators and operator estimators over all walkers in crowd.
   *  Not all estimators make use of all these arguments
   *  \param[in]     walkers         walkers in crowd
   *  \param[in]     psets           walker particle sets
   *  \param[in]     wfns            walker wavefunctions
   *  \param[inout]  rng             crowd scope RandomGenerator
   *
   *  walkers is especially questionable since its really just hiding the full sweep hamiltonian values from
   *  the most recent (maybe) QMCHamiltonian evaluate which are written into it by the QMCHamiltonians
   *  previous to the accumulate.
   *  walkers might additionally be useful because they hold another copy of the dynamic (electron) particle sets coords
   *  that could be inconsistent with psets.
   *
   *  As soon as the legacy Estimators are dropped this API should be reviewed with an eye to disentangling
   *  ParticleSet, Walker, and QMCHamiltonian.
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecRealType>& rng);

  ScalarEstimatorBase& get_main_estimator() { return *main_estimator_; }
  RefVector<ScalarEstimatorBase> get_scalar_estimators() { return convertUPtrToRefVector(scalar_estimators_); }
  RefVector<qmcplusplus::OperatorEstBase> get_operator_estimators() { return convertUPtrToRefVector(operator_ests_); }

  RealType get_block_num_samples() const { return block_num_samples_; }
  RealType get_block_weight() const { return block_weight_; }

  /** This registers the crowd lever estimators that require listeners into the QMCHamiltonianMultiWalkerResources
   *  We really only need a QMCHamiltonian leader but resource acquisition and release works better this way.
   */
  void registerListeners(const RefVectorWithLeader<QMCHamiltonian>& ham_list);

private:
  ///number of samples accumulated in a block
  RealType block_num_samples_;
  ///total weight accumulated in a block
  RealType block_weight_;

  UPtr<ScalarEstimatorBase> main_estimator_;
  ///estimators of simple scalars
  UPtrVector<ScalarEstimatorBase> scalar_estimators_;

  UPtrVector<OperatorEstBase> operator_ests_;
};

} // namespace qmcplusplus

#endif
