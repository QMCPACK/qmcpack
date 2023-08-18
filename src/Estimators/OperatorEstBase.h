//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorBase.h
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 */
#ifndef QMCPLUSPLUS_OPERATORESTBASE_H
#define QMCPLUSPLUS_OPERATORESTBASE_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "type_traits/DataLocality.h"
#include "hdf/hdf_archive.h"
#include <bitset>

namespace qmcplusplus
{
class TrialWaveFunction;
namespace testing
{
class OEBAccessor;
}
/** @ingroup Estimators
 * @brief An abstract class for gridded estimators
 *
 */
class OperatorEstBase
{
public:
  using QMCT      = QMCTraits;
  using FullPrecRealType = QMCT::FullPrecRealType;
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  using Data = std::vector<QMCT::RealType>;

  ///constructor
  OperatorEstBase(DataLocality dl);
  /** Shallow copy constructor!
   *  This alows us to keep the default copy constructors for derived classes which
   *  is quite useful to the spawnCrowdClone design.
   *  Data is likely to be quite large and since the OperatorEstBase design is that the children 
   *  reduce to the parent it is infact undesirable for them to copy the data the parent has.
   *  Initialization of Data (i.e. call to resize) if any is the responsibility of the derived class.
   */
  OperatorEstBase(const OperatorEstBase& oth);
  ///virtual destructor
  virtual ~OperatorEstBase() = default;

  /** Accumulate whatever it is you are accumulating with respect to walkers
   * 
   *  This method is assumed to be called from the crowd context
   *  It provides parallelism with respect to computational effort of the estimator
   *  without causing a global sync.
   *  Depending on data locality the accumlation of the result may be different from
   *  the single thread write directly into the OperatorEstimator data.
   *  \param[in]      walkers
   *  \param[inout]   pset_target   crowd scope target pset (should be returned to starting state after call)
   *  \param[in]      psets         per walker psets
   *  \param[in]      wnfs          per walker TrialWaveFunction
   *  \param[inout]   rng           crowd scope RandomGenerator
   */
  virtual void accumulate(const RefVector<MCPWalker>& walkers,
                          const RefVector<ParticleSet>& psets,
                          const RefVector<TrialWaveFunction>& wfns,
                          RandomBase<FullPrecRealType>& rng) = 0;

  /** Reduce estimator result data from crowds to rank
   *
   *  This is assumed to be called from only from one thread per crowds->rank
   *  reduction. Implied is this is during a global sync or there is a guarantee
   *  that the crowd operator estimators accumulation data is not being written to.
   *
   *  There could be concurrent operations inside the scope of the collect call.
   */
  virtual void collect(const RefVector<OperatorEstBase>& oebs);

  virtual void normalize(QMCT::RealType invToWgt);

  virtual void startBlock(int steps) = 0;

  std::vector<QMCT::RealType>& get_data() { return data_; }

  /*** create and tie OperatorEstimator's observable_helper hdf5 wrapper to stat.h5 file
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerOperatorEstimator(hdf_archive& file) {}

  virtual std::unique_ptr<OperatorEstBase> spawnCrowdClone() const = 0;

  /** Write to previously registered observable_helper hdf5 wrapper.
   *
   *  if you haven't registered Operator Estimator 
   *  this will do nothing.
   */
  void write(hdf_archive& file);

  /** zero data appropriately for the DataLocality
   */
  void zero();

  /** Return the total walker weight for this block
   */
  QMCT::FullPrecRealType get_walkers_weight() const { return walkers_weight_; }

  const std::string& get_my_name() const { return my_name_; }

  /** Register 0-many listeners with a leading QMCHamiltonian instance i.e. a QMCHamiltonian
   *  that has acquired the crowd scope QMCHamiltonianMultiWalkerResource.
   *  This must be called for each crowd scope estimator that listens to register listeners into
   *  the crowd scope QMCHamiltonianMultiWalkerResource.
   *
   *  Many estimators don't need per particle values so the default implementation is no op.
   */
  virtual void registerListeners(QMCHamiltonian& ham_leader){};

  bool isListenerRequired() { return requires_listener_; }

  DataLocality get_data_locality() const { return data_locality_; }

protected:
  /** locality for accumulation of estimator data.
   *  This designates the memory scheme used for the estimator
   *  The default is:
   *  DataLocality::Crowd, each crowd and the rank level estimator have a full representation of the data
   *  Memory Savings Schemes:
   *  One:
   *  DataLocality::Rank,  This estimator has the full representation of the data but its crowd spawn will have
   *  One per crowd:
   *  DataLocality::Queue  This estimator accumulates queue of values to collect to the Rank estimator data
   *  DataLocality::?      Another way to reduce memory use on thread/crowd local estimators.
   */
  DataLocality data_locality_;

  ///name of this object -- only used for debugging and h5 output
  std::string my_name_;

  QMCT::FullPrecRealType walkers_weight_;

  // convenient Descriptors hdf5 for Operator Estimators only populated for rank scope OperatorEstimator
  std::vector<ObservableHelper> h5desc_;

  Data data_;

  bool requires_listener_ = false;

  friend testing::OEBAccessor;
};
} // namespace qmcplusplus
#endif
