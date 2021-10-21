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
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "type_traits/DataLocality.h"
#include <bitset>

namespace qmcplusplus
{
class DistanceTableData;
class TrialWaveFunction;

/** @ingroup Estimators
 * @brief An abstract class for gridded estimators
 *
 */
class OperatorEstBase
{
public:
  using QMCT      = QMCTraits;
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  /** the type in this variant changes based on data locality
   */
  using Data = UPtr<std::vector<QMCT::RealType>>;

  /// locality for accumulation data. FIXME full documentation of this state machine.
  DataLocality data_locality_;

  ///name of this object
  std::string myName;

  QMCT::FullPrecRealType get_walkers_weight() const { return walkers_weight_; }
  ///constructor
  OperatorEstBase(DataLocality dl);
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
                          RandomGenerator_t& rng) = 0;

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

  std::vector<QMCT::RealType>& get_data_ref() { return *data_; }

  Data& get_data() { return data_; };

  /*** create and tie OperatorEstimator's observable_helper hdf5 wrapper to stat.h5 file
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerOperatorEstimator(hid_t gid) {}


  virtual OperatorEstBase* clone() = 0;

  /** Write to previously registered observable_helper hdf5 wrapper.
   *
   *  if you haven't registered Operator Estimator 
   *  this will do nothing.
   */
  void write();

  /** zero data appropriately for the DataLocality
   */
  void zero();

  /** Return the total walker weight for this block
   */
  QMCT::FullPrecRealType get_walkers_weight() { return walkers_weight_; }

protected:
  QMCT::FullPrecRealType walkers_weight_;

  // convenient Descriptors hdf5 for Operator Estimators only populated for rank scope OperatorEstimator
  UPtrVector<ObservableHelper> h5desc_;

  /** create the typed data block for the Operator.
   *
   *  this is only slightly better than a byte buffer
   *  it allows easy porting of the legacy implementations
   *  Which wrote into a shared buffer per walker.
   *  And it make's datalocality fairly easy but
   *  more descriptive and safe data structures would be better
   */
  static Data createLocalData(size_t size, DataLocality data_locality);

  Data data_;
};
} // namespace qmcplusplus
#endif
