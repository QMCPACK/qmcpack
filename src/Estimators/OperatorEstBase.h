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

#include <variant>

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/observable_helper.h"
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

  /// locality for accumulation data
  DataLocality data_locality_;

  ///name of this object
  std::string myName;

  QMCT::FullPrecRealType walkers_weight_;

  QMCT::FullPrecRealType get_walkers_weight() const { return walkers_weight_; }
  ///constructor
  OperatorEstBase();
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
   */
  virtual void accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets) = 0;

  /** Reduce estimator result data from crowds to rank
   *
   *  This is assumed to be called from only from one thread per crowds->rank
   *  reduction. Implied is this is during a global sync or there is a guarantee
   *  that the crowd operator estimators accumulation data is not being written to.
   *
   *  There could be concurrent operations inside the scope of the collect call.
   */
  virtual void collect(const RefVector<OperatorEstBase>& oebs);

  std::vector<QMCT::RealType>& get_data_ref() { return *data_; }

  Data& get_data() { return data_; };
  /*** add to OperatorEstimator descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerOperatorEstimator(std::vector<observable_helper*>& h5desc, hid_t gid) const {}

  virtual OperatorEstBase* clone() = 0;

protected:
  /** data management
   */
  static Data createLocalData(size_t size, DataLocality data_locality);

  Data data_;
};
} // namespace qmcplusplus
#endif
