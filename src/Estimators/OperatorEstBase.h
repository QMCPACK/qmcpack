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
  using Data = std::variant<std::unique_ptr<std::vector<QMCT::RealType>>, std::shared_ptr<std::vector<QMCT::RealType>>>;

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
  virtual ~OperatorEstBase() {}

  /** Accumulate whatever it is you are accumulating with respect to walkers
   * 
   *  This method is assumed to be called from a concurrent context
   *  This method is responsible for being thread safe.
   */
  virtual void accumulate(RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets) = 0;

  /** This is assumed to be called from only one thread with respect to oeb
   *
   *  oeb's Data will be written to in a non thread safe manner.
   */
  virtual void collect(const OperatorEstBase& oeb) = 0;

  std::vector<QMCT::RealType>& get_data_ref() { return std::visit([](auto& data) -> std::vector<QMCT::RealType>& {
        return *data; }, data_); }
  
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
