//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#define QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/RecordProperty.h"
#include "Estimators/accumulators.h"
#include "Particle/Walker.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#endif

namespace qmcplusplus
{
class ObservableHelper;

/** Abstract class for an estimator of a scalar operator.
 *
 * ScalarEstimators derived from ScalarEstimatorBase  implement three main functions
 * - accumulate : measure and accumulate its value and the square of the value
 * - add2Record : \todo document this
 * - registerObservables : \todo document this
 * - clone : because all types must be erased
 * ScalarEstimatorBase and its derived classes do not perform any I/O function.
 */
struct ScalarEstimatorBase
{
  using RealType         = QMCTraits::FullPrecRealType;
  using accumulator_type = accumulator_set<RealType>;
  using Walker_t         = MCWalkerConfiguration::Walker_t;
  using MCPWalker        = Walker<QMCTraits, PtclOnLatticeTraits>;
  using WalkerIterator   = MCWalkerConfiguration::const_iterator;
  using RecordListType   = RecordNamedProperty<RealType>;

  ///first index within an record of the first element handled by an object
  int FirstIndex;
  ///last index within an record of the first element handled by an object
  int LastIndex;
  ///scalars to be measured
  std::vector<accumulator_type> scalars;
  ///scalars saved
  std::vector<accumulator_type> scalars_saved;
  //     RealType NSTEPS;

  inline ScalarEstimatorBase() : FirstIndex(0), LastIndex(0) {}

  virtual ~ScalarEstimatorBase() {}

  /// Is this estimator a main estimator i.e. the estimator required for a particular driver.
  virtual bool isMainEstimator() const { return false; }
  virtual std::string getName() const = 0;
  
  ///return average of the
  inline RealType average(int i = 0) const { return scalars_saved[i].mean(); }
  ///return a variance
  inline RealType variance(int i = 0) const { return scalars_saved[i].variance(); }
  ///retrun mean and variance
  inline std::pair<RealType, RealType> operator[](int i) const { return scalars[i].mean_and_variance(); }

  ///return the size of scalars it manages
  virtual inline int size() const { return scalars.size(); }

  ///clear the scalars to collect
  inline void clear()
  {
    for (int i = 0; i < scalars.size(); i++)
      scalars[i].clear();
  }

  /** take block average and write to a common container */
  template<typename IT>
  inline void takeBlockAverage(IT first)
  {
    first += FirstIndex;
    for (int i = 0; i < scalars.size(); i++)
    {
      *first++         = scalars[i].mean();
      scalars_saved[i] = scalars[i]; //save current block
      scalars[i].clear();
    }
  }

  /** take block average and write to common containers for values and squared values
   * @param first starting iterator of values
   * @param first_sq starting iterator of squared values
   */
  template<typename IT>
  inline void takeBlockAverage(IT first, IT first_sq)
  {
    first += FirstIndex;
    first_sq += FirstIndex;
    for (int i = 0; i < scalars.size(); i++)
    {
      *first++    = scalars[i].mean();
      *first_sq++ = scalars[i].mean2();
      // For mixed precision this is where data goes from Actual QMCT::RealType
      // to the local RealType which is hard coded to QMCTFullPrecRealType.
      // I think it is a bad idea to have RealType to have a changing meaning,
      // I also feel like having the floating point expension needed to write always to
      // double in hdf5 is poor form especially since I had to figure this out in many
      // years later.
      scalars_saved[i] = scalars[i]; //save current block
      scalars[i].clear();
    }
  }

  /** add the block accumulated scalars
   * @param first starting iterator of values
   */
  template<typename IT>
  inline void addAccumulated(IT first)
  {
    first += FirstIndex;
    for (int i = 0; i < scalars.size(); i++)
    {
      *first++ += scalars[i].result();
      scalars[i].clear();
    }
  }

  /** a virtual function to accumulate observables or collectables
   * @param W const MCWalkerConfiguration
   * @param first const_iterator for the first walker
   * @param last const_iterator for the last walker
   * @param wgt weight
   *
   * Pass W along with the iterators so that the properties of W can be utilized.
   */
  virtual void accumulate(const MCWalkerConfiguration& W, WalkerIterator first, WalkerIterator last, RealType wgt) = 0;

  /** a virtual function to accumulate observables or collectables
   * @param global_walkers_ walkers per ranks or walkers total?
   * @param RefVector of MCPWalkers
   * @param wgt weight or maybe norm
   *
   */
  virtual void accumulate(const RefVector<MCPWalker>&) = 0;

  /** add the content of the scalar estimator to the record
   * @param record scalar data list
   *
   * Each ScalarEstimatorBase object adds 1 to many accumulator_type
   */
  virtual void add2Record(RecordNamedProperty<RealType>& record) = 0;

  /** add descriptors of observables to utilize hdf5
   * @param h5desc descriptor of a data stored in a h5 group
   * @param file file to which each statistical data will be stored
   */
  virtual void registerObservables(std::vector<ObservableHelper>& h5dec, hdf_archive& file) = 0;

  ///clone the object
  virtual ScalarEstimatorBase* clone() = 0;

  /// String representation of the derived type of the ScalarEstimator
  virtual const std::string& getSubTypeStr() const = 0;
};
} // namespace qmcplusplus

#endif
