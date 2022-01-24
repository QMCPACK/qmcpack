//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: EstimatorManagerBase.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATORMANAGERNEW_H
#define QMCPLUSPLUS_ESTIMATORMANAGERNEW_H

#include <memory>

#include "Configuration.h"
#include "Utilities/Timer.h"
#include "Pools/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "OperatorEstBase.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"
#include "type_traits/template_types.hpp"
#include <bitset>

namespace qmcplusplus
{
class QMCHamiltonian;
class WaveFunctionFactory;
class CollectablesEstimator;
class hdf_archive;

namespace testing
{
class EstimatorManagerNewTest;
} // namespace testing


/** Class to manage a set of ScalarEstimators
 * As a manager, this class handles the aggregation of data from crowds, MPI ranks and I/O logics.
 * The actually per-crowd data accumulation is done by EstimatorManagerCrowd.
 */
class EstimatorManagerNew
{
public:
  /// This is to deal with vague expression of precision in legacy code. Don't use in new code.
  using RealType         = QMCTraits::FullPrecRealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;

  using QMCT          = QMCTraits;
  using EstimatorType = ScalarEstimatorBase;
  using FPRBuffer     = std::vector<FullPrecRealType>;
  using MCPWalker     = Walker<QMCTraits, PtclOnLatticeTraits>;

  ///default constructor
  EstimatorManagerNew(Communicate* c);
  ///copy constructor, deleted
  EstimatorManagerNew(EstimatorManagerNew& em) = delete;
  ///destructor
  ~EstimatorManagerNew();

  ///return the number of ScalarEstimators
  inline int size() const { return Estimators.size(); }

  /** add a "non" physical operator estimator 
   *
   *  this is a dratically reduced version of OperatorBase right now it just supports
   *  what the SpinDensityNew estimator needs
   *
   *  What is actually important is that it has its own locality aware data and
   *  EstimatorManagerNew doesn't know about or manage that data.
   */
  int addEstOperator(OperatorEstBase& op_est);

  ///process xml tag associated with estimators
  bool put(QMCHamiltonian& H,
           const ParticleSet& pset,
           const TrialWaveFunction& twf,
           const WaveFunctionFactory& wf_factory,
           xmlNodePtr cur);

  /** Start the manager at the beginning of a driver run().
   * Open files. Setting zeros.
   * @param blocks number of blocks
   * @param record if true, will write to a file
   *
   * Replace reportHeader and reset functon.
   */
  void startDriverRun();

  /** Stop the manager at the end of a driver run().
   * Flush/close files.
   */
  void stopDriverRun();

  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  /** unified: stop a block
   * @param accept acceptance rate of this block
   * \param[in] accept
   * \param[in] reject
   * \param[in] block_weight
   */
  void stopBlock(unsigned long accept, unsigned long reject, RealType block_weight);

  /** At end of block collect the scalar estimators for the entire rank
   *   
   *  \todo remove assumption of one ScalarEstimator per crowd.
   *  see how OperatorEstimators are handled
   *
   *  Each is currently accumulates on for crowd of 1 or more walkers
   */
  void collectScalarEstimators(const RefVector<ScalarEstimatorBase>& scalar_estimators);

  /** Reduces OperatorEstimator data from Crowds to the manager's OperatorEstimator data
   *
   *  \param[in] op_ests - vector of each crowds vector of references to their OperatorEstimators.
   *
   *  A particular OperatorEstimators reduction via a call to collect may be straight forward
   *  if the crowd context OperatorEstimator holds a copy of the estimator data structure
   *  or more complex if it just collects for instance a list of writes to locations
   *  in the data structure.
   */
  void collectOperatorEstimators(const std::vector<RefVector<OperatorEstBase>>& op_ests);

  /** get the average of per-block energy and variance of all the blocks
   * Note: this is not weighted average. It can be the same as weighted average only when block weights are identical.
   */
  void getApproximateEnergyVariance(RealType& e, RealType& var);

  auto& get_AverageCache() { return AverageCache; }

private:
  /** reset the estimator
   */
  void reset();

  /** add an Estimator
   * @param newestimator New Estimator
   * @param aname name of the estimator
   * @return locator of newestimator
   */
  int add(std::unique_ptr<EstimatorType> newestimator, const std::string& aname);

  ///return a pointer to the estimator aname
  EstimatorType* getEstimator(const std::string& a);

  /// collect data and write
  void makeBlockAverages(unsigned long accept, unsigned long reject);

  /// write scalars to scalar.dat and h5
  void writeScalarH5();

  /** do the rank wise reduction of the OperatorEstimators
   *
   *  Why do this here?
   *  1. Operator estimators don't know about the concurrency model
   *  2. EstimatorManager owns the resources:
   *       send & receive buffers
   *  3. The operation is generic as long as OperatorEstimator satisfies
   *     the requirement that get_data_ref() returns a reference to
   *     std::vector<RealType>
   *
   *  Implementation makes the assumption that sending each OperatorEstimator
   *  separately is the correct memory use vs. mpi message balance.
   */
  void reduceOperatorEstimators();
  /** Write OperatorEstimator data to *.stat.h5
   *
   *  Note that OperatorEstimator owns its own observable_helpers
   */
  void writeOperatorEstimators();
  /** OperatorEstimators need to be zeroed out after the block is finished.
   */
  void zeroOperatorEstimators();

  ///name of the primary estimator name
  std::string MainEstimatorName;
  ///number of records in a block
  int RecordCount;
  ///index for the block weight PropertyCache(weightInd)
  int weightInd;
  ///index for the block cpu PropertyCache(cpuInd)
  int cpuInd;
  ///index for the accept counter PropertyCache(acceptInd)
  int acceptRatioInd;
  ///hdf5 handler
  std::unique_ptr<hdf_archive> h_file;
  ///file handler to write data
  std::unique_ptr<std::ofstream> Archive;
  ///file handler to write data for debugging
  std::unique_ptr<std::ofstream> DebugArchive;
  ///communicator to handle communication
  Communicate* my_comm_;
  /** pointer to the CollectablesEstimator
   *
   * Do not need to clone: owned by the master thread
   */
  CollectablesEstimator* Collectables;
  /** accumulator for the energy
   *
   * @todo expand it for all the scalar observables to report the final results
   */
  ScalarEstimatorBase::accumulator_type energyAccumulator;
  /** accumulator for the variance **/
  ScalarEstimatorBase::accumulator_type varAccumulator;
  ///cached block averages of the values
  Vector<RealType> AverageCache;
  ///cached block averages of properties, e.g. BlockCPU
  Vector<RealType> PropertyCache;
  ///manager of scalar data
  RecordNamedProperty<RealType> BlockAverages;
  ///manager of property data
  RecordNamedProperty<RealType> BlockProperties;
  ///column map
  std::map<std::string, int> EstimatorMap;
  ///estimators of simple scalars
  std::vector<std::unique_ptr<EstimatorType>> Estimators;
  ///convenient descriptors for hdf5
  std::vector<ObservableHelper> h5desc;
  /** OperatorEst Observables
   *
   * since the operator estimators are also a close set at compile time
   * they could be treated just like the inputs.
   * However the idea of a shared interface is much more straight forward for
   * them.
   */
  std::vector<std::unique_ptr<OperatorEstBase>> operator_ests_;

  ///block timer
  Timer block_timer_;

  ///number of maximum data for a scalar.dat
  int max4ascii;

  ///add header to an std::ostream
  void addHeader(std::ostream& o);
  size_t FieldWidth;

  friend class EstimatorManagerCrowd;
  friend class qmcplusplus::testing::EstimatorManagerNewTest;
};
} // namespace qmcplusplus
#endif
