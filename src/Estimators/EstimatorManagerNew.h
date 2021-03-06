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
#include "Utilities/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/EstimatorManagerInterface.h"
#include "OperatorEstBase.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"
#include "type_traits/template_types.hpp"
#include <bitset>

namespace qmcplusplus
{
class QMCHamiltonian;
class CollectablesEstimator;

namespace testing
{
class EstimatorManagerNewTest;
} // namespace testing


/** Class to manage a set of ScalarEstimators */
class EstimatorManagerNew
{
public:
  /// This is to deal with vague expression of precision in legacy code. Don't use in new code.
  typedef QMCTraits::FullPrecRealType RealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;

  using QMCT = QMCTraits;
  typedef ScalarEstimatorBase EstimatorType;
  using FPRBuffer = std::vector<FullPrecRealType>;
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  ///default constructor
  EstimatorManagerNew(Communicate* c = 0);
  ///copy constructor, deleted
  EstimatorManagerNew(EstimatorManagerNew& em) = delete;
  ///destructor
  ~EstimatorManagerNew();

  ///return the number of ScalarEstimators
  inline int size() const { return Estimators.size(); }

  /** add a property with a name
   * @param aname name of the column
   * @return the property index so that its value can be set by setProperty(i)
   *
   * Append a named column. BlockProperties do not contain any meaning data
   * but manages the name to index map for PropertyCache.
   */
  inline int addProperty(const char* aname) { return BlockProperties.add(aname); }

  /** set the value of the i-th column with a value v
   * @param i column index
   * @param v value
   */
  inline void setProperty(int i, RealType v) { PropertyCache[i] = v; }

  inline RealType getProperty(int i) const { return PropertyCache[i]; }

  int addObservable(const char* aname);

  inline RealType getObservable(int i) const { return TotalAverages[i]; }

  void getData(int i, std::vector<RealType>& values);

  /** add an Estimator
   * @param newestimator New Estimator
   * @param aname name of the estimator
   * @return locator of newestimator
   */
  int add(EstimatorType* newestimator, const std::string& aname);

  /** add a main estimator
   * @param newestimator New Estimator
   * @return locator of newestimator
   */
  int add(EstimatorType* newestimator) { return add(newestimator, MainEstimatorName); }

  /** add a "non" physical operator estimator 
   *
   *  this is a dratically reduced version of OperatorBase right now it just supports
   *  what the SpinDensityNew estimator needs
   *
   *  What is actually important is that it has its own locality aware data and
   *  EstimatorManagerNew doesn't know about or manage that data.
   */
  int addEstOperator(OperatorEstBase& op_est);

  ///return a pointer to the estimator aname
  EstimatorType* getEstimator(const std::string& a);

  ///return the average for estimator i
  inline RealType average(int i) const { return Estimators[i]->average(); }

  ///returns a variance for estimator i
  inline RealType variance(int i) const { return Estimators[i]->variance(); }

  ///process xml tag associated with estimators
  bool put(QMCHamiltonian& H, const ParticleSet& pset, xmlNodePtr cur);

  /** reset the estimator
   */
  void reset();

  /** start a run
   * @param blocks number of blocks
   * @param record if true, will write to a file
   *
   * Replace reportHeader and reset functon.
   */
  void start(int blocks, bool record = true);

  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  /** unified: stop a block
   * @param accept acceptance rate of this block
   * \param[in] accept
   * \param[in] reject
   * \param[in] block_weight
   * \param[in] cpu_block_time Timer returns double so this is not altered by "mixed" precision
   */
  void stopBlock(unsigned long accept, unsigned long reject, RealType block_weight, double cpu_block_time);

  /** At end of block collect the scalar estimators for the entire rank
   *   
   *  \todo remove assumption of one ScalarEstimator per crowd.
   *  see how OperatorEstimators are handled
   *
   *  Each is currently accumulates on for crowd of 1 or more walkers
   *  returns the total weight across all crowds. 
   */
  RealType collectScalarEstimators(const RefVector<ScalarEstimatorBase>& scalar_estimators);

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

  template<class CT>
  void write(CT& anything, bool doappend)
  {
    anything.write(h_file, doappend);
  }

  auto& get_AverageCache() { return AverageCache; }
  auto& get_SquaredAverageCache() { return SquaredAverageCache; }

private:
  friend class EstimatorManagerCrowd;
  ///name of the primary estimator name
  std::string MainEstimatorName;
  //  TODO: fix needless use of bitset instead of clearer more visible booleans
  std::bitset<8> Options;
  ///size of the message buffer
  int BufferSize;
  ///number of records in a block
  int RecordCount;
  ///index for the block weight PropertyCache(weightInd)
  int weightInd;
  ///index for the block cpu PropertyCache(cpuInd)
  int cpuInd;
  ///index for the accept counter PropertyCache(acceptInd)
  int acceptRatioInd;
  ///hdf5 handler
  hid_t h_file;
  ///total weight accumulated in a block
  RealType BlockWeight;
  ///file handler to write data
  std::ofstream* Archive;
  ///file handler to write data for debugging
  std::ofstream* DebugArchive;
  ///communicator to handle communication
  Communicate* my_comm_;
  /** pointer to the primary ScalarEstimatorBase
   */
  ScalarEstimatorBase* MainEstimator;
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
  ///cached block averages of the squared values
  Vector<RealType> SquaredAverageCache;
  ///cached block averages of properties, e.g. BlockCPU
  Vector<RealType> PropertyCache;
  ///manager of scalar data
  RecordNamedProperty<RealType> BlockAverages;
  ///manager of property data
  RecordNamedProperty<RealType> BlockProperties;
  ///block averages: name to value
  RecordNamedProperty<RealType> TotalAverages;
  ///data accumulated over the blocks
  Matrix<RealType> TotalAveragesData;
  ///index mapping between BlockAverages and TotalAverages
  std::vector<int> Block2Total;
  ///column map
  std::map<std::string, int> EstimatorMap;
  ///estimators of simple scalars
  std::vector<EstimatorType*> Estimators;
  ///convenient descriptors for hdf5
  std::vector<observable_helper*> h5desc;
  /** OperatorEst Observables
   *
   * since the operator estimators are also a close set at compile time
   * they could be treated just like the inputs.
   * However the idea of a shared interface is much more straight forward for
   * them.
   */
  std::vector<std::unique_ptr<OperatorEstBase>> operator_ests_;

  /////estimators of composite data
  //CompositeEstimatorSet* CompEstimators;
  ///Timer
  Timer MyTimer;

private:
  ///number of maximum data for a scalar.dat
  int max4ascii;

  /// collect data and write
  void makeBlockAverages(unsigned long accept, unsigned long reject);

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

  ///add header to an std::ostream
  void addHeader(std::ostream& o);
  size_t FieldWidth;

  friend class qmcplusplus::testing::EstimatorManagerNewTest;
};
} // namespace qmcplusplus
#endif
