//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file EstimatorManagerBase.h
 * @brief Manager class of scalar estimators
 */
#ifndef QMCPLUSPLUS_ESTIMATORMANAGERBASE_H
#define QMCPLUSPLUS_ESTIMATORMANAGERBASE_H

#include "Configuration.h"
#include "Utilities/Timer.h"
#include "Pools/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "io/hdf/hdf_archive.h"
#include <bitset>

namespace qmcplusplus
{
class MCWalkerConfiguration;
class QMCHamiltonian;
class CollectablesEstimator;

namespace testing
{
class EstimatorManagerBaseTest;
} // namespace testing


/** Class to manage a set of ScalarEstimators */
class EstimatorManagerBase
{
public:
  using RealType         = QMCTraits::FullPrecRealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;

  using EstimatorType = ScalarEstimatorBase;
  using BufferType    = std::vector<RealType>;
  using MCPWalker     = Walker<QMCTraits, PtclOnLatticeTraits>;

  ///default constructor
  EstimatorManagerBase(Communicate* c = 0);
  ///copy constructor
  EstimatorManagerBase(EstimatorManagerBase& em);
  ///destructor
  virtual ~EstimatorManagerBase();

  /** set the communicator */
  void setCommunicator(Communicate* c);

  /** return the communicator
   */
  Communicate* getCommunicator() { return myComm; }

  /** return true if the rank == 0
   */
  inline bool is_manager() const { return !myComm->rank(); }

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
  int add(std::unique_ptr<EstimatorType> newestimator, const std::string& aname);
  //int add(CompositeEstimatorBase* newestimator, const std::string& aname);

  /** add a main estimator
   * @param newestimator New Estimator
   * @return locator of newestimator
   */
  int add(std::unique_ptr<EstimatorType> newestimator) { return add(std::move(newestimator), main_estimator_name_); }

  ///return a pointer to the estimator aname
  EstimatorType* getEstimator(const std::string& a);

  ///return a pointer to the estimator
  EstimatorType* getMainEstimator();

  void setCollectionMode(bool collect);
  //void setAccumulateMode (bool setAccum) {AccumulateBlocks = setAccum;};

  ///process xml tag associated with estimators
  //bool put(xmlNodePtr cur);
  bool put(QMCHamiltonian& H, xmlNodePtr cur);

  void resetTargetParticleSet(ParticleSet& p);

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
  /** stop a qmc run
   *
   * Replace finalize();
   */
  void stop();
  /** stop a qmc run
   */
  void stop(const std::vector<EstimatorManagerBase*> m);


  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  /** stop a block
   * @param accept acceptance rate of this block
   */
  void stopBlock(RealType accept, bool collectall = true);

  /** stop a block
   * @param m list of estimator which has been collecting data independently
   */
  void stopBlock(const std::vector<EstimatorManagerBase*>& m);

  /** accumulate the measurements
   * @param W walkers
   */
  void accumulate(MCWalkerConfiguration& W);

  /** accumulate the measurements for a subset of walkers [it,it_end)
   * @param W walkers
   * @param it first walker
   * @param it_end last walker
   */
  void accumulate(MCWalkerConfiguration& W, MCWalkerConfiguration::iterator it, MCWalkerConfiguration::iterator it_end);

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

protected:
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
  ///index for the acceptance rate PropertyCache(acceptInd)
  int acceptInd;
  ///hdf5 handler
  hdf_archive h_file;
  ///total weight accumulated in a block
  RealType BlockWeight;
  ///file handler to write data
  std::unique_ptr<std::ofstream> Archive;
#if defined(DEBUG_ESTIMATOR_ARCHIVE)
  ///file handler to write data for debugging
  std::unique_ptr<std::ofstream> DebugArchive;
#endif
  ///communicator to handle communication
  Communicate* myComm;
  /** pointer to the primary ScalarEstimatorBase
   */
  ScalarEstimatorBase* MainEstimator;
  /** pointer to the CollectablesEstimator
   *
   * Do not need to clone: owned by the master thread
   */
  std::unique_ptr<CollectablesEstimator> Collectables;
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
  UPtrVector<EstimatorType> Estimators;
  ///convenient descriptors for hdf5
  std::vector<ObservableHelper> h5desc;
  /////estimators of composite data
  //CompositeEstimatorSet* CompEstimators;
  ///Timer
  Timer MyTimer;

private:
  ///name of the primary estimator name
  std::string main_estimator_name_;

  ///number of maximum data for a scalar.dat
  int max_output_scalar_dat_;

  //Data for communication
  std::vector<std::unique_ptr<BufferType>> RemoteData;

  ///collect data and write
  void collectBlockAverages();

  ///add header to an std::ostream
  void addHeader(std::ostream& o);

  ///largest name in BlockAverages adding 2 characters
  size_t max_block_avg_name_;

  friend class qmcplusplus::testing::EstimatorManagerBaseTest;
};
} // namespace qmcplusplus
#endif
