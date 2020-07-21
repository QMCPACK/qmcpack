//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
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
#include "Utilities/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/EstimatorManagerInterface.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{
class MCWalkerConifugration;
class QMCHamiltonian;
class CollectablesEstimator;

/** Thread local portion of EstimatorManager
 *
 *  Stepping away from the CloneManger + clones design which creates EstimatorManagers
 *  Which operate differently based on internal switches.
 */
class EstimatorManagerCrowd : public EstimatorManagerInterface
{
public:
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
  //enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

  ///name of the primary estimator name
  std::string MainEstimatorName;
  ///the root file name
  std::string RootName;
  ///energy
  TinyVector<RealType, 4> RefEnergy;
  ///default constructor
  EstimatorManagerCrowd() = delete;
  ///copy constructor
  EstimatorManagerCrowd(EstimatorManagerNew& em);
  ///destructor
  virtual ~EstimatorManagerCrowd(){};

  inline bool is_manager() const { return false; }

  ///return the number of ScalarEstimators
  inline int size() const { return scalar_estimators_.size(); }

  /** start a run
   * @param blocks number of blocks
   * @param record if true, will write to a file
   *
   * Replace reportHeader and reset functon.
   */
  void start(int blocks, bool record = true) {}

  /** stop a qmc run
   *
   * Replace finalize();
   */
  void stop() {}

  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  void stopBlock();

  void setNumberOfBlocks(int blocks)
  {
    for (int i = 0; i < scalar_estimators_.size(); i++)
      scalar_estimators_[i]->setNumberOfBlocks(blocks);
  }

  void accumulate(int global_walkers, RefVector<MCPWalker>& walkers, RefVector<ParticleSet>& psets)
  {
    block_weight_ += walkers.size();
    //Don't normalize we only divide once after reduction.
    //RealType norm             = 1.0 / global_walkers;
    int num_scalar_estimators = scalar_estimators_.size();
    for (int i = 0; i < num_scalar_estimators; ++i)
      scalar_estimators_[i]->accumulate(global_walkers, walkers, 1);
  }

  RefVector<EstimatorType> get_scalar_estimators() { return convertPtrToRefVector(scalar_estimators_); }
  RealType get_block_weight() const { return block_weight_; }
  double get_cpu_block_time() const { return cpu_block_time_; }

protected:
  ///use bitset to handle options
  std::bitset<8> Options;
  ///number of records in a block
  int RecordCount;
  ///index for the block weight PropertyCache(weightInd)
  int weightInd;
  ///index for the block cpu PropertyCache(cpuInd)
  int cpuInd;
  ///index for the acceptance rate PropertyCache(acceptInd)
  int acceptInd;
  ///total weight accumulated in a block
  RealType block_weight_;

  ///file handler to write data
  std::ofstream* Archive;
  ///file handler to write data for debugging
  std::ofstream* DebugArchive;
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
  std::vector<EstimatorType*> scalar_estimators_;

private:
  // This is needed for "efficiency" measure
  Timer crowd_estimator_timer_;
  double cpu_block_time_ = 0;

  ///number of maximum data for a scalar.dat
  int max4ascii;
  ///collect data and write
  void collectBlockAverages(int num_threads);
  ///add header to an std::ostream
  void addHeader(std::ostream& o);
  size_t FieldWidth;
};

} // namespace qmcplusplus

#endif
