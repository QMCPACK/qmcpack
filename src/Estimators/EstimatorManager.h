//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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
    
    


/** @file EstimatorManager.h
 * @brief Manager class of scalar estimators
 */
#ifndef QMCPLUSPLUS_ESTIMATORMANAGER_H
#define QMCPLUSPLUS_ESTIMATORMANAGER_H

#include "Configuration.h"
#include "Utilities/Timer.h"
#include "Utilities/PooledData.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"
#include <bitset>

namespace qmcplusplus
{

class MCWalkerConifugration;
class QMCHamiltonian;
class CollectablesEstimator;

/**Class to manage a set of ScalarEstimators */
class EstimatorManager
{

public:

  typedef QMCTraits::EstimatorRealType  RealType;
  typedef ScalarEstimatorBase           EstimatorType;
  typedef std::vector<RealType>              BufferType;
  //enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

  ///name of the primary estimator name
  std::string MainEstimatorName;
  ///the root file name
  std::string RootName;
  ///energy
  TinyVector<RealType,4> RefEnergy;
  // //Cummulative energy, weight and variance
  // TinyVector<RealType,4>  EPSum;
  ///default constructor
  EstimatorManager(Communicate* c=0);
  ///copy constructor
  EstimatorManager(EstimatorManager& em);
  ///destructor
  virtual ~EstimatorManager();

  /** set the communicator */
  void setCommunicator(Communicate* c);

  /** return the communicator
   */
  Communicate* getCommunicator()
  {
    return myComm;
  }

  /** return true if the rank == 0
   */
  inline bool is_manager() const
  {
    return !myComm->rank();
  }

  ///return the number of ScalarEstimators
  inline int size() const
  {
    return Estimators.size();
  }

  /** add a property with a name
   * @param aname name of the column
   * @return the property index so that its value can be set by setProperty(i)
   *
   * Append a named column. BlockProperties do not contain any meaning data
   * but manages the name to index map for PropertyCache.
   */
  inline int addProperty(const char* aname)
  {
    return BlockProperties.add(aname);
  }

  /** set the value of the i-th column with a value v
   * @param i column index
   * @param v value
   */
  inline void setProperty(int i, RealType v)
  {
    PropertyCache[i]=v;
  }

  inline RealType getProperty(int i) const
  {
    return PropertyCache[i];
  }

  int addObservable(const char* aname);

  inline RealType getObservable(int i) const
  {
    return  TotalAverages[i];
  }

  void getData(int i, std::vector<RealType>& values);

  /** add an Estimator
   * @param newestimator New Estimator
   * @param aname name of the estimator
   * @return locator of newestimator
   */
  int add(EstimatorType* newestimator, const std::string& aname);
  //int add(CompositeEstimatorBase* newestimator, const std::string& aname);

  /** add a main estimator
   * @param newestimator New Estimator
   * @return locator of newestimator
   */
  int add(EstimatorType* newestimator)
  {
    return add(newestimator,MainEstimatorName);
  }

  ///return a pointer to the estimator aname
  EstimatorType* getEstimator(const std::string& a);

  ///return a pointer to the estimator
  EstimatorType* getMainEstimator();

  ///return the average for estimator i
  inline RealType average(int i) const
  {
    return Estimators[i]->average();
  }

  ///returns a variance for estimator i
  inline RealType variance(int i) const
  {
    return Estimators[i]->variance();
  }

  void setCollectionMode(bool collect);
  //void setAccumulateMode (bool setAccum) {AccumulateBlocks = setAccum;};

  ///process xml tag associated with estimators
  //bool put(xmlNodePtr cur);
  bool put(MCWalkerConfiguration& W, QMCHamiltonian& H, xmlNodePtr cur);

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
  void start(int blocks, bool record=true);
  /** stop a qmc run
   *
   * Replace finalize();
   */
  void stop();
  /** stop a qmc run
   */
  void stop(const std::vector<EstimatorManager*> m);


  /** start  a block
   * @param steps number of steps in a block
   */
  void startBlock(int steps);

  void setNumberOfBlocks(int blocks)
  {
    for(int i=0; i<Estimators.size(); i++)
      Estimators[i]->setNumberOfBlocks(blocks);
  }

  /** stop a block
   * @param accept acceptance rate of this block
   */
  void stopBlock(RealType accept, bool collectall=true);
  /** stop a block
   * @param m list of estimator which has been collecting data independently
   */
  void stopBlock(const std::vector<EstimatorManager*>& m);

  /** accumulate the measurements
   * @param W walkers
   */
  void accumulate(MCWalkerConfiguration& W);

  /** accumulate the measurements for a subset of walkers [it,it_end)
   * @param W walkers
   * @param it first walker
   * @param it_end last walker
   */
  void accumulate(MCWalkerConfiguration& W, MCWalkerConfiguration::iterator it,
                  MCWalkerConfiguration::iterator it_end);

//     /** accumulate the FW observables
//      */
//     void accumulate(HDF5_FW_observables& OBS, HDF5_FW_weights& WGTS, std::vector<int>& Dims);

  ///** set the cummulative energy and weight
  void getEnergyAndWeight(RealType& e, RealType& w, RealType& var);

  void getCurrentStatistics(MCWalkerConfiguration& W, RealType& eavg, RealType& var);

  template<class CT>
  void write(CT& anything, bool doappend)
  {
    anything.write(h_file,doappend);
  }

protected:
  ///use bitset to handle options
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
  hid_t h_file;
  ///total weight accumulated in a block
  RealType BlockWeight;
  ///file handler to write data
  std::ofstream* Archive;
  ///file handler to write data for debugging
  std::ofstream* DebugArchive;
  ///communicator to handle communication
  Communicate* myComm;
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
  std::map<std::string,int> EstimatorMap;
  ///estimators of simple scalars
  std::vector<EstimatorType*> Estimators;
  ///convenient descriptors for hdf5
  std::vector<observable_helper*> h5desc;
  /////estimators of composite data
  //CompositeEstimatorSet* CompEstimators;
  ///Timer
  Timer MyTimer;
private:
  ///number of maximum data for a scalar.dat
  int max4ascii;
  ///number of requests
  int pendingRequests;
  //Data for communication
  std::vector<BufferType*> RemoteData;
  //storage for MPI_Request
  std::vector<Communicate::request> myRequest;
  ///collect data and write
  void collectBlockAverages(int num_threads);
  ///add header to an std::ostream
  void addHeader(std::ostream& o);
  size_t FieldWidth;
};
}
#endif
