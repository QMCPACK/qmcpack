//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ScalarEstimatorManager.h
 * @brief Manager class of scalar estimators
 * @authors J. Kim and J. Vincent
 */
#ifndef QMCPLUSPLUS_SCALAR_ESTIMATORMANAGER_H
#define QMCPLUSPLUS_SCALAR_ESTIMATORMANAGER_H

#include "Configuration.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus {

  class MCWalkerConifugration;
  class QMCHamiltonian;

  /**Class to manage a set of ScalarEstimators */
  class ScalarEstimatorManager: public QMCTraits
  {

  public:

    typedef ScalarEstimatorBase           EstimatorType;
    typedef EstimatorType::BufferType     BufferType;
    //enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

    ///name of the primary estimator name
    std::string MainEstimatorName;
    ///the root file name
    std::string RootName;

    ScalarEstimatorManager(QMCHamiltonian& h, Communicate* c=0);
    ScalarEstimatorManager(ScalarEstimatorManager& em, QMCHamiltonian& h);
    virtual ~ScalarEstimatorManager();

    /** set the communicator */
    void setCommunicator(Communicate* c);

    /** return the communicator 
     */
    Communicate* getCommunicator()  
    { 
      return myComm;
    }

    ///return the number of ScalarEstimators
    inline int size() const { return Estimators.size();}

    /** add a column with the name
     * @param aname name of the column
     * @return the column index
     *
     * Append a named column. BlockProperties do not contain any meaning data
     * but manages the name to index map for PropertyCache.
     */
    inline int addColumn(const char* aname) {
      return BlockProperties.add(aname);
    }

    /** set the value of the i-th column with a value v
     * @param i column index
     * @param v value 
     */
    inline void setColumn(int i, RealType v) {
      PropertyCache(RecordCount,i)=v;
    }

    inline RealType getColumn(int i) const {
      return PropertyCache(RecordCount,i);
    }

    int addObservable(const char* aname);

    inline RealType getObservable(int i) const {
      return  TotalAverages[i];
    }

    void getData(int i, vector<RealType>& values);

    /** add an Estimator 
     * @param newestimator New Estimator
     * @param aname name of the estimator
     * @return locator of newestimator
     */
    int add(EstimatorType* newestimator, const string& aname);

    /** add a main estimator
     * @param newestimator New Estimator
     * @return locator of newestimator
     */
    int add(EstimatorType* newestimator)
    {
      return add(newestimator,MainEstimatorName);
    }

    ///return a pointer to the estimator aname
    EstimatorType* getEstimator(const string& a);

    ///return a pointer to the estimator 
    EstimatorType* getMainEstimator();
 
    ///return the average for estimator i
    inline RealType average(int i) const { 
      return Estimators[i]->average(); 
    }

    ///returns a variance for estimator i
    inline RealType variance(int i) const { 
      return Estimators[i]->variance();
    }

    void setCollectionMode(bool collect);
    //void setAccumulateMode (bool setAccum) {AccumulateBlocks = setAccum;};

    ///process xml tag associated with estimators
    bool put(xmlNodePtr cur);

    ///** reset the estimator
    // * @param aname root file name
    // * @param append if yes, the data should be appended.
    // *
    // * Replace resetReportSettings. 
    // */
    //void reset(const string& aname, bool append);
    
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
    void stop(const vector<ScalarEstimatorManager*> m);

    /** start  a block
     * @param steps number of steps in a block
     */
    inline void startBlock(int steps) { MyTimer.restart();}
    /** stop a block
     * @param accept acceptance rate of this block
     */
    void stopBlock(RealType accept);
    /** stop a block
     * @param m list of estimator which has been collecting data independently
     */
    void stopBlock(const vector<ScalarEstimatorManager*> m);

    void accumulate(MCWalkerConfiguration::iterator it,
        MCWalkerConfiguration::iterator it_end);

    /** accumulate the scalar observables
     */
    void accumulate(MCWalkerConfiguration& W)
    {
      accumulate(W.begin(),W.end());
    }

    /** accumulate the scalar observables
     */
    void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker);

    /** set the cummulative energy and weight
     */
    void getEnergyAndWeight(RealType& e, RealType& w);


  protected:
    ///if true, responsible for reduction operation, broadcast of EPSum, and printout text file
    bool Manager;
    ///if true, the averages are collected over mpi nodes
    bool CollectSum;
    ///if true, record is appened to an exisiting data
    bool AppendRecord;
    ///if true, record is collected. 
    bool Collected;
    ///size of multiple estimator
    int ThreadCount;
    ///number of records in a block
    int RecordCount;
    ///index for the block cpu PropertyCache(*,cpuInd) 
    int cpuInd;
    ///index for the acceptance rate PropertyCache(*,acceptInd) 
    int acceptInd;
    ///hdf5 file handler
    hid_t h_file;
    ///observables handler
    hid_t h_obs;
    ///communicator to handle communication
    Communicate* myComm;
    //Cummulative energy and weight
    TinyVector<RealType,2>  EPSum;
    ///ostream for the output
    QMCHamiltonian& H;
    ///pointer to the primary ScalarEstimatorBase
    ScalarEstimatorBase* MainEstimator;
    ///save the weights
    Vector<RealType> TotalWeight;
    ///save block averages (scalar data) to Cache 
    Matrix<RealType> AverageCache;
    ///save property data to Cache 
    Matrix<RealType> PropertyCache;
    ///manager of scalar data
    RecordNamedProperty<RealType> BlockAverages;
    ///manager of property data
    RecordNamedProperty<RealType> BlockProperties;
    ///block averages: name to value
    RecordNamedProperty<RealType> TotalAverages;
    ///data accumulated over the blocks
    Matrix<RealType> TotalAveragesData;
    ///index mapping between BlockAverages and TotalAverages
    vector<int> Block2Total;
    ///column map
    std::map<string,int> EstimatorMap;
    ///estimators
    vector<EstimatorType*> Estimators;
    ///Timer
    Timer MyTimer;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
