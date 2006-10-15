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
#include "Estimators/ScalarEstimatorBase.h"
#include "Message/Communicate.h"
#include "Utilities/Timer.h"

namespace qmcplusplus {

  class MCWalkerConifugration;
  class QMCHamiltonian;
  class SimpleFixedNodeBranch;

  /**Class to manage a set of ScalarEstimators */
  class ScalarEstimatorManager: public QMCTraits {

  public:

    typedef ScalarEstimatorBase<RealType> EstimatorType;
    typedef EstimatorType::BufferType     BufferType;
    enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

    ScalarEstimatorManager(QMCHamiltonian& h);
    virtual ~ScalarEstimatorManager();

    ///return the number of ScalarEstimators
    inline int size() const { return Estimators.size();}


    /** add a column with the name
     * @param aname name of the column
     * @return the column index
     *
     * Append a named column. 
     */
    inline int addColumn(const char* aname) {
      return BlockAverages.add(aname);
    }

    /** set the value of the i-th column with a value v
     * @param i column index
     * @param v value 
     */
    inline void setColumn(int i, RealType v) {
      BlockAverages[i] = v;
    }

    inline RealType getColumn(int i) const {
      return  BlockAverages[i];
    }

    int addObservable(const char* aname);

    inline RealType getObservable(int i) const {
      return  TotalAverages[i];
    }

    void getData(int i, vector<RealType>& values);

    /*!\return the index of the newestimator
     *\brief add a new estimator with name aname
     */
    int add(EstimatorType* newestimator, const string& aname);

    ///set the stride for all the estimators
    void setPeriod(int p) { Period=p;}

    ///return a pointer to the estimator aname
    EstimatorType* getEstimator(const string& a);
 
    ///return the average for estimator i
    inline RealType average(int i) const { 
      return Estimators[i]->average(); 
    }

    ///returns a variance for estimator i
    inline RealType variance(int i) const { 
      return Estimators[i]->variance();
    }
  
    /** set the total weight
     *@param w the weight
    inline void setWeight(RealType w) {
      MyData[WEIGHT_INDEX] = w;
    }
     */

    /// start the timer
    inline void startBlock() { MyTimer.restart();}
    /** stop the timer and set the block cpu time
     * @param iter current step
     * @param accept acceptance rate of this block
     */
    inline void stopBlock(RealType accept) { 
      MyData[BLOCK_CPU_INDEX]=MyTimer.elapsed();
      MyData[ACCEPT_RATIO_INDEX]=accept;
      flush();
    }

    void setCollectionMode(bool collect);

    void setAccumulateMode (bool setAccum) {AccumulateBlocks = setAccum;};

    ///process xml tag associated with estimators
    virtual bool put(xmlNodePtr cur);
    virtual void accumulate(MCWalkerConfiguration& W);
    virtual void resetReportSettings(const string& aname, bool append);
    virtual void reportHeader(bool append);
    virtual void flushreport(int iter);
    virtual void report(int iter);
    virtual void flush();
    virtual void finalize();
    virtual void reset();
    virtual void finalize(SimpleFixedNodeBranch& branchEngine);

  protected:

    ///if yes, this estimator will write the data to a file
    bool FileManager;
    ///if yes, the averages are collected over mpi nodes
    bool CollectSum;
    ///if yes, the block averages  are stored in TotalAverages
    bool AccumulateBlocks;
    ///boolean to determine whether a header should be printed.
    bool firstReport;
    ///if true, needs to collect data using irecv.
    bool ManagerNode;
    int BinSize;
    int index;

    //this should be encapsulated
    int myNodeID;
    int numNodes;
    int myGroupID;
    int msgBufferSize;
    ///storage for MPI_Request
    vector<MPI_Request> myRequest;
    ///storage for MPI_Status
    vector<MPI_Status> myStatus;
    ///Data for communication
    vector<BufferType*> RemoteData;

    ///Period to write 
    IndexType Period;
    ///the root file name
    string RootName;
    ///Weight for global collection
    RealType NodeWeight;
    ///Index map for the data maintained by this object
    TinyVector<IndexType,TOTAL_INDEX>  MyIndex;
    ///Data maintained by this object
    TinyVector<RealType,TOTAL_INDEX>  MyData;
    ///Accumulated energy and population
    TinyVector<RealType,2> EPSum;
    ///ostream for the output
    ostream* OutStream;
    ///ostream for the output
    QMCHamiltonian& H;
    ///block averages
    RecordNamedProperty<RealType> BlockAverages;

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
