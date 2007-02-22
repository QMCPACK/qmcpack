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

namespace qmcplusplus {

  class MCWalkerConifugration;
  class QMCHamiltonian;

  /**Class to manage a set of ScalarEstimators */
  class ScalarEstimatorManager: public QMCTraits
  {

  public:

    typedef ScalarEstimatorBase           EstimatorType;
    typedef EstimatorType::BufferType     BufferType;
    enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

    ///name of the primary estimator name
    std::string MainEstimatorName;

    ScalarEstimatorManager(QMCHamiltonian& h, Communicate* c=0);
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
    bool put(xmlNodePtr cur);
    void accumulate(MCWalkerConfiguration& W);
    void accumulate(ParticleSet& P, MCWalkerConfiguration::Walker_t& awalker);
    void resetReportSettings(const string& aname, bool append);
    void resetReportSettings(bool append);
    void reportHeader(bool append);
    void flushreport(int iter);
    void report(int iter);
    void flush();
    void reset();
    void getEnergyAndWeight(RealType& e, RealType& w);
    void finalize();

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
    ///if true, needs to cancel it
    bool HasIrecvIssued;
    int BinSize;
    int index;

    /** communicator to handle communication
     */
    Communicate* myComm;
    //this should be encapsulated
    int myNodeID;
    int numNodes;
    int msgBufferSize;
    ///storage for MPI_Request
    vector<Communicate::mpi_request_type> myRequest;
    ///storage for MPI_Status
    vector<Communicate::mpi_status_type> myStatus;
    ///Data for communication
    vector<BufferType*> RemoteData;

    ///the root file name
    string RootName;
    ///Weight for global collection
    RealType NodeWeight;
    ///Index map for the data maintained by this object
    TinyVector<IndexType,TOTAL_INDEX>  MyIndex;
    ///Data maintained by this object
    TinyVector<RealType,TOTAL_INDEX>  MyData;
    ///Cummulative energy and weight
    TinyVector<RealType,2>  EPSum;
    ///ostream for the output
    ostream* OutStream;
    ///ostream for the output
    QMCHamiltonian& H;
    ///pointer to the primary ScalarEstimatorBase
    ScalarEstimatorBase* MainEstimator;
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
