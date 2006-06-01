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
#include "Utilities/Timer.h"

namespace qmcplusplus {

  class MCWalkerConifugration;
  class QMCHamiltonian;

  /**Class to manage a set of ScalarEstimators */
  class ScalarEstimatorManager: public QMCTraits {

  public:

    typedef ScalarEstimatorBase<RealType> EstimatorType;
    enum { WEIGHT_INDEX=0, BLOCK_CPU_INDEX, ACCEPT_RATIO_INDEX, TOTAL_INDEX};

    ScalarEstimatorManager(QMCHamiltonian& h);
    ~ScalarEstimatorManager();

    ///return the number of ScalarEstimators
    inline int size() const { return Estimators.size();}

    ///process xml tag associated with estimators
    bool put(xmlNodePtr cur);

    void accumulate(MCWalkerConfiguration& W);
    void resetReportSettings(const string& aname, bool append);
    void reportHeader(bool append);
    void flushreport(int iter);
    void report(int iter);
    void flush();
    void finalize();
    void reset();
  
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


  private:

    ///the root file name
    string RootName;
    ///if yes, this estimator will write the data to a file
    bool FileManager;
    ///if yes, the averages are collected over mpi nodes
    bool CollectSum;
    ///Period to write 
    IndexType Period;
    ///Weight for global collection
    RealType NodeWeight;
    ///Index map for the data maintained by this object
    TinyVector<IndexType,TOTAL_INDEX>  MyIndex;
    ///Data maintained by this object
    TinyVector<RealType,TOTAL_INDEX>  MyData;
    ///ostream for the output
    ostream* OutStream;
    ///ostream for the output
    QMCHamiltonian& H;
    ///block averages
    RecordNamedProperty<RealType> BlockAverages;
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
