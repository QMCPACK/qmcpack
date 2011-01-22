//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file QMCLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"
#include "QMCApp/WaveFunctionPool.h"


namespace qmcplusplus
{

///forward declaration of a cost function
class QMCCostFunctionBase;
class HamiltonianPool;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCLinearOptimize: public QMCDriver
{
public:

    ///Constructor.
    QMCLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                      QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);
                      
    ///Destructor
    virtual ~QMCLinearOptimize();

    ///Run the Optimization algorithm.
    virtual bool run()=0;
    ///process xml node
    virtual bool put(xmlNodePtr cur);
    void resetComponents(xmlNodePtr cur);
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);
    void setWaveFunctionNode(xmlNodePtr cur)
    {
        wfNode=cur;
    }
    
    vector<RealType> optdir, optparm;
    ///index to denote the partition id
    int PartID;
    ///total number of partitions that will share a set of configuratons
    int NumParts;
    ///total number of Warmup Blocks
    int WarmupBlocks;
    ///total number of Warmup Blocks
    int NumOfVMCWalkers;
    ///Number of iterations maximum before generating new configurations.
    int Max_iterations;
    ///need to know HamiltonianPool to use OMP
    HamiltonianPool& hamPool;
    ///target cost function to optimize
    QMCCostFunctionBase* optTarget;
    ///Dimension of matrix and number of parameters
    int N,numParams;
    ///vmc engine
    QMCDriver* vmcEngine;
    ///xml node to be dumped
    xmlNodePtr wfNode;
    ///xml node for optimizer
    xmlNodePtr optNode;
    ///list of files storing configurations
    vector<string> ConfigFile;
    
    RealType param_tol;
    
    inline bool tooLow(RealType safeValue, RealType CurrentValue)
    {
      RealType lowestCostAllowed=std::min(safeValue-0.03*std::abs(safeValue),safeValue-3.0);
      if (CurrentValue<lowestCostAllowed) return true;
      else return false;
    }
    
    void start();
    void finish();
    //asymmetric generalized EV
    RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev);
    //asymmetric EV
    RealType getLowestEigenvector(Matrix<RealType>& A, vector<RealType>& ev);
    RealType getSplitEigenvectors(int first, int last, Matrix<RealType>& FullLeft, Matrix<RealType>& FullRight, vector<RealType>& FullEV, vector<RealType>& LocalEV, string CSF_Option, bool& CSF_scaled);
    void getNonLinearRange(int& first, int& last);
    bool nonLinearRescale( vector<RealType>& dP, Matrix<RealType> S);
    RealType getNonLinearRescale( vector<RealType>& dP, Matrix<RealType> S);
    void generateSamples();
    void add_timers(vector<NewTimer*>& timers);
    vector<NewTimer*> myTimers;
    Timer t1;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
