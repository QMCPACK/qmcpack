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
#include "Optimize/NRCOptimization.h"

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

class QMCLinearOptimize: public QMCDriver, private NRCOptimization<QMCTraits::RealType>
{
public:

    ///Constructor.
    QMCLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                      QMCHamiltonian& h, HamiltonianPool& hpool);

    ///Destructor
    ~QMCLinearOptimize();

    ///Run the Optimization algorithm.
    bool run();
    ///process xml node
    bool put(xmlNodePtr cur);
    void resetComponents(xmlNodePtr cur);
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);
    RealType Func(Return_t dl);
    void setWaveFunctionNode(xmlNodePtr cur)
    {
        wfNode=cur;
    }

private:

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
    ///stop optimizing if cost function decrease is less than this
    RealType costgradtol;
    ///yes/no applicable only first time
    string SkipSampleGeneration;
    ///need to know HamiltonianPool to use OMP
    HamiltonianPool& hamPool;
    ///target cost function to optimize
    //QMCCostFunction* optTarget;
    QMCCostFunctionBase* optTarget;
    /// switch to control whether NRCOptimization::lineoptimization() is used or somethign else
    string MinMethod, GEVtype, StabilizerMethod, GEVSplit;

    vector<RealType> optdir, optparm;
    RealType allowedCostDifference, stabilizerScale, bigChange, exp0, exp1, savedQuadstep;
    int nstabilizers;
    /// number of previous steps to orthogonalize to.
    int eigCG;
    /// total number of cg steps per iterations
    int  TotalCGSteps;
    /// percent variance or H2 to mix in
    RealType w_beta;
    ///Dimension of matrix and number of parameters
    int N,numParams;
    ///vmc engine
    QMCDriver* vmcEngine;
    ///xml node to be dumped
    xmlNodePtr wfNode;
    ///xml node for optimizer
    xmlNodePtr optNode;
    ///method for optimization, default conjugate gradient
    string optmethod;
    ///list of files storing configurations
    vector<string> ConfigFile;
    ///Copy Constructor (disabled).
    QMCLinearOptimize(const QMCLinearOptimize& a): QMCDriver(a),hamPool(a.hamPool) { }
    ///Copy operator (disabled).
    QMCLinearOptimize& operator=(const QMCLinearOptimize&)
    {
        return *this;
    }
    bool ValidCostFunction(bool valid);
    
    inline bool tooLow(RealType safeValue, RealType CurrentValue)
    {
      RealType lowestCostAllowed=std::min(safeValue-0.03*std::abs(safeValue),safeValue-3.0);
      if (CurrentValue<lowestCostAllowed) return true;
      else return false;
    }
    
    void start();
    void finish();
    
    RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev);
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
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCLinearOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
