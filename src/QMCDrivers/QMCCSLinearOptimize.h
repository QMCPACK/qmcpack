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
/** @file QMCCSLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/VMC/VMCLinearOptOMP.h"
#include "QMCDrivers/QMCCSLinearOptimizeWFmanagerOMP.h"
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

class QMCCSLinearOptimize: public QMCDriver
{
public:

    ///Constructor.
    QMCCSLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                      QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool);

    ///Destructor
    ~QMCCSLinearOptimize();

    ///Run the Optimization algorithm.
    bool run();
    ///process xml node
    bool put(xmlNodePtr cur);
    void resetComponents(xmlNodePtr cur);
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);
    void setWaveFunctionNode(xmlNodePtr cur)
    {
        wfNode=cur;
    }

private:

    int NumOfVMCWalkers;
    ///Number of iterations maximum before generating new configurations.
    int Max_iterations;
    ///need to know HamiltonianPool to use OMP
    HamiltonianPool& hamPool;
    WaveFunctionPool& psipool;
    ///target cost function to optimize
//     QMCCostFunction* optTarget;
    QMCCSLinearOptimizeWFmanagerOMP* optTarget;
    /// switch to control whether NRCOptimization::lineoptimization() is used or somethign else
    string MinMethod, GEVtype;

    vector<RealType> optdir, optparm;
    RealType stabilizerScale, bigChange, exp0, stepsize;
    RealType Lambda;
    int nstabilizers;
    /// percent variance or H2 to mix in
    RealType w_beta;
    ///Dimension of matrix and number of parameters
    int N,numParams;
    ///vmc engine
    VMCLinearOptOMP* vmcEngine;
    ///xml node to be dumped
    xmlNodePtr wfNode;
    ///xml node for optimizer
    xmlNodePtr optNode;
    ///method for optimization, default best
    string optmethod;
    ///list of files storing configurations
    vector<string> ConfigFile;
    ///Copy Constructor (disabled).
    QMCCSLinearOptimize(const QMCCSLinearOptimize& a): QMCDriver(a),hamPool(a.hamPool), psipool(a.psipool) { }
    ///Copy operator (disabled).
    QMCCSLinearOptimize& operator=(const QMCCSLinearOptimize&)
    {
        return *this;
    }
    bool ValidCostFunction(bool valid){return true;};
    
    void start();
    void finish();
    
    inline int CubicFormula (double a, double b, double c, double d,
            double &x1, double &x2, double &x3)
    {
      double A = b/a;
      double B = c/a;
      double C = d/a;
      double Q = (A*A - 3.0*B)/9.0;
      double R = (2.0*A*A*A - 9.0*A*B + 27.0*C)/54.0;
      //cerr << "Q = " << Q << " R = " << R << "\n";
      if ((R*R) < (Q*Q*Q))
        {
          double theta = std::acos(R/std::sqrt(Q*Q*Q));
          double twosqrtQ = 2.0*std::sqrt(Q);
          double third = 1.0/3.0;
          double thirdA = third * A;
          x1 = -twosqrtQ*std::cos(third*theta) - thirdA;
          x2 = -twosqrtQ*std::cos(third*(theta + 2.0*M_PI)) - thirdA;
          x3 = -twosqrtQ*std::cos(third*(theta - 2.0*M_PI)) - thirdA;
          return 3;
        }
      else 
      {
        double D = -Q*Q*Q + R*R;
        double u = cbrt(-R + std::sqrt(D));
        double v = cbrt(-R - std::sqrt(D));
        double y1 = u+v;
        x1 = y1 - A/3.0;
        return 1;
      }
    }
  
    inline RealType QuarticMinimum (vector<RealType> &coefs)
    {
      double a, b, c, d;
      a = 4.0*coefs[4];
      b = 3.0*coefs[3];
      c = 2.0*coefs[2];
      d = coefs[1];
      double x1, x2, x3;
      int numroots = CubicFormula (a, b, c, d, x1, x2, x3);
      if (numroots == 1)
        return x1;
      else {
        double v1 = coefs[0] + coefs[1]*x1 + coefs[2]*x1*x1 + coefs[3]*x1*x1*x1
    + coefs[4]*x1*x1*x1*x1;
        double v2 = coefs[0] + coefs[1]*x2 + coefs[2]*x2*x2 + coefs[3]*x2*x2*x2
    + coefs[4]*x2*x2*x2*x2;
        double v3 = coefs[0] + coefs[1]*x3 + coefs[2]*x3*x3 + coefs[3]*x3*x3*x3
    + coefs[4]*x3*x3*x3*x3;
        if (v1 < v2 && v1 < v3)
    return x1;
        if (v2 < v1 && v2 < v3)
    return x2;
        if (v3 < v1 && v3 < v2)
    return x3;
        return x1;
      }
    }
    
    RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev);
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
 * $Id: QMCCSLinearOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
