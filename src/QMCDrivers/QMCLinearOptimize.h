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
#include "Numerics/LinearFit.h"


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
    
    inline bool fitMappedStabilizers(vector<std::pair<RealType,RealType> >& mappedStabilizers, RealType& XS, RealType val=0)
    {
      int nms(0);
      for (int i=0; i<mappedStabilizers.size(); i++) if (mappedStabilizers[i].second==mappedStabilizers[i].second) nms++;
      bool SuccessfulFit(false);
      if (nms>=5)
      {//Quartic fit the stabilizers we have tried and try to choose the best we can
        vector<RealType>  Y(nms), Coefs(5);
        Matrix<RealType> X(nms,5);
        for (int i=0; i<nms; i++)
          if(mappedStabilizers[i].second==mappedStabilizers[i].second)
          {
            X(i,0)=1.0;
            X(i,1)=mappedStabilizers[i].first;
            X(i,2)=X(i,1)*X(i,1);
            X(i,3)=X(i,2)*X(i,1);
            X(i,4)=X(i,3)*X(i,1);
            Y[i]=mappedStabilizers[i].second;
          }
        LinearFit(Y,X,Coefs);

        XS = QuarticMinimum(Coefs);
        val=0;
        for (int i=0; i<5; i++) val+=std::pow(XS,i)*Coefs[i];
        
        SuccessfulFit=true;
        for (int i=0; i<nms; i++)
          if(mappedStabilizers[i].second==mappedStabilizers[i].second)
            if (val>mappedStabilizers[i].second) SuccessfulFit=false;
      }
      else if (nms>=3)
      {//Quadratic fit the stabilizers we have tried and try to choose the best we can
        std::sort(mappedStabilizers.begin(),mappedStabilizers.end());
        vector<RealType>  Y(nms), Coefs(3);
        Matrix<RealType> X(nms,3);
        for (int i=0; i<nms; i++)
          if(mappedStabilizers[i].second==mappedStabilizers[i].second)
          {
            X(i,0)=1.0;
            X(i,1)=mappedStabilizers[i].first;
            X(i,2)=X(i,1)*X(i,1);
            Y[i]=mappedStabilizers[i].second;
          }
        LinearFit(Y,X,Coefs);
        
        //extremum really.
        RealType quadraticMinimum(-0.5*Coefs[1]/Coefs[2]);
        val=quadraticMinimum*quadraticMinimum*Coefs[2]+quadraticMinimum*Coefs[1]+Coefs[0];
//                 RealType dltaBest=std::max(stabilityBase, quadraticMinimum);
//               app_log()<<"smallest XS:      "<<X(0,1)<<endl;
//               app_log()<<"quadraticMinimum: "<<quadraticMinimum<<endl;
        SuccessfulFit=true;
        if (val<Y[0])
          XS = quadraticMinimum;
        else
          SuccessfulFit=false;
      }
    }
    
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
    
    void start();
    void finish();
    //asymmetric generalized EV
    RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev);
    //asymmetric EV
    RealType getLowestEigenvector(Matrix<RealType>& A, vector<RealType>& ev);
    RealType getSplitEigenvectors(int first, int last, Matrix<RealType>& FullLeft, Matrix<RealType>& FullRight, vector<RealType>& FullEV, vector<RealType>& LocalEV, string CSF_Option, bool& CSF_scaled);
    void getNonLinearRange(int& first, int& last);
    bool nonLinearRescale( vector<RealType>& dP, Matrix<RealType>& S);
    RealType getNonLinearRescale( vector<RealType>& dP, Matrix<RealType>& S);
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
