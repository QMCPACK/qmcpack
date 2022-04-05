//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H

#include <memory>
#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Numerics/LinearFit.h"
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#endif


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

class QMCLinearOptimize : public QMCDriver
{
public:
  ///Constructor.
  QMCLinearOptimize(MCWalkerConfiguration& w,
                    TrialWaveFunction& psi,
                    QMCHamiltonian& h,
                    Communicate* comm,
                    const std::string& QMC_driver_type = "QMCLinearOptimize");

  ///Destructor
  ~QMCLinearOptimize() override = default;

  ///Run the Optimization algorithm.
  bool run() override = 0;
  ///process xml node
  bool put(xmlNodePtr cur) override;
  ///add a configuration file to the list of files
  void addConfiguration(const std::string& a);
  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  std::vector<RealType> optdir, optparam;
  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  ///total number of VMC walkers
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///Dimension of matrix and number of parameters
  int N, numParams;
  ///vmc engine
  std::unique_ptr<QMCDriver> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;
  ///list of files storing configurations
  std::vector<std::string> ConfigFile;

  RealType param_tol;

  inline bool tooLow(RealType safeValue, RealType CurrentValue)
  {
    RealType lowestCostAllowed = std::min(safeValue - 0.03 * std::abs(safeValue), safeValue - 3.0);
    if (CurrentValue < lowestCostAllowed)
      return true;
    else
      return false;
  }

  bool fitMappedStabilizers(std::vector<std::pair<RealType, RealType>>& mappedStabilizers,
                            RealType& XS,
                            RealType& val,
                            RealType tooBig);

  inline int CubicFormula(double a, double b, double c, double d, double& x1, double& x2, double& x3)
  {
    double A = b / a;
    double B = c / a;
    double C = d / a;
    double Q = (A * A - 3.0 * B) / 9.0;
    double R = (2.0 * A * A * A - 9.0 * A * B + 27.0 * C) / 54.0;
    //cerr << "Q = " << Q << " R = " << R << "\n";
    if ((R * R) < (Q * Q * Q))
    {
      double theta    = std::acos(R / std::sqrt(Q * Q * Q));
      double twosqrtQ = 2.0 * std::sqrt(Q);
      double third    = 1.0 / 3.0;
      double thirdA   = third * A;
      x1              = -twosqrtQ * std::cos(third * theta) - thirdA;
      x2              = -twosqrtQ * std::cos(third * (theta + 2.0 * M_PI)) - thirdA;
      x3              = -twosqrtQ * std::cos(third * (theta - 2.0 * M_PI)) - thirdA;
      return 3;
    }
    else
    {
      double D  = -Q * Q * Q + R * R;
      double u  = cbrt(-R + std::sqrt(D));
      double v  = cbrt(-R - std::sqrt(D));
      double y1 = u + v;
      x1        = y1 - A / 3.0;
      return 1;
    }
  }

  inline RealType QuarticMinimum(std::vector<RealType>& coefs)
  {
    double a, b, c, d;
    a = 4.0 * coefs[4];
    b = 3.0 * coefs[3];
    c = 2.0 * coefs[2];
    d = coefs[1];
    double x1, x2, x3;
    int numroots = CubicFormula(a, b, c, d, x1, x2, x3);
    if (numroots == 1)
      return x1;
    else
    {
      double v1 =
          coefs[0] + coefs[1] * x1 + coefs[2] * x1 * x1 + coefs[3] * x1 * x1 * x1 + coefs[4] * x1 * x1 * x1 * x1;
      double v2 =
          coefs[0] + coefs[1] * x2 + coefs[2] * x2 * x2 + coefs[3] * x2 * x2 * x2 + coefs[4] * x2 * x2 * x2 * x2;
      double v3 =
          coefs[0] + coefs[1] * x3 + coefs[2] * x3 * x3 + coefs[3] * x3 * x3 * x3 + coefs[4] * x3 * x3 * x3 * x3;
      if (v1 < v2 && v1 < v3)
        return x1;
      if (v2 < v1 && v2 < v3)
        return x2;
      if (v3 < v1 && v3 < v2)
        return x3;
      return x1;
    }
  }

  ///common operation to start optimization, used by the derived classes
  void start();
#ifdef HAVE_LMY_ENGINE
  void engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                    DescentEngine& descentEngineObj,
                    std::string MinMethod);
#endif
  ///common operation to finish optimization, used by the derived classes
  void finish();
  //asymmetric generalized EV
  RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, std::vector<RealType>& ev);
  //asymmetric EV
  RealType getLowestEigenvector(Matrix<RealType>& A, std::vector<RealType>& ev);
  void getNonLinearRange(int& first, int& last);
  void orthoScale(std::vector<RealType>& dP, Matrix<RealType>& S);
  bool nonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S);
  RealType getNonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S);
  void generateSamples();

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }
  NewTimer& generate_samples_timer_;
  NewTimer& initialize_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;
  Timer t1;
};
} // namespace qmcplusplus
#endif
