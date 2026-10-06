//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCFixedSampleLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H

#include "NRCOptimization.h"
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/Optimizers/OptimizerTypes.h"
#include "OutputMatrix.h"

namespace qmcplusplus
{

///forward declaration of a cost function
class QMCCostFunctionBase;
class GradientTest;
class VMC;
namespace testing
{
class QMCFixedSampleLinearOptimizeInputTest;
}

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCFixedSampleLinearOptimize : public QMCDriver, private NRCOptimization<QMCTraits::RealType>
{
public:
  ///Constructor.
  QMCFixedSampleLinearOptimize(const ProjectData& project_data,
                               MCWalkerConfiguration& w,
                               TrialWaveFunction& psi,
                               QMCHamiltonian& h,
                               Communicate*);

  ///Destructor
  ~QMCFixedSampleLinearOptimize() override;

  ///Run the Optimization algorithm.
  void run() override;
  ///preprocess xml node
  bool put(xmlNodePtr cur) override;
  ///process xml node value (parameters for both VMC and OPT) for the actual optimization
  bool processOptXML(xmlNodePtr cur, const std::string& vmcMove, bool reportH5);

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }

private:
  friend class testing::QMCFixedSampleLinearOptimizeInputTest;

  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log() << " Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization "
                   "or reduce the number of cycles. "
                << std::endl;
    return valid;
  }

  // perform the single-shift update, no sample regeneration
  void one_shift_run();


  // Perform test of parameter gradients
  void test_run();

  std::unique_ptr<GradientTest> testEngineObj;


  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0, exp1, stepsize, savedQuadstep;
  std::string StabilizerMethod;
  RealType w_beta;
  /// number of previous steps to orthogonalize to.
  int eigCG;
  /// total number of cg steps per iterations
  int TotalCGSteps;
  /// the previous best identity shift
  RealType bestShift_i;
  /// the previous best overlap shift
  RealType bestShift_s;
  /// current shift_i, shift_s input values
  RealType shift_i_input, shift_s_input;
  /// accept history, remember the last 2 iterations, value 00, 01, 10, 11
  std::bitset<2> accept_history;
  /// Shift_s adjustment base
  RealType shift_s_base;
  //Variables for alternatives to linear method

  //name of the current optimization method, updated by processOptXML before run
  std::string MinMethod;

  //type of the current optimization method, updated by processOptXML before run
  OptimizerType current_optimizer_type_;

  bool doGradientTest;

  // Output Hamiltonian and overlap matrices
  bool do_output_matrices_;

  // Flag to open the files on first pass and print header line
  bool output_matrices_initialized_;

  OutputMatrix output_hamiltonian_;
  OutputMatrix output_overlap_;

  // Freeze variational parameters.  Do not update them during each step.
  bool freeze_parameters_;

  std::vector<RealType> optdir, optparam;
  ///total number of VMC walkers
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///vmc engine
  std::unique_ptr<VMC> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;

  RealType param_tol;

  ///common operation to start optimization, used by the derived classes
  void start();
  ///common operation to finish optimization, used by the derived classes
  void finish();
  void generateSamples();

  NewTimer& generate_samples_timer_;
  NewTimer& initialize_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& involvmat_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;
  Timer t1;
};
} // namespace qmcplusplus
#endif
