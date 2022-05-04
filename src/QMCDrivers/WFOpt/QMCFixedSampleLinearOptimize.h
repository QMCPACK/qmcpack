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

#include "Optimize/NRCOptimization.h"
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#endif
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "QMCDrivers/Optimizers/HybridEngine.h"
#include "OutputMatrix.h"
#include "LinearMethod.h"

namespace qmcplusplus
{

class GradientTest;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCFixedSampleLinearOptimize : public QMCDriver, public LinearMethod, private NRCOptimization<QMCTraits::RealType>
{
public:
  ///Constructor.
  QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm);

  ///Destructor
  ~QMCFixedSampleLinearOptimize() override;

  ///Run the Optimization algorithm.
  bool run() override;
  ///preprocess xml node
  bool put(xmlNodePtr cur) override;
  ///process xml node value (parameters for both VMC and OPT) for the actual optimization
  bool processOptXML(xmlNodePtr cur, const std::string& vmcMove, bool reportH5, bool useGPU);

  RealType Func(RealType dl) override;

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }

private:
  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log() << " Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization "
                   "or reduce the number of cycles. "
                << std::endl;
    return valid;
  }

  // check if the proposed new cost function value is the best available
  bool is_best_cost(const int ii,
                    const std::vector<RealType>& cv,
                    const std::vector<double>& sh,
                    const RealType ic) const;

  // perform the adaptive three-shift update
  bool adaptive_three_shift_run();

  // perform the single-shift update, no sample regeneration
  bool one_shift_run();

  // perform optimization using a gradient descent algorithm
  bool descent_run();

#ifdef HAVE_LMY_ENGINE
  // use hybrid approach of descent and blocked linear method for optimization
  bool hybrid_run();
#endif

  // Perform test of parameter gradients
  bool test_run();

  std::unique_ptr<GradientTest> testEngineObj;


  void solveShiftsWithoutLMYEngine(const std::vector<double>& shifts_i,
                                   const std::vector<double>& shiffts_s,
                                   std::vector<std::vector<RealType>>& parameterDirections);

#ifdef HAVE_LMY_ENGINE
  formic::VarDeps vdeps;
  cqmc::engine::LMYEngine<ValueType>* EngineObj;
#endif

  //engine for running various gradient descent based algorithms for optimization
  std::unique_ptr<DescentEngine> descentEngineObj;

  //engine for controlling a optimization using a hybrid combination of linear method and descent
  std::unique_ptr<HybridEngine> hybridEngineObj;

  // prepare a vector of shifts to try
  std::vector<double> prepare_shifts(const double central_shift) const;

  // previous update
#ifdef HAVE_LMY_ENGINE
  std::vector<formic::ColVec<double>> previous_update;
#endif

  void print_cost_summary_header();
  void print_cost_summary(const double si,
                          const double ss,
                          const RealType mc,
                          const RealType cv,
                          const int ind,
                          const int bi,
                          const bool gu);

  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0, exp1, stepsize, savedQuadstep;
  std::string GEVtype, StabilizerMethod, GEVSplit;
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
  /// number of shifts we will try
  int num_shifts;
  /// the maximum relative change in the cost function for the adaptive three-shift scheme
  RealType max_relative_cost_change;
  ///max amount a parameter may change relative to current wave function weight
  RealType max_param_change;
  /// the tolerance to cost function increases when choosing the best shift in the adaptive shift method
  RealType cost_increase_tol;
  /// the shift_i value that the adaptive shift method should aim for
  RealType target_shift_i;
  ///whether we are targeting an excited state
  std::string targetExcitedStr;
  ///whether we are targeting an excited state
  bool targetExcited;
  ///whether we are doing block algorithm
  std::string block_lmStr;
  ///whether we are doing block algorithm
  bool block_lm;
  ///number of blocks used in block algorithm
  int nblocks;
  ///number of old updates kept
  int nolds;
  ///number of directions kept
  int nkept;
  ///number of samples to do in correlated sampling part
  int nsamp_comp;
  ///the shift to use when targeting an excited state
  RealType omega_shift;
  ///whether to do the first part of block lm
  bool block_first;
  ///whether to do the second part of block lm
  bool block_second;
  ///whether to do the third part of block lm
  bool block_third;

  //Variables for alternatives to linear method

  //name of the current optimization method, updated by processOptXML before run
  std::string MinMethod;

  //type of the previous optimization method, updated by processOptXML before run
  OptimizerType previous_optimizer_type_;

  //type of the current optimization method, updated by processOptXML before run
  OptimizerType current_optimizer_type_;

  //whether to use hybrid method
  bool doHybrid;

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
  std::unique_ptr<QMCDriver> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;

  RealType param_tol;

  ///common operation to start optimization, used by the derived classes
  void start();
#ifdef HAVE_LMY_ENGINE
  void engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                    DescentEngine& descentEngineObj,
                    std::string MinMethod);
#endif
  ///common operation to finish optimization, used by the derived classes
  void finish();
  void generateSamples();

  NewTimer& generate_samples_timer_;
  NewTimer& initialize_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;
  Timer t1;
};
} // namespace qmcplusplus
#endif
