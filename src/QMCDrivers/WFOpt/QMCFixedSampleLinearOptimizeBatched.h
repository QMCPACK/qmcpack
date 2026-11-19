//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCFixedSampleLinearOptimizeBatched.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_BATCHED_H
#define QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_BATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "NRCOptimization.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "QMCDrivers/Optimizers/OptimizerTypes.h"
#include "OutputMatrix.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */


///forward declaration of a cost function
class QMCCostFunctionBase;
class VMCBatched;
class GradientTest;
namespace testing
{
class QMCFixedSampleLinearOptimizeInputTest;
}


class QMCFixedSampleLinearOptimizeBatched : public QMCDriverNew
{
public:
  ///Constructor.
  QMCFixedSampleLinearOptimizeBatched(const ProjectData& project_data,
                                      QMCDriverInput&& qmcdriver_input,
                                      VMCDriverInput&& vmcdriver_input,
                                      WalkerConfigurations& wc,
                                      MCPopulation&& population,
                                      const RefVector<RandomBase<FullPrecRealType>>& rng_refs,
                                      SampleStack& samples,
                                      Communicate& comm);

  ///Destructor
  ~QMCFixedSampleLinearOptimizeBatched() override;

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  ///Run the Optimization algorithm.
  void run() override;
  ///preprocess xml node
  void process(xmlNodePtr cur) override;
  ///process xml node value (parameters for both VMC and OPT) for the actual optimization
  bool processOptXML(xmlNodePtr cur, const std::string& vmcMove, bool reportH5);

  ///common operation to start optimization
  void start();

  using ValueType = QMCTraits::ValueType;
  void descent_start();


  ///common operation to finish optimization, used by the derived classes
  void finish();

  void generateSamples();


private:
  friend class testing::QMCFixedSampleLinearOptimizeInputTest;

  NRCOptimization<RealType> nrc_opt_;

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

  // simple stochastic reconfig
  void stochastic_reconfiguration_conjugate_gradient();

  // perform optimization using a gradient descent algorithm
  void descent_run();

  // Previous linear optimizers ("quartic" and "rescale")
  void previous_linear_methods_run();


  // Perform test of gradients
  void test_run();

  std::unique_ptr<GradientTest> testEngineObj;


  //engine for running various gradient descent based algorithms for optimization
  std::unique_ptr<DescentEngine> descentEngineObj;

  // ------------------------------------
  // Used by legacy linear method algos

  std::vector<RealType> optdir, optparam;

  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;

  RealType param_tol;
  //-------------------------------------

  /// Choice of eigenvalue solver
  std::string eigensolver_;

  ///Number of iterations maximum before generating new configurations.
  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0;
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
  /// SR projection timestep
  RealType sr_tau;
  /// SR regularization parameter
  RealType sr_regularization;
  /// tolerance for CG solution in SR
  RealType sr_tolerance;

  ///name of the current optimization method, updated by processOptXML before run
  std::string MinMethod;
  OptimizerType current_optimizer_type_ = OptimizerType::NONE;

  // Test parameter gradients
  bool doGradientTest;

  // Output Hamiltonian and overlap matrices
  bool do_output_matrices_csv_;

  // Output Hamiltonian and overlap matrices in HDF format
  bool do_output_matrices_hdf_;

  // Flag to open the files on first pass and print header line
  bool output_matrices_initialized_;

  OutputMatrix output_hamiltonian_;
  OutputMatrix output_overlap_;

  // Freeze variational parameters.  Do not update them during each step.
  bool freeze_parameters_;

  bool use_line_search_;

  NewTimer& initialize_timer_;
  NewTimer& generate_samples_timer_;
  NewTimer& build_olv_ham_timer_;
  NewTimer& invert_olvmat_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;
  NewTimer& sr_solver_timer_;

  ///xml node to be dumped
  xmlNodePtr wfNode;

  ParameterSet m_param;

  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;

  ///vmc engine
  std::unique_ptr<VMCBatched> vmcEngine;

  VMCDriverInput vmcdriver_input_;
  SampleStack& samples_;

  /// This is retained in order to construct and reconstruct the vmcEngine.
  const std::optional<EstimatorManagerInput> global_emi_;
};
} // namespace qmcplusplus
#endif
