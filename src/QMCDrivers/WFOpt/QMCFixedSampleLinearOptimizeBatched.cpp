//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCFixedSampleLinearOptimizeBatched.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"
#include "QMCDrivers/WFOpt/GradientTest.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCDrivers/WFOpt/ConjugateGradient.h"
#include "Concurrency/Info.hpp"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "EstimatorInputDelegates.h"
#include "Message/UniformCommunicateError.h"
#include "Numerics/DeterminantOperators.h"
#include "LinearMethod.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <stdexcept>


namespace qmcplusplus
{
using MatrixOperators::product;


QMCFixedSampleLinearOptimizeBatched::QMCFixedSampleLinearOptimizeBatched(
    const ProjectData& project_data,
    QMCDriverInput&& qmcdriver_input,
    VMCDriverInput&& vmcdriver_input,
    WalkerConfigurations& wc,
    MCPopulation&& population,
    const RefVector<RandomBase<FullPrecRealType>>& rng_refs,
    SampleStack& samples,
    Communicate* comm)
    : QMCDriverNew(
          project_data,
          std::move(qmcdriver_input),
          nullptr, // this class is not a real QMCDriverNew as far as I can tell so we don't give it the actual EM
          wc,
          std::move(population),
          rng_refs,
          "QMCLinearOptimizeBatched::",
          comm,
          "QMCLinearOptimizeBatched"),
      Max_iterations(1),
      param_tol(1e-4),
      nstabilizers(3),
      stabilizerScale(2.0),
      bigChange(50),
      exp0(-16),
      bestShift_i(-1.0),
      bestShift_s(-1.0),
      shift_i_input(0.01),
      shift_s_input(1.00),
      accept_history(3),
      shift_s_base(4.0),
      sr_tau(0.01),
      sr_regularization(0.01),
      sr_tolerance(1e-6),
      MinMethod("OneShiftOnly"),
      do_output_matrices_csv_(false),
      do_output_matrices_hdf_(false),
      output_matrices_initialized_(false),
      freeze_parameters_(false),
      use_line_search_(false),
      initialize_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Initialize", timer_level_medium)),
      generate_samples_timer_(createGlobalTimer("QMCLinearOptimizeBatched::generateSamples", timer_level_medium)),
      build_olv_ham_timer_(createGlobalTimer("QMCLinearOptimizedBatched::build_ovl_ham_matrix", timer_level_medium)),
      invert_olvmat_timer_(createGlobalTimer("QMCLinearOptimizedBatched::invert_ovlMat", timer_level_medium)),
      eigenvalue_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Eigenvalue", timer_level_medium)),
      line_min_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Line_Minimization", timer_level_medium)),
      cost_function_timer_(createGlobalTimer("QMCLinearOptimizeBatched::CostFunction", timer_level_medium)),
      sr_solver_timer_(createGlobalTimer("QMCLinearOptimizeBatched::StochasticReconfiguration", timer_level_medium)),
      wfNode(NULL),
      vmcdriver_input_(vmcdriver_input),
      samples_(samples)
{
  //set the optimization flag
  qmc_driver_mode_.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  m_param.add(MinMethod, "MinMethod");
  m_param.add(Max_iterations, "max_its");
  m_param.add(nstabilizers, "nstabilizers");
  m_param.add(stabilizerScale, "stabilizerscale");
  m_param.add(bigChange, "bigchange");
  m_param.add(exp0, "exp0");
  m_param.add(param_tol, "alloweddifference");
  m_param.add(shift_i_input, "shift_i");
  m_param.add(shift_s_input, "shift_s");
  m_param.add(sr_tau, "sr_tau");
  m_param.add(sr_regularization, "sr_regularization");
  m_param.add(sr_tolerance, "sr_tolerance");
}

/** Clean up the vector */
QMCFixedSampleLinearOptimizeBatched::~QMCFixedSampleLinearOptimizeBatched() = default;

void QMCFixedSampleLinearOptimizeBatched::start()
{
  generateSamples();

  optTarget->setWaveFunctionNode(wfNode);

  {
    app_log() << std::endl
              << "*****************************" << std::endl
              << "Compute parameter derivatives" << std::endl
              << "*****************************" << std::endl
              << std::endl;
    ScopedTimer local(initialize_timer_);
    Timer t_deriv;
    optTarget->getConfigurations("");
    optTarget->setRng(rngs_);
    NullEngineHandle handle;
    if (current_optimizer_type_ == OptimizerType::STOCHASTIC_RECONFIGURATION_CG)
      optTarget->checkConfigurationsSR(handle);
    else
      optTarget->checkConfigurations(handle);
    app_log() << "  Execution time (derivatives) = " << std::setprecision(4) << t_deriv.elapsed() << std::endl;
  }
}

void QMCFixedSampleLinearOptimizeBatched::descent_start()
{
  app_log() << "entering descent_start function" << std::endl;
  DescentEngineHandle handle(*descentEngineObj);

  // generate samples
  generate_samples_timer_.start();
  generateSamples();
  generate_samples_timer_.stop();

  // store active number of walkers
  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;

  // reset the root name
  optTarget->setRootName(get_root_name());
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "     Reading configurations from h5FileRoot " << std::endl;

  // get configuration from the previous run
  Timer t1;
  initialize_timer_.start();
  optTarget->getConfigurations("");
  optTarget->setRng(rngs_);
  optTarget->checkConfigurations(handle);

  initialize_timer_.stop();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;
  app_log() << R"(<opt stage="main" walkers=")" << optTarget->getNumSamples() << "\">" << std::endl;
}


void QMCFixedSampleLinearOptimizeBatched::finish()
{
  if (optTarget->reportH5)
    optTarget->reportParametersH5();
  optTarget->reportParameters();
}

void QMCFixedSampleLinearOptimizeBatched::generateSamples()
{
  ScopedTimer local(generate_samples_timer_);
  app_log() << std::endl
            << "******************" << std::endl
            << "Generating samples" << std::endl
            << "******************" << std::endl
            << std::endl;
  samples_.resetSampleCount();

  Timer t_gen;
  vmcEngine->run();
  app_log() << "  Execution time (sampling) = " << std::setprecision(4) << t_gen.elapsed() << std::endl;

  //reset the rootname
  h5_file_root_ = get_root_name();
  optTarget->setRootName(get_root_name());
}

void QMCFixedSampleLinearOptimizeBatched::run()
{
  if (do_output_matrices_csv_ && !output_matrices_initialized_)
  {
    const int numParams = optTarget->getNumParams();
    const int N         = numParams + 1;
    output_overlap_.init_file(get_root_name(), "ovl", N);
    output_hamiltonian_.init_file(get_root_name(), "ham", N);
    output_matrices_initialized_ = true;
  }

  if (doGradientTest)
  {
    app_log() << "Doing gradient test run" << std::endl;
    test_run();
  }
  else if (current_optimizer_type_ == OptimizerType::DESCENT)
    descent_run();
  else if (current_optimizer_type_ == OptimizerType::ONESHIFTONLY)
    one_shift_run();
  else if (current_optimizer_type_ == OptimizerType::STOCHASTIC_RECONFIGURATION_CG)
    stochastic_reconfiguration_conjugate_gradient();
  else
    previous_linear_methods_run();
}

void QMCFixedSampleLinearOptimizeBatched::test_run()
{
  // generate samples and compute weights, local energies, and derivative vectors
  start();

  testEngineObj->run(*optTarget, get_root_name());

  finish();
}

void QMCFixedSampleLinearOptimizeBatched::previous_linear_methods_run()
{
  start();
  bool Valid(true);
  int Total_iterations(0);
  //size of matrix
  const int numParams = optTarget->getNumParams();
  const int N         = numParams + 1;
  //   where we are and where we are pointing
  std::vector<RealType> currentParameterDirections(N, 0);
  std::vector<RealType> currentParameters(numParams, 0);
  std::vector<RealType> bestParameters(numParams, 0);
  for (int i = 0; i < numParams; i++)
    bestParameters[i] = currentParameters[i] = std::real(optTarget->Params(i));
  //   proposed direction and new parameters
  optdir.resize(numParams, 0);
  optparam.resize(numParams, 0);

  auto costfunc_evaluator = [this](RealType dl) {
    for (int i = 0; i < optparam.size(); i++)
      optTarget->Params(i) = optparam[i] + dl * optdir[i];
    auto effective_weight = optTarget->correlatedSampling(false);
    nrc_opt_.validFuncVal = optTarget->isEffectiveWeightValid(effective_weight);
    return optTarget->computedCost();
  };

  while (Total_iterations < Max_iterations)
  {
    Total_iterations += 1;
    app_log() << "Iteration: " << Total_iterations << "/" << Max_iterations << std::endl;
    if (!ValidCostFunction(Valid))
      continue;
    //this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
    //     reset params if necessary
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currentParameters[i];
    cost_function_timer_.start();
    auto effective_weight = optTarget->correlatedSampling(true);
    RealType lastCost(optTarget->computedCost());
    cost_function_timer_.stop();
    //     if cost function is currently invalid continue
    Valid = optTarget->isEffectiveWeightValid(effective_weight);
    if (!ValidCostFunction(Valid))
      continue;
    RealType newCost(lastCost);
    RealType startCost(lastCost);
    Matrix<RealType> Left(N, N);
    Matrix<RealType> Right(N, N);
    Matrix<RealType> S(N, N);
    //     stick in wrong matrix to reduce the number of matrices we need by 1.( Left is actually stored in Right, & vice-versa)
    optTarget->fillOverlapHamiltonianMatrices(Right, Left);
    S.copy(Left);
    bool apply_inverse(true);
    if (apply_inverse)
    {
      Matrix<RealType> RightT(Left);
      invert_matrix(RightT, false);
      Left = 0;
      product(RightT, Right, Left);
      //       Now the left matrix is the Hamiltonian with the inverse of the overlap applied ot it.
    }
    //Find largest off-diagonal element compared to diagonal element.
    //This gives us an idea how well conditioned it is, used to stabilize.
    RealType od_largest(0);
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        od_largest = std::max(std::max(od_largest, std::abs(Left(i, j)) - std::abs(Left(i, i))),
                              std::abs(Left(i, j)) - std::abs(Left(j, j)));
    app_log() << "od_largest " << od_largest << std::endl;
    //if(od_largest>0)
    //  od_largest = std::log(od_largest);
    //else
    //  od_largest = -1e16;
    //if (od_largest<stabilityBase)
    //  stabilityBase=od_largest;
    //else
    //  stabilizerScale = std::max( 0.2*(od_largest-stabilityBase)/nstabilizers, stabilizerScale);
    app_log() << "  stabilityBase " << stabilityBase << std::endl;
    app_log() << "  stabilizerScale " << stabilizerScale << std::endl;
    int failedTries(0);
    bool acceptedOneMove(false);
    for (int stability = 0; stability < nstabilizers; stability++)
    {
      bool goodStep(true);
      //       store the Hamiltonian matrix in Right
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
          Right(i, j) = Left(j, i);
      RealType XS(stabilityBase + stabilizerScale * (failedTries + stability));
      for (int i = 1; i < N; i++)
        Right(i, i) += std::exp(XS);
      app_log() << "  Using XS:" << XS << " " << failedTries << " " << stability << std::endl;
      {
        ScopedTimer local(eigenvalue_timer_);
        LinearMethod::getLowestEigenvector(Right, currentParameterDirections);
        nrc_opt_.Lambda = LinearMethod::getNonLinearRescale(currentParameterDirections, S, *optTarget);
      }
      //       biggest gradient in the parameter direction vector
      RealType bigVec(0);
      for (int i = 0; i < numParams; i++)
        bigVec = std::max(bigVec, std::abs(currentParameterDirections[i + 1]));
      //       this can be overwritten during the line minimization
      RealType evaluated_cost(startCost);
      if (MinMethod == "rescale")
      {
        if (std::abs(nrc_opt_.Lambda * bigVec) > bigChange)
        {
          goodStep = false;
          app_log() << "  Failed Step. Magnitude of largest parameter change: " << std::abs(nrc_opt_.Lambda * bigVec)
                    << std::endl;
          if (stability == 0)
          {
            failedTries++;
            stability--;
          }
          else
            stability = nstabilizers;
        }
        for (int i = 0; i < numParams; i++)
          optTarget->Params(i) = currentParameters[i] + nrc_opt_.Lambda * currentParameterDirections[i + 1];
      }
      else
      {
        for (int i = 0; i < numParams; i++)
          optparam[i] = currentParameters[i];
        for (int i = 0; i < numParams; i++)
          optdir[i] = currentParameterDirections[i + 1];
        nrc_opt_.TOL              = param_tol / bigVec;
        nrc_opt_.AbsFuncTol       = true;
        nrc_opt_.largeQuarticStep = bigChange / bigVec;
        nrc_opt_.LambdaMax        = 0.5 * nrc_opt_.Lambda;
        line_min_timer_.start();
        if (MinMethod == "quartic")
        {
          int npts(7);
          nrc_opt_.quadstep         = nrc_opt_.stepsize * nrc_opt_.Lambda;
          nrc_opt_.largeQuarticStep = bigChange / bigVec;
          Valid                     = nrc_opt_.lineoptimization3(costfunc_evaluator, npts, evaluated_cost);
        }
        else
          Valid = nrc_opt_.lineoptimization2(costfunc_evaluator);
        line_min_timer_.stop();
        RealType biggestParameterChange = bigVec * std::abs(nrc_opt_.Lambda);
        if (biggestParameterChange > bigChange)
        {
          goodStep = false;
          failedTries++;
          app_log() << "  Failed Step. Largest LM parameter change:" << biggestParameterChange << std::endl;
          if (stability == 0)
            stability--;
          else
            stability = nstabilizers;
        }
        else
        {
          for (int i = 0; i < numParams; i++)
            optTarget->Params(i) = optparam[i] + nrc_opt_.Lambda * optdir[i];
          app_log() << "  Good Step. Largest LM parameter change:" << biggestParameterChange << std::endl;
        }
      }

      if (goodStep)
      {
        // 	this may have been evaluated already
        // 	newCost=evaluated_cost;
        //get cost at new minimum
        auto effective_weight = optTarget->correlatedSampling(false);
        newCost               = optTarget->computedCost();
        app_log() << " OldCost: " << lastCost << " NewCost: " << newCost << " Delta Cost:" << (newCost - lastCost)
                  << std::endl;
        optTarget->printEstimates();
        //                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
        // mmorales
        Valid = optTarget->isEffectiveWeightValid(effective_weight);
        //if (MinMethod!="rescale" && !ValidCostFunction(Valid))
        if (!ValidCostFunction(Valid))
        {
          goodStep = false;
          app_log() << "  Good Step, but cost function invalid" << std::endl;
          failedTries++;
          if (stability > 0)
            stability = nstabilizers;
          else
            stability--;
        }
        if (newCost < lastCost && goodStep)
        {
          //Move was acceptable
          for (int i = 0; i < numParams; i++)
            bestParameters[i] = std::real(optTarget->Params(i));
          lastCost        = newCost;
          acceptedOneMove = true;
          if (std::abs(newCost - lastCost) < 1e-4)
          {
            failedTries++;
            stability = nstabilizers;
            continue;
          }
        }
        else if (stability > 0)
        {
          failedTries++;
          stability = nstabilizers;
          continue;
        }
      }
      app_log().flush();
      if (failedTries > 20)
        break;
      //APP_ABORT("QMCFixedSampleLinearOptimizeBatched::run TOO MANY FAILURES");
    }

    if (acceptedOneMove)
    {
      app_log() << "Setting new Parameters" << std::endl;
      for (int i = 0; i < numParams; i++)
        optTarget->Params(i) = bestParameters[i];
    }
    else
    {
      app_log() << "Reverting to old Parameters" << std::endl;
      for (int i = 0; i < numParams; i++)
        optTarget->Params(i) = currentParameters[i];
    }
    app_log().flush();
  }

  finish();

}

/** Parses the xml input file for parameter definitions for the wavefunction
* optimization.
* @param q current xmlNode
* @return true if successful
*/
void QMCFixedSampleLinearOptimizeBatched::process(xmlNodePtr q)
{
  std::string vmcMove("pbyp");
  std::string ReportToH5("no");
  std::string OutputMatrices("no");
  std::string OutputMatricesHDF("no");
  std::string FreezeParameters("no");
  std::string UseLineSearch("no");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(vmcMove, "move");
  oAttrib.add(ReportToH5, "hdf5");

  m_param.add(OutputMatrices, "output_matrices_csv", {"no", "yes"});
  m_param.add(OutputMatricesHDF, "output_matrices_hdf", {"no", "yes"});
  m_param.add(FreezeParameters, "freeze_parameters", {"no", "yes"});
  m_param.add(UseLineSearch, "line_search", {"no", "yes"});

  m_param.add(eigensolver_, "eigensolver",
              {
                  "inverse", // Inverse + nonsymmetric eigenvalue solver
                  "general"  // General eigenvalue problem solver
              });

  oAttrib.put(q);
  m_param.put(q);

  do_output_matrices_csv_ = (OutputMatrices == "yes");
  do_output_matrices_hdf_ = (OutputMatricesHDF == "yes");
  freeze_parameters_      = (FreezeParameters == "yes");
  use_line_search_        = (UseLineSearch == "yes");

  // Use freeze_parameters with output_matrices to generate multiple lines in the output with
  // the same parameters so statistics can be computed in post-processing.

  if (freeze_parameters_)
  {
    app_log() << std::endl;
    app_warning() << "  The option 'freeze_parameters' is enabled.  Variational parameters will not be updated.  This "
                     "run will not perform variational parameter optimization!"
                  << std::endl;
    app_log() << std::endl;
  }


  doGradientTest = false;
  processChildren(q, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "optimize")
    {
      const std::string att(getXMLAttributeValue(element, "method"));
      if (!att.empty() && att == "gradient_test")
      {
        GradientTestInput test_grad_input;
        test_grad_input.readXML(element);
        if (!testEngineObj)
          testEngineObj = std::make_unique<GradientTest>(std::move(test_grad_input));
        doGradientTest = true;
        MinMethod      = "gradient_test";
      }
      else
      {
        std::stringstream error_msg;
        app_log() << "Unknown or missing 'method' attribute in optimize tag: " << att << "\n";
        throw UniformCommunicateError(error_msg.str());
      }
    }
  });


  processOptXML(q, vmcMove, ReportToH5 == "yes");
}

bool QMCFixedSampleLinearOptimizeBatched::processOptXML(xmlNodePtr opt_xml,
                                                        const std::string& vmcMove,
                                                        bool reportH5)
{
  m_param.put(opt_xml);

  auto iter = OptimizerNames.find(MinMethod);
  if (iter == OptimizerNames.end())
    throw std::runtime_error("Unknown MinMethod!\n");
  current_optimizer_type_ = OptimizerNames.at(MinMethod);

  if (current_optimizer_type_ == OptimizerType::DESCENT && !descentEngineObj)
    descentEngineObj = std::make_unique<DescentEngine>(myComm, opt_xml);

  // check shift sanity
  if (shift_i_input <= 0.0)
    throw std::runtime_error("shift_i must be positive in QMCFixedSampleLinearOptimizeBatched::put");
  if (shift_s_input <= 0.0)
    throw std::runtime_error("shift_s must be positive in QMCFixedSampleLinearOptimizeBatched::put");

  // if this is the first time this function has been called, set the initial shifts
  if (current_optimizer_type_ == OptimizerType::ONESHIFTONLY ||
      current_optimizer_type_ == OptimizerType::STOCHASTIC_RECONFIGURATION_CG)
    bestShift_i = shift_i_input;
  if (bestShift_s < 0.0)
    bestShift_s = shift_s_input;

  xmlNodePtr qsave = opt_xml;
  xmlNodePtr cur   = qsave->children;
  int pid          = OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    cur = cur->next;
  }

  // Destroy old object to stop timer to correctly order timer with object lifetime scope
  vmcEngine.reset(nullptr);

  // Explicitly copy the driver input objects since they will be used to instantiate the VMCEngine repeatedly.
  QMCDriverInput qmcdriver_input_copy = qmcdriver_input_;
  VMCDriverInput vmcdriver_input_copy = vmcdriver_input_;
  qmcdriver_input_copy.readXML(opt_xml);
  vmcdriver_input_copy.readXML(opt_xml);


  // create VMC engine
  vmcEngine =
      std::make_unique<VMCBatched>(project_data_, std::move(qmcdriver_input_copy), nullptr,
                                   std::move(vmcdriver_input_copy), walker_configs_ref_,
                                   MCPopulation(myComm->size(), myComm->rank(), population_.get_golden_electrons(),
                                                population_.get_golden_twf(), population_.get_golden_hamiltonian()),
                                   rngs_, samples_, myComm);

  vmcEngine->setUpdateMode(vmcMove[0] == 'p');


  bool AppendRun = false;
  vmcEngine->setStatus(get_root_name(), h5_file_root_, AppendRun);
  vmcEngine->process(qsave);

  vmcEngine->enable_sample_collection();

  auto& qmcdriver_input = vmcEngine->getQMCDriverInput();
  QMCDriverNew::AdjustedWalkerCounts awc =
      adjustGlobalWalkerCount(*myComm, walker_configs_ref_.getActiveWalkers(), qmcdriver_input_.get_total_walkers(),
                              qmcdriver_input_.get_walkers_per_rank(), 1.0,
                              determineNumCrowds(qmcdriver_input_.get_num_crowds(), rngs_.size()));


  bool success = true;
  //allways reset optTarget
  optTarget = std::make_unique<QMCCostFunctionBatched>(population_.get_golden_electrons(), population_.get_golden_twf(),
                                                       population_.get_golden_hamiltonian(), samples_,
                                                       awc.walkers_per_crowd, myComm);
  optTarget->setStream(&app_log());
  if (reportH5)
    optTarget->reportH5 = true;
  success = optTarget->put(qsave);

  return success;
}
void QMCFixedSampleLinearOptimizeBatched::one_shift_run()
{
  // ensure the cost function is set to compute derivative vectors
  optTarget->setneedGrads(true);

  // generate samples and compute weights, local energies, and derivative vectors
  start();

  // get number of optimizable parameters
  const int numParams = optTarget->getNumParams();

  // get dimension of the linear method matrices
  const int N = numParams + 1;

  // prepare vectors to hold the initial and current parameters
  std::vector<RealType> currentParameters(numParams, 0.0);

  // initialize the initial and current parameter vectors
  for (int i = 0; i < numParams; i++)
    currentParameters.at(i) = std::real(optTarget->Params(i));

  // prepare vectors to hold the parameter update directions for each shift
  std::vector<RealType> parameterDirections;
  parameterDirections.assign(N, 0.0);

  // compute the initial cost
  const RealType initCost = optTarget->computedCost();

  // allocate the matrices we will need
  Matrix<RealType> ovlMat(N, N);
  ovlMat = 0.0;
  Matrix<RealType> hamMat(N, N);
  hamMat = 0.0;
  Matrix<RealType> invMat(N, N);
  invMat = 0.0;
  Matrix<RealType> prdMat(N, N);
  prdMat = 0.0;

  // for outputing matrices and eigenvalue/vectors to disk
  hdf_archive hout;

  {
    ScopedTimer local(build_olv_ham_timer_);
    Timer t_build_matrices;
    // say what we are doing
    app_log() << std::endl
              << "*****************************************" << std::endl
              << "Building overlap and Hamiltonian matrices" << std::endl
              << "*****************************************" << std::endl;

    // build the overlap and hamiltonian matrices
    optTarget->fillOverlapHamiltonianMatrices(hamMat, ovlMat);

    if (do_output_matrices_csv_)
    {
      output_overlap_.output(ovlMat);
      output_hamiltonian_.output(hamMat);
    }

    if (do_output_matrices_hdf_ && is_manager())
    {
      std::string newh5 = get_root_name() + ".linear_matrices.h5";
      hout.create(newh5, H5F_ACC_TRUNC);
      hout.write(ovlMat, "overlap");
      hout.write(hamMat, "Hamiltonian");
      hout.write(bestShift_i, "bestShift_i");
      hout.write(bestShift_s, "bestShift_s");
    }
    app_log() << "  Execution time (building matrices) = " << std::setprecision(4) << t_build_matrices.elapsed()
              << std::endl;
  }

  // Solve eigenvalue problem on one rank.
  if (is_manager())
  {
    ScopedTimer local(eigenvalue_timer_);
    Timer t_eigen;
    app_log() << std::endl
              << "**************************" << std::endl
              << "Solving eigenvalue problem" << std::endl
              << "**************************" << std::endl;

    invMat.copy(ovlMat);
    // apply the identity shift
    for (int i = 1; i < N; i++)
    {
      hamMat(i, i) += bestShift_i;
      if (invMat(i, i) == 0)
        invMat(i, i) = bestShift_i * bestShift_s;
    }

    // apply the overlap shift
    for (int i = 1; i < N; i++)
      for (int j = 1; j < N; j++)
        hamMat(i, j) += bestShift_s * ovlMat(i, j);

    RealType lowestEV;
    // compute the lowest eigenvalue and the corresponding eigenvector
    if (eigensolver_ == "general")
    {
      app_log() << "  Using generalized eigenvalue solver (ggev)" << std::endl;
      lowestEV = LinearMethod::getLowestEigenvector_Gen(hamMat, invMat, parameterDirections);
    }
    else if (eigensolver_ == "inverse")
    {
      app_log() << "  Using inverse + regular eigenvalue solver (geev)" << std::endl;
      lowestEV = LinearMethod::getLowestEigenvector_Inv(hamMat, invMat, parameterDirections);
    }
    else if (eigensolver_ == "arpack")
    {
      app_log() << "ARPACK not compiled into this QMCPACK executable" << std::endl;
      throw std::runtime_error("ARPACK not present (QMC_USE_ARPACK not set)");
    }
    else
    {
      throw std::runtime_error("Unknown eigenvalue solver: " + eigensolver_);
    }

    app_log() << "  Execution time (eigenvalue) = " << std::setprecision(4) << t_eigen.elapsed() << std::endl;

    // compute the scaling constant to apply to the update
    auto lambda = LinearMethod::getNonLinearRescale(parameterDirections, ovlMat, *optTarget);

    if (do_output_matrices_hdf_)
    {
      hout.write(lowestEV, "lowest_eigenvalue");
      hout.write(parameterDirections, "scaled_eigenvector");
      hout.write(lambda, "non_linear_rescale");
      hout.close();
    }

    // scale the update by the scaling constant
    for (int i = 0; i < numParams; i++)
      parameterDirections.at(i + 1) *= lambda;
  }
  myComm->bcast(parameterDirections);

  // now that we are done building the matrices, prevent further computation of derivative vectors
  optTarget->setneedGrads(false);

  // prepare to use the middle shift's update as the guiding function for a new sample
  if (!freeze_parameters_)
  {
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currentParameters.at(i) + parameterDirections.at(i + 1);
  }

  RealType largestChange(0);
  int max_element = 0;
  for (int i = 0; i < numParams; i++)
    if (std::abs(parameterDirections.at(i + 1)) > largestChange)
    {
      largestChange = std::abs(parameterDirections.at(i + 1));
      max_element   = i;
    }
  app_log() << std::endl
            << "Among totally " << numParams << " optimized parameters, "
            << "largest LM parameter change : " << largestChange << " at parameter " << max_element << std::endl;

  // compute the new cost
  auto effective_weight  = optTarget->correlatedSampling(false);
  const RealType newCost = optTarget->computedCost();


  app_log() << std::endl
            << "******************************************************************************" << std::endl
            << "Init Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << initCost
            << "    New Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << newCost
            << "  Delta Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4)
            << newCost - initCost << std::endl
            << "******************************************************************************" << std::endl;

  if (!optTarget->isEffectiveWeightValid(effective_weight) || qmcplusplus::isnan(newCost))
  {
    app_log() << std::endl << "The new set of parameters is not valid. Revert to the old set!" << std::endl;
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currentParameters.at(i);
    bestShift_s = bestShift_s * shift_s_base;
    if (accept_history[0] == true && accept_history[1] == false) // rejected the one before last and accepted the last
    {
      shift_s_base = std::sqrt(shift_s_base);
      app_log() << "Update shift_s_base to " << shift_s_base << std::endl;
    }
    accept_history <<= 1;
  }
  else
  {
    if (bestShift_s > 1.0e-2)
      bestShift_s = bestShift_s / shift_s_base;
    // say what we are doing
    app_log() << std::endl << "The new set of parameters is valid. Updating the trial wave function!" << std::endl;
    accept_history <<= 1;
    accept_history.set(0, true);
  }

  app_log() << std::endl
            << "*****************************************************************************" << std::endl
            << "Applying the update for shift_i = " << std::scientific << std::right << std::setw(12)
            << std::setprecision(4) << bestShift_i << "     and shift_s = " << std::scientific << std::right
            << std::setw(12) << std::setprecision(4) << bestShift_s << std::endl
            << "*****************************************************************************" << std::endl;

  // perform some finishing touches for this linear method iteration
  finish();


}

void QMCFixedSampleLinearOptimizeBatched::stochastic_reconfiguration_conjugate_gradient()
{
  app_log() << std::endl
            << "*****************************************************************************" << std::endl
            << "                   Running Stochastic Reconfiguration                   " << std::endl
            << "*****************************************************************************" << std::endl;
  // ensure the cost function is set to compute derivative vectors
  optTarget->setneedGrads(true);

  // generate samples and compute weights, local energies, and derivative vectors
  // Note: this has a switch for checkConfigurations or checkConfigurationsSR to do stochastic reconfiguration
  // The SR version avoids calculating the dhpsioverpsi terms and only does dlogpsi
  start();

  // get number of optimizable parameters
  const int numParams = optTarget->getNumParams();

  // get dimension of the linear method matrices
  const int N = numParams + 1;

  // prepare vectors to hold the initial and current parameters
  std::vector<RealType> currentParameters(numParams, 0.0);

  // initialize the initial and current parameter vectors
  for (int i = 0; i < numParams; i++)
    currentParameters.at(i) = std::real(optTarget->Params(i));

  // prepare vectors to hold the parameter update directions for each shift
  std::vector<RealType> parameterDirections;
  parameterDirections.assign(N, 0.0);

  // compute the initial cost
  const RealType initCost = optTarget->computedCost();

  std::vector<RealType> ham(N, 0);

  // for outputing matrices and eigenvalue/vectors to disk
  hdf_archive hout;

  {
    ScopedTimer local(build_olv_ham_timer_);
    Timer t_build_matrices;
    // say what we are doing
    app_log() << std::endl
              << "********************************************************" << std::endl
              << "Building <Psi_i/Psi_0 Psi_j/Psi_0> and <Psi_i/Psi_0 E_L>" << std::endl
              << "********************************************************" << std::endl;

    //This constructs \langle \psi_i/\Psi_0 * E_L \rangle
    optTarget->fillHamVec(ham);

    {
      ScopedTimer local(sr_solver_timer_);
      Timer t_eigen;
      app_log() << std::endl
                << "*********************" << std::endl
                << "Solving linear system" << std::endl
                << "*********************" << std::endl;


      std::vector<RealType> param_update;
      std::vector<RealType> bvec(numParams, 0);
      for (int i = 0; i < numParams; i++)
        bvec[i] = -sr_tau * ham[i + 1];
      ConjugateGradient cg(sr_tolerance, sr_regularization);
      int iterations = cg.run(*optTarget, bvec, param_update);
      app_log() << "Solved iterative krylov in " << iterations << " iterations" << std::endl;
      // compute the scaling constant to apply to the update
      for (int i = 0; i < numParams; i++)
        parameterDirections[i + 1] = param_update[i];
      nrc_opt_.Lambda = cg.getNonLinearRescale(*optTarget);
    }
  }

  //We get the parameter direction from the SR solve above using CG algorithm
  //Then, we can either use a line search with correlated sampling to find the best update along that direction,
  //or we can use a simple approach where we just accept the move based on the size of the step...sr_tau in this case.
  //The line search with correlated sampling converges faster, but can have issues if the weight from correlated
  //sampling gets small and stays small. Otherwise, just taking a small sr_tau will work, but can take a lot of iterations
  //
  //im sure there are better ways to do this
  if (use_line_search_)
  {
    optTarget->setneedGrads(false);

    optdir.resize(numParams, 0);
    optparam.resize(numParams, 0);

    auto costfunc_evaluator = [this](RealType dl) {
      for (int i = 0; i < optparam.size(); i++)
        optTarget->Params(i) = optparam[i] + dl * optdir[i];
      auto effective_weight = optTarget->correlatedSampling(false);
      nrc_opt_.validFuncVal = optTarget->isEffectiveWeightValid(effective_weight);
      return optTarget->computedCost();
    };

    //set up line search stuff
    for (int i = 0; i < numParams; i++)
      optparam[i] = currentParameters[i];
    for (int i = 0; i < numParams; i++)
      optdir[i] = parameterDirections[i + 1];

    RealType bigVec(0);
    for (int i = 0; i < numParams; i++)
      bigVec = std::max(bigVec, std::abs(parameterDirections[i + 1]));

    //Settings for line search, taken from previous_linear_methods_run
    nrc_opt_.TOL              = param_tol / bigVec;
    nrc_opt_.AbsFuncTol       = true;
    nrc_opt_.largeQuarticStep = bigChange / bigVec;
    nrc_opt_.LambdaMax        = 0.5 * nrc_opt_.Lambda;
    bool Valid                = true;
    {
      ScopedTimer local(line_min_timer_);
      Valid = nrc_opt_.lineoptimization2(costfunc_evaluator);
    }

    if (Valid || (!Valid && std::abs(nrc_opt_.Lambda) > 0.0))
    {
      for (int i = 0; i < numParams; i++)
        optTarget->Params(i) = optparam[i] + nrc_opt_.Lambda * optdir[i];
    }
    else
    {
      for (int i = 0; i < numParams; i++)
        optTarget->Params(i) = currentParameters.at(i) + nrc_opt_.Lambda * parameterDirections.at(i + 1);
    }
  }
  else
  {
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currentParameters.at(i) + nrc_opt_.Lambda * parameterDirections.at(i + 1);
  }

  // say what we are doing
  app_log() << std::endl << "The new set of parameters is valid. Updating the trial wave function!" << std::endl;
  accept_history <<= 1;
  accept_history.set(0, true);

  app_log() << std::endl
            << "*****************************************************************************" << std::endl
            << "Applying the update for shift_i = " << std::scientific << std::right << std::setw(12)
            << std::setprecision(4) << bestShift_i << "     and shift_s = " << std::scientific << std::right
            << std::setw(12) << std::setprecision(4) << bestShift_s << std::endl
            << "*****************************************************************************" << std::endl;

  // perform some finishing touches for this linear method iteration
  finish();

  // return whether the cost function's report counter is positive

}

//Function for optimizing using gradient descent
void QMCFixedSampleLinearOptimizeBatched::descent_run()
{
  descent_start();

  int descent_num = descentEngineObj->getDescentNum();

  if (descent_num == 0)
    descentEngineObj->setupUpdate(optTarget->getOptVariables());

  //Store the derivatives and then compute parameter updates
  descentEngineObj->storeDerivRecord();

  descentEngineObj->updateParameters();

  std::vector<ValueType> results = descentEngineObj->retrieveNewParams();


  for (int i = 0; i < results.size(); i++)
  {
    optTarget->Params(i) = std::real(results[i]);
  }

  finish();

}

} // namespace qmcplusplus
