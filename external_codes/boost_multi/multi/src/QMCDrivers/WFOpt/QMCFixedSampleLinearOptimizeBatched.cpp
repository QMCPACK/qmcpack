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
#include "Concurrency/Info.hpp"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "EstimatorInputDelegates.h"
#include "Message/UniformCommunicateError.h"
#include <cassert>
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/random.h"
#include "formic/utils/lmyengine/var_dependencies.h"
#endif
#include <iostream>
#include <fstream>
#include <stdexcept>


namespace qmcplusplus
{
using MatrixOperators::product;


QMCFixedSampleLinearOptimizeBatched::QMCFixedSampleLinearOptimizeBatched(
    const ProjectData& project_data,
    QMCDriverInput&& qmcdriver_input,
    const std::optional<EstimatorManagerInput>& global_emi,
    VMCDriverInput&& vmcdriver_input,
    WalkerConfigurations& wc,
    MCPopulation&& population,
    SampleStack& samples,
    Communicate* comm)
    : QMCDriverNew(
          project_data,
          std::move(qmcdriver_input),
          std::
              nullopt, // this class is not a real QMCDriverNew as far as I can tell so we don't give it the actual global_emi_
          wc,
          std::move(population),
          "QMCLinearOptimizeBatched::",
          comm,
          "QMCLinearOptimizeBatched"),
      objFuncWrapper_(*this),
#ifdef HAVE_LMY_ENGINE
      vdeps(1, std::vector<double>()),
#endif
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
      MinMethod("OneShiftOnly"),
      do_output_matrices_csv_(false),
      do_output_matrices_hdf_(false),
      output_matrices_initialized_(false),
      freeze_parameters_(false),
      generate_samples_timer_(createGlobalTimer("QMCLinearOptimizeBatched::GenerateSamples", timer_level_medium)),
      initialize_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Initialize", timer_level_medium)),
      eigenvalue_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Eigenvalue", timer_level_medium)),
      involvmat_timer_(createGlobalTimer("QMCLinearOptimizedBatched::invert_matrix", timer_level_medium)),
      line_min_timer_(createGlobalTimer("QMCLinearOptimizeBatched::Line_Minimization", timer_level_medium)),
      cost_function_timer_(createGlobalTimer("QMCLinearOptimizeBatched::CostFunction", timer_level_medium)),
      wfNode(NULL),
      vmcdriver_input_(vmcdriver_input),
      samples_(samples),
      global_emi_(global_emi)
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
  // options_LMY_
  m_param.add(options_LMY_.targetExcited, "options_LMY_.targetExcited");
  m_param.add(options_LMY_.block_lm, "options_LMY_.block_lm");
  m_param.add(options_LMY_.nblocks, "options_LMY_.nblocks");
  m_param.add(options_LMY_.nolds, "options_LMY_.nolds");
  m_param.add(options_LMY_.nkept, "options_LMY_.nkept");
  m_param.add(options_LMY_.nsamp_comp, "options_LMY_.nsamp_comp");
  m_param.add(options_LMY_.omega_shift, "omega");
  m_param.add(options_LMY_.max_relative_cost_change, "options_LMY_.max_relative_cost_change");
  m_param.add(options_LMY_.max_param_change, "options_LMY_.max_param_change");
  m_param.add(options_LMY_.num_shifts, "options_LMY_.num_shifts");
  m_param.add(options_LMY_.cost_increase_tol, "options_LMY_.cost_increase_tol");
  m_param.add(options_LMY_.target_shift_i, "options_LMY_.target_shift_i");
  m_param.add(options_LMY_.filter_param, "filter_param");
  m_param.add(options_LMY_.ratio_threshold, "deriv_threshold");
  m_param.add(options_LMY_.store_samples, "store_samples");
  m_param.add(options_LMY_.filter_info, "filter_info");


#ifdef HAVE_LMY_ENGINE
  //app_log() << "construct QMCFixedSampleLinearOptimizeBatched" << endl;
  std::vector<double> shift_scales(3, 1.0);
  EngineObj = new cqmc::engine::LMYEngine<ValueType>(&vdeps,
                                                     false, // exact sampling
                                                     true,  // ground state?
                                                     false, // variance correct,
                                                     true,
                                                     true,  // print matrices,
                                                     true,  // build matrices
                                                     false, // spam
                                                     false, // use var deps?
                                                     true,  // chase lowest
                                                     false, // chase closest
                                                     false, // eom
                                                     false,
                                                     false,  // eom related
                                                     false,  // eom related
                                                     false,  // use block?
                                                     120000, // number of samples
                                                     0,      // number of parameters
                                                     60,     // max krylov iter
                                                     0,      // max spam inner iter
                                                     1,      // spam appro degree
                                                     0,      // eom related
                                                     0,      // eom related
                                                     0,      // eom related
                                                     0.0,    // omega
                                                     0.0,    // var weight
                                                     1.0e-6, // convergence threshold
                                                     0.99,   // minimum S singular val
                                                     0.0, 0.0,
                                                     10.0, // max change allowed
                                                     1.00, // identity shift
                                                     1.00, // overlap shift
                                                     0.3,  // max parameter change
                                                     shift_scales, app_log());
#endif
}

/** Clean up the vector */
QMCFixedSampleLinearOptimizeBatched::~QMCFixedSampleLinearOptimizeBatched()
{
#ifdef HAVE_LMY_ENGINE
  delete EngineObj;
#endif
}

QMCFixedSampleLinearOptimizeBatched::RealType QMCFixedSampleLinearOptimizeBatched::costFunc(RealType dl)
{
  for (int i = 0; i < optparam.size(); i++)
    optTarget->Params(i) = optparam[i] + dl * optdir[i];
  QMCFixedSampleLinearOptimizeBatched::RealType c = optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
  //    if (validFuncVal)
  objFuncWrapper_.validFuncVal = optTarget->IsValid;
  return c;
}

void QMCFixedSampleLinearOptimizeBatched::start()
{
  //close files automatically generated by QMCDriver
  //     branchEngine->finalize();
  //generate samples
  generate_samples_timer_.start();
  generateSamples();
  generate_samples_timer_.stop();
  //store active number of walkers
  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  //reset the rootname
  optTarget->setRootName(get_root_name());
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "   Reading configurations from h5FileRoot " << std::endl;
  //get configuration from the previous run
  Timer t1;
  initialize_timer_.start();
  optTarget->getConfigurations("");
  optTarget->setRng(vmcEngine->getRngRefs());
  NullEngineHandle handle;
  optTarget->checkConfigurations(handle);
  initialize_timer_.stop();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << R"(<opt stage="main" walkers=")" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}

#ifdef HAVE_LMY_ENGINE
void QMCFixedSampleLinearOptimizeBatched::engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                                                       DescentEngine& descentEngineObj,
                                                       std::string MinMethod)
{
  app_log() << "entering engine_start function" << std::endl;

  std::unique_ptr<EngineHandle> handle;
  if (MinMethod == "descent")
    handle = std::make_unique<DescentEngineHandle>(descentEngineObj);
  else if (MinMethod == "adaptive")
    handle = std::make_unique<LMYEngineHandle>(*EngineObj);
  else
    handle = std::make_unique<NullEngineHandle>();


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
  optTarget->setRng(vmcEngine->getRngRefs());
  optTarget->checkConfigurations(*handle);

  initialize_timer_.stop();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << R"(<opt stage="main" walkers=")" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}
#endif


void QMCFixedSampleLinearOptimizeBatched::finish()
{
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;

  if (optTarget->reportH5)
    optTarget->reportParametersH5();
  optTarget->reportParameters();

  app_log() << "</opt>" << std::endl;
  app_log() << "</optimization-report>" << std::endl;
}

void QMCFixedSampleLinearOptimizeBatched::generateSamples()
{
  app_log() << "<optimization-report>" << std::endl;
  t1.restart();
  //     W.reset();
  samples_.resetSampleCount();

  vmcEngine->run();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "</vmc>" << std::endl;
  h5_file_root_ = get_root_name();
}

bool QMCFixedSampleLinearOptimizeBatched::run()
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
    return test_run();
  }
#ifdef HAVE_LMY_ENGINE
  if (options_LMY_.doHybrid)
  {
    app_log() << "Doing hybrid run" << std::endl;
    return hybrid_run();
  }

  // if requested, perform the update via the adaptive three-shift or single-shift method
  if (options_LMY_.current_optimizer_type == OptimizerType::ADAPTIVE)
    return adaptive_three_shift_run();

  if (options_LMY_.current_optimizer_type == OptimizerType::DESCENT)
    return descent_run();

#endif

  if (options_LMY_.current_optimizer_type == OptimizerType::ONESHIFTONLY)
    return one_shift_run();

  return previous_linear_methods_run();
}

bool QMCFixedSampleLinearOptimizeBatched::test_run()
{
  // generate samples and compute weights, local energies, and derivative vectors
  start();

  testEngineObj->run(*optTarget, get_root_name());

  finish();

  return true;
}

bool QMCFixedSampleLinearOptimizeBatched::previous_linear_methods_run()
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
    RealType lastCost(optTarget->Cost(true));
    cost_function_timer_.stop();
    //     if cost function is currently invalid continue
    Valid = optTarget->IsValid;
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
        getLowestEigenvector(Right, currentParameterDirections);
        objFuncWrapper_.Lambda = getNonLinearRescale(currentParameterDirections, S, *optTarget);
      }
      //       biggest gradient in the parameter direction vector
      RealType bigVec(0);
      for (int i = 0; i < numParams; i++)
        bigVec = std::max(bigVec, std::abs(currentParameterDirections[i + 1]));
      //       this can be overwritten during the line minimization
      RealType evaluated_cost(startCost);
      if (MinMethod == "rescale")
      {
        if (std::abs(objFuncWrapper_.Lambda * bigVec) > bigChange)
        {
          goodStep = false;
          app_log() << "  Failed Step. Magnitude of largest parameter change: "
                    << std::abs(objFuncWrapper_.Lambda * bigVec) << std::endl;
          if (stability == 0)
          {
            failedTries++;
            stability--;
          }
          else
            stability = nstabilizers;
        }
        for (int i = 0; i < numParams; i++)
          optTarget->Params(i) = currentParameters[i] + objFuncWrapper_.Lambda * currentParameterDirections[i + 1];
        optTarget->IsValid = true;
      }
      else
      {
        for (int i = 0; i < numParams; i++)
          optparam[i] = currentParameters[i];
        for (int i = 0; i < numParams; i++)
          optdir[i] = currentParameterDirections[i + 1];
        objFuncWrapper_.TOL              = param_tol / bigVec;
        objFuncWrapper_.AbsFuncTol       = true;
        objFuncWrapper_.largeQuarticStep = bigChange / bigVec;
        objFuncWrapper_.LambdaMax        = 0.5 * objFuncWrapper_.Lambda;
        line_min_timer_.start();
        if (MinMethod == "quartic")
        {
          int npts(7);
          objFuncWrapper_.quadstep         = objFuncWrapper_.stepsize * objFuncWrapper_.Lambda;
          objFuncWrapper_.largeQuarticStep = bigChange / bigVec;
          Valid                            = objFuncWrapper_.lineoptimization3(npts, evaluated_cost);
        }
        else
          Valid = objFuncWrapper_.lineoptimization2();
        line_min_timer_.stop();
        RealType biggestParameterChange = bigVec * std::abs(objFuncWrapper_.Lambda);
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
            optTarget->Params(i) = optparam[i] + objFuncWrapper_.Lambda * optdir[i];
          app_log() << "  Good Step. Largest LM parameter change:" << biggestParameterChange << std::endl;
        }
      }

      if (goodStep)
      {
        // 	this may have been evaluated already
        // 	newCost=evaluated_cost;
        //get cost at new minimum
        newCost = optTarget->Cost(false);
        app_log() << " OldCost: " << lastCost << " NewCost: " << newCost << " Delta Cost:" << (newCost - lastCost)
                  << std::endl;
        optTarget->printEstimates();
        //                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
        // mmorales
        Valid = optTarget->IsValid;
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
  return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction
* optimization.
* @param q current xmlNode
* @return true if successful
*/
void QMCFixedSampleLinearOptimizeBatched::process(xmlNodePtr q)
{
  std::string useGPU("yes");
  std::string vmcMove("pbyp");
  std::string ReportToH5("no");
  std::string OutputMatrices("no");
  std::string OutputMatricesHDF("no");
  std::string FreezeParameters("no");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU, "gpu");
  oAttrib.add(vmcMove, "move");
  oAttrib.add(ReportToH5, "hdf5");

  m_param.add(OutputMatrices, "output_matrices_csv", {"no", "yes"});
  m_param.add(OutputMatricesHDF, "output_matrices_hdf", {"no", "yes"});
  m_param.add(FreezeParameters, "freeze_parameters", {"no", "yes"});

  oAttrib.put(q);
  m_param.put(q);

  do_output_matrices_csv_ = (OutputMatrices == "yes");
  do_output_matrices_hdf_ = (OutputMatricesHDF == "yes");
  freeze_parameters_      = (FreezeParameters == "yes");

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


  options_LMY_.doHybrid = false;
  if (MinMethod == "hybrid")
  {
    options_LMY_.doHybrid = true;
    if (!hybridEngineObj)
      hybridEngineObj = std::make_unique<HybridEngine>(myComm, q);

    hybridEngineObj->incrementStepCounter();

    processOptXML(hybridEngineObj->getSelectedXML(), vmcMove, ReportToH5 == "yes", useGPU == "yes");
  }
  else
  {
    processOptXML(q, vmcMove, ReportToH5 == "yes", useGPU == "yes");
  }
}

bool QMCFixedSampleLinearOptimizeBatched::processOptXML(xmlNodePtr opt_xml,
                                                        const std::string& vmcMove,
                                                        bool reportH5,
                                                        bool useGPU)
{
  m_param.put(opt_xml);

  auto iter = OptimizerNames.find(MinMethod);
  if (iter == OptimizerNames.end())
    throw std::runtime_error("Unknown MinMethod!\n");
  options_LMY_.previous_optimizer_type = options_LMY_.current_optimizer_type;
  options_LMY_.current_optimizer_type  = OptimizerNames.at(MinMethod);

  if (options_LMY_.current_optimizer_type == OptimizerType::DESCENT && !descentEngineObj)
    descentEngineObj = std::make_unique<DescentEngine>(myComm, opt_xml);

  // sanity check
  if (options_LMY_.targetExcited && options_LMY_.current_optimizer_type != OptimizerType::ADAPTIVE &&
      options_LMY_.current_optimizer_type != OptimizerType::DESCENT)
    APP_ABORT("options_LMY_.targetExcited = \"yes\" requires that MinMethod = \"adaptive or descent");

#ifdef _OPENMP
  if (options_LMY_.current_optimizer_type == OptimizerType::ADAPTIVE && (omp_get_max_threads() > 1))
  {
    // throw std::runtime_error("OpenMP threading not enabled with AdaptiveThreeShift optimizer. Use MPI for parallelism instead, and set OMP_NUM_THREADS to 1.");
    app_log() << "test version of OpenMP threading with AdaptiveThreeShift optimizer" << std::endl;
  }
#endif

  // check parameter change sanity
  if (options_LMY_.max_param_change <= 0.0)
    throw std::runtime_error(
        "options_LMY_.max_param_change must be positive in QMCFixedSampleLinearOptimizeBatched::put");

  // check cost change sanity
  if (options_LMY_.max_relative_cost_change <= 0.0)
    throw std::runtime_error(
        "options_LMY_.max_relative_cost_change must be positive in QMCFixedSampleLinearOptimizeBatched::put");

  // check shift sanity
  if (shift_i_input <= 0.0)
    throw std::runtime_error("shift_i must be positive in QMCFixedSampleLinearOptimizeBatched::put");
  if (shift_s_input <= 0.0)
    throw std::runtime_error("shift_s must be positive in QMCFixedSampleLinearOptimizeBatched::put");

  // check cost increase tolerance sanity
  if (options_LMY_.cost_increase_tol < 0.0)
    throw std::runtime_error(
        "options_LMY_.cost_increase_tol must be non-negative in QMCFixedSampleLinearOptimizeBatched::put");

  // if this is the first time this function has been called, set the initial shifts
  if (bestShift_i < 0.0 && (options_LMY_.current_optimizer_type == OptimizerType::ADAPTIVE || options_LMY_.doHybrid))
    bestShift_i = shift_i_input;
  if (options_LMY_.current_optimizer_type == OptimizerType::ONESHIFTONLY)
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
  //Overwriting input information is also done here to account for the hybrid method
  QMCDriverInput qmcdriver_input_copy = qmcdriver_input_;
  VMCDriverInput vmcdriver_input_copy = vmcdriver_input_;

  if (MinMethod == "hybrid")
  {
    qmcdriver_input_copy.readXML(hybridEngineObj->getSelectedXML());
    vmcdriver_input_copy.readXML(hybridEngineObj->getSelectedXML());
  }
  else
  {
    qmcdriver_input_copy.readXML(opt_xml);
    vmcdriver_input_copy.readXML(opt_xml);
  }


  // create VMC engine
  vmcEngine =
      std::make_unique<VMCBatched>(project_data_, std::move(qmcdriver_input_copy), global_emi_,
                                   std::move(vmcdriver_input_copy), walker_configs_ref_,
                                   MCPopulation(myComm->size(), myComm->rank(), &population_.get_golden_electrons(),
                                                &population_.get_golden_twf(), &population_.get_golden_hamiltonian()),
                                   samples_, myComm);

  vmcEngine->setUpdateMode(vmcMove[0] == 'p');


  bool AppendRun = false;
  vmcEngine->setStatus(get_root_name(), h5_file_root_, AppendRun);
  vmcEngine->process(qsave);

  vmcEngine->enable_sample_collection();

  auto& qmcdriver_input = vmcEngine->getQMCDriverInput();
  QMCDriverNew::AdjustedWalkerCounts awc =
      adjustGlobalWalkerCount(*myComm, walker_configs_ref_.getActiveWalkers(), qmcdriver_input_.get_total_walkers(),
                              qmcdriver_input_.get_walkers_per_rank(), 1.0, qmcdriver_input_.get_num_crowds());


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

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a vector of three shift values centered around the provided shift.
///
/// \param[in]      central_shift  the central shift
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> QMCFixedSampleLinearOptimizeBatched::prepare_shifts(const double central_shift) const
{
  std::vector<double> retval(options_LMY_.num_shifts);

  // check to see whether the number of shifts is odd
  if (options_LMY_.num_shifts % 2 == 0)
    throw std::runtime_error("number of shifts must be odd in QMCFixedSampleLinearOptimizeBatched::prepare_shifts");

  // decide the central shift index
  int central_index = options_LMY_.num_shifts / 2;

  for (int i = 0; i < options_LMY_.num_shifts; i++)
  {
    if (i < central_index)
      retval.at(i) = central_shift / (4.0 * (central_index - i));
    else if (i > central_index)
      retval.at(i) = central_shift * (4.0 * (i - central_index));
    else if (i == central_index)
      retval.at(i) = central_shift;
    //retval.at(i) = central_shift
    //retval.at(0) = central_shift * 4.0;
    //retval.at(1) = central_shift;
    //retval.at(2) = central_shift / 4.0;
  }
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  prints a header for the summary of each shift's result
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimizeBatched::print_cost_summary_header()
{
  app_log() << "   " << std::right << std::setw(12) << "shift_i";
  app_log() << "   " << std::right << std::setw(12) << "shift_s";
  app_log() << "   " << std::right << std::setw(20) << "max param change";
  app_log() << "   " << std::right << std::setw(20) << "cost function value";
  app_log() << std::endl;
  app_log() << "   " << std::right << std::setw(12) << "------------";
  app_log() << "   " << std::right << std::setw(12) << "------------";
  app_log() << "   " << std::right << std::setw(20) << "--------------------";
  app_log() << "   " << std::right << std::setw(20) << "--------------------";
  app_log() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  prints a summary of the computed cost for the given shift
///
/// \param[in]      si             the identity shift
/// \param[in]      ss             the overlap shift
/// \param[in]      mc             the maximum parameter change
/// \param[in]      cv             the cost function value
/// \param[in]      ind            the shift index: -1 (for initial state), 0, 1, or 2
/// \param[in]      bi             index of the best shift
/// \param[in]      gu             flag telling whether it was a good update
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimizeBatched::print_cost_summary(const double si,
                                                             const double ss,
                                                             const RealType mc,
                                                             const RealType cv,
                                                             const int ind,
                                                             const int bi,
                                                             const bool gu)
{
  if (ind >= 0)
  {
    if (gu)
    {
      app_log() << "   " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << si;
      app_log() << "   " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << ss;
      app_log() << "   " << std::scientific << std::right << std::setw(20) << std::setprecision(4) << mc;
      app_log() << "   " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv;
      //app_log() << "   " << std::right << std::setw(12) << ( ind == 0 ? "big shift" : ( ind == 1 ? "medium shift" : "small shift" ) );
    }
    else
    {
      app_log() << "   " << std::right << std::setw(12) << "N/A";
      app_log() << "   " << std::right << std::setw(12) << "N/A";
      app_log() << "   " << std::right << std::setw(20) << "N/A";
      app_log() << "   " << std::right << std::setw(20) << "N/A";
      app_log() << "   " << std::right << std::setw(12) << "bad update";
    }
  }
  else
  {
    app_log() << "   " << std::right << std::setw(12) << "N/A";
    app_log() << "   " << std::right << std::setw(12) << "N/A";
    app_log() << "   " << std::right << std::setw(20) << "N/A";
    app_log() << "   " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv;
    app_log() << "   " << std::right << std::setw(12) << "initial";
  }
  if (ind == bi)
    app_log() << "  <--";
  app_log() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns whether the proposed new cost is the best compared to the others.
///
/// \param[in]      ii             index of the proposed best cost
/// \param[in]      cv             vector of new costs
/// \param[in]      sh             vector of identity shifts (shift_i values)
/// \param[in]      ic             the initial cost
///
///////////////////////////////////////////////////////////////////////////////////////////////////
bool QMCFixedSampleLinearOptimizeBatched::is_best_cost(const int ii,
                                                       const std::vector<RealType>& cv,
                                                       const std::vector<double>& sh,
                                                       const RealType ic) const
{
  //app_log() << "determining best cost with options_LMY_.cost_increase_tol = " << options_LMY_.cost_increase_tol << " and options_LMY_.target_shift_i = " << options_LMY_.target_shift_i << std::endl;

  // initialize return value
  bool retval = true;

  //app_log() << "retval = " << retval << std::endl;

  // compare to other costs
  for (int i = 0; i < cv.size(); i++)
  {
    // don't compare to yourself
    if (i == ii)
      continue;

    // we only worry about the other value if it is within the maximum relative change threshold and not too high
    const bool other_is_valid =
        ((ic == 0.0 ? 0.0 : std::abs((cv.at(i) - ic) / ic)) < options_LMY_.max_relative_cost_change &&
         cv.at(i) < ic + options_LMY_.cost_increase_tol);
    if (other_is_valid)
    {
      // if we are using a target shift and the cost is not too much higher, then prefer this cost if its shift is closer to the target shift
      if (options_LMY_.target_shift_i > 0.0)
      {
        const bool closer_to_target =
            (std::abs(sh.at(ii) - options_LMY_.target_shift_i) < std::abs(sh.at(i) - options_LMY_.target_shift_i));
        const bool cost_is_similar    = (std::abs(cv.at(ii) - cv.at(i)) < options_LMY_.cost_increase_tol);
        const bool cost_is_much_lower = (!cost_is_similar && cv.at(ii) < cv.at(i) - options_LMY_.cost_increase_tol);
        if (cost_is_much_lower || (closer_to_target && cost_is_similar))
          retval = (retval && true);
        else
          retval = false;

        // if we are not using a target shift, then prefer this cost if it is lower
      }
      else
      {
        retval = (retval && cv.at(ii) <= cv.at(i));
      }
    }

    //app_log() << "cv.at(ii)   = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(ii) << " <= "
    //          << "cv.at(i)    = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(i)  << " ?" << std::endl;
    //app_log() << "retval = " << retval << std::endl;
  }

  // new cost can only be the best cost if it is less than (or not too much higher than) the initial cost
  retval = (retval && cv.at(ii) < ic + options_LMY_.cost_increase_tol);
  //app_log() << "cv.at(ii)   = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(ii) << " <= "
  //          << "ic          = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << ic        << " ?" << std::endl;
  //app_log() << "retval = " << retval << std::endl;

  // new cost is only best if it's relative change from the initial cost is not too large ( or if the initial cost is exactly zero )
  retval = (retval && (ic == 0.0 ? 0.0 : std::abs((cv.at(ii) - ic) / ic)) < options_LMY_.max_relative_cost_change);
  //app_log() << "std::abs( ( cv.at(ii) - ic ) / ic ) = " << std::fixed << std::right << std::setw(20) << std::setprecision(12)
  //          << std::abs( ( cv.at(ii) - ic ) / ic ) << " <= " << this->options_LMY_.max_relative_cost_change << " ? " << std::endl;
  //app_log() << "retval = " << retval << std::endl;
  //app_log() << std::endl;

  // return whether the proposed cost is actually the best
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  For each set of shifts, solves the linear method eigenproblem by building and
///         diagonalizing the matrices.
///
/// \param[in]      shfits_i              vector of identity shifts
/// \param[in]      shfits_s              vector of overlap shifts
/// \param[out]     parameterDirections   on exit, the update directions for the different shifts
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimizeBatched::solveShiftsWithoutLMYEngine(
    const std::vector<double>& shifts_i,
    const std::vector<double>& shifts_s,
    std::vector<std::vector<RealType>>& parameterDirections)
{
  // get number of shifts to solve
  const int nshifts = shifts_i.size();

  // get number of optimizable parameters
  const int numParams = optTarget->getNumParams();

  // get dimension of the linear method matrices
  const int N = numParams + 1;

  // prepare vectors to hold the parameter updates
  parameterDirections.resize(nshifts);
  for (int i = 0; i < parameterDirections.size(); i++)
    parameterDirections.at(i).assign(N, 0.0);

  // allocate the matrices we will need
  Matrix<RealType> ovlMat(N, N);
  ovlMat = 0.0;
  Matrix<RealType> hamMat(N, N);
  hamMat = 0.0;
  Matrix<RealType> invMat(N, N);
  invMat = 0.0;
  Matrix<RealType> sftMat(N, N);
  sftMat = 0.0;
  Matrix<RealType> prdMat(N, N);
  prdMat = 0.0;

  // build the overlap and hamiltonian matrices
  optTarget->fillOverlapHamiltonianMatrices(hamMat, ovlMat);

  //// print the hamiltonian matrix
  //app_log() << std::endl;
  //app_log() << "printing H matrix:" << std::endl;
  //for (int i = 0; i < hamMat.rows(); i++) {
  //  for (int j = 0; j < hamMat.cols(); j++)
  //    app_log() << " " << std::scientific << std::right << std::setw(14) << std::setprecision(5) << hamMat(i,j);
  //  app_log() << std::endl;
  //}
  //app_log() << std::endl;

  //// print the overlap matrix
  //app_log() << std::endl;
  //app_log() << "printing S matrix:" << std::endl;
  //for (int i = 0; i < ovlMat.rows(); i++) {
  //  for (int j = 0; j < ovlMat.cols(); j++)
  //    app_log() << " " << std::scientific << std::right << std::setw(14) << std::setprecision(5) << ovlMat(i,j);
  //  app_log() << std::endl;
  //}
  //app_log() << std::endl;

  // compute the inverse of the overlap matrix
  invMat.copy(ovlMat);
  invert_matrix(invMat, false);

  // compute the update for each shift
  for (int shift_index = 0; shift_index < nshifts; shift_index++)
  {
    // prepare to shift the hamiltonain matrix
    sftMat.copy(hamMat);

    // apply the identity shift
    for (int i = 1; i < N; i++)
      sftMat(i, i) += shifts_i.at(shift_index);

    // apply the overlap shift
    for (int i = 1; i < N; i++)
      for (int j = 1; j < N; j++)
        sftMat(i, j) += shifts_s.at(shift_index) * ovlMat(i, j);

    // multiply the shifted hamiltonian matrix by the inverse of the overlap matrix
    qmcplusplus::MatrixOperators::product(invMat, sftMat, prdMat);

    // transpose the result (why?)
    for (int i = 0; i < N; i++)
      for (int j = i + 1; j < N; j++)
        std::swap(prdMat(i, j), prdMat(j, i));

    // compute the lowest eigenvalue of the product matrix and the corresponding eigenvector
    getLowestEigenvector(prdMat, parameterDirections.at(shift_index));

    // compute the scaling constant to apply to the update
    objFuncWrapper_.Lambda = getNonLinearRescale(parameterDirections.at(shift_index), ovlMat, *optTarget);

    // scale the update by the scaling constant
    for (int i = 0; i < numParams; i++)
      parameterDirections.at(shift_index).at(i + 1) *= objFuncWrapper_.Lambda;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Performs one iteration of the linear method using an adaptive scheme that tries three
///         different shift magnitudes and picks the best one.
///         The scheme is adaptive in that it saves the best shift to use as a starting point
///         in the next iteration.
///         Note that the best shift is chosen based on a different sample than that used to
///         construct the linear method matrices in order to avoid over-optimizing on a particular
///         sample.
///
/// \return  ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LMY_ENGINE
bool QMCFixedSampleLinearOptimizeBatched::adaptive_three_shift_run()
{
  EngineObj->setStoringSamples(options_LMY_.store_samples);

  //Set whether LM will only update a filtered set of parameters
  EngineObj->setFiltering(options_LMY_.filter_param);
  EngineObj->setFilterInfo(options_LMY_.filter_info);

  if (options_LMY_.filter_param && !options_LMY_.store_samples)
    myComm->barrier_and_abort(" Error: Parameter Filtration requires storing the samples. \n");

  if (options_LMY_.filter_param)
    EngineObj->setThreshold(options_LMY_.ratio_threshold);

  // remember what the cost function grads flag was
  const bool saved_grads_flag = optTarget->getneedGrads();

  // remember the initial number of samples
  const int init_num_samp = optTarget->getNumSamples();

  // the index of central shift
  const int central_index = options_LMY_.num_shifts / 2;

  // get number of optimizable parameters
  const int numParams = optTarget->getNumParams();

  // prepare the shifts that we will try
  const std::vector<double> shifts_i = prepare_shifts(bestShift_i);
  const std::vector<double> shifts_s = prepare_shifts(bestShift_s);
  std::vector<double> shift_scales(shifts_i.size(), 1.0);
  for (int i = 0; i < shift_scales.size(); i++)
    shift_scales.at(i) = shifts_i.at(i) / shift_i_input;

  // ensure the cost function is set to compute derivative vectors
  optTarget->setneedGrads(true);

  // prepare previous updates
  int count = 0;
  while (options_LMY_.block_lm && previous_update.size() < options_LMY_.nolds)
  {
    previous_update.push_back(formic::ColVec<double>(numParams));
    for (int i = 0; i < numParams; i++)
      previous_update.at(count).at(i) = 2.0 * (formic::random_number<double>() - 0.5);
    count++;
  }

  if (!EngineObj->full_init())
  {
    // prepare a variable dependency object with no dependencies
    formic::VarDeps real_vdeps(numParams, std::vector<double>());
    vdeps = real_vdeps;
    EngineObj->get_param(&vdeps,
                         false, // exact sampling
                         !options_LMY_.targetExcited,
                         false, // variable deps use?
                         false, // eom
                         false, // ssquare
                         options_LMY_.block_lm, 12000, numParams, options_LMY_.omega_shift,
                         options_LMY_.max_relative_cost_change, shifts_i.at(central_index), shifts_s.at(central_index),
                         options_LMY_.max_param_change, shift_scales);
  }

  //Reset parameter number for vdeps to the total number in case filtration happened on a previous iteration
  if (options_LMY_.filter_param)
  {
    formic::VarDeps tmp_vdeps(numParams, std::vector<double>());
    vdeps = tmp_vdeps;
    EngineObj->var_deps_ptr_update(&vdeps);
  }

  // update shift
  EngineObj->shift_update(shift_scales);

  // turn on wavefunction update mode
  EngineObj->turn_on_update();

  //The initial intialization of the LM engine is handled differently if parameters are being filtered
  if (!options_LMY_.filter_param)
  {
    // initialize the engine if we do not use block lm or it's the first part of block lm
    EngineObj->initialize(options_LMY_.nblocks, 0, options_LMY_.nkept, previous_update, false);
    EngineObj->reset();
  }
  else
  {
    app_log() << "Skipping initialization at first" << std::endl;
    EngineObj->store_blocked_lm_info(options_LMY_.nblocks, options_LMY_.nkept);
  }


  // reset the engine
  EngineObj->reset();

  // generate samples and compute weights, local energies, and derivative vectors
  engine_start(EngineObj, *descentEngineObj, MinMethod);

  int new_num = 0;

  //To handle different cases for the LM's mode of operation, first check if samples are being stored
  if (options_LMY_.store_samples)
  {
    //Need to clear lists from previous iter
    EngineObj->reset();

    //If samples are being stored, check for the subcase where parameters are also being filtered
    if (options_LMY_.filter_param)
    {
      EngineObj->selectParameters();

      for (int i = 0; i < numParams; i++)
        if (EngineObj->getParameterSetting(i))
          new_num++;

      formic::VarDeps real_vdeps(new_num, std::vector<double>());
      vdeps = real_vdeps;
      EngineObj->var_deps_ptr_update(&vdeps);

      //Also need to check if Blocked LM is being used
      if (EngineObj->use_blm())
      {
        //If so, the old update vectors need to be trimmed to remove the filtered out parameters
        std::vector<formic::ColVec<double>> trimmed_old_updates(previous_update.size());

        //Check if this Blocked LM step is part of a hybrid optimization
        if (EngineObj->getOnHybrid())
        {
          //If so, get the old update vectors from the descent engine
          std::vector<std::vector<ValueType>> hybridBLM_Input = descentEngineObj->retrieveHybridBLM_Input();


          app_log() << "Blocked LM is part of hybrid run. Need to filter vectors from descent. " << std::endl;

          //This section handles the trimming of the old update vectors from descent
          for (int i = 0; i < hybridBLM_Input.size(); i++)
          {
            std::vector<ValueType> full_vec = hybridBLM_Input[i];
            std::vector<ValueType> filtered_vec;

            formic::ColVec<double> reduced_vector(new_num, 0.0);
            int count = 0;

            for (int j = 0; j < full_vec.size(); j++)
              if (EngineObj->getParameterSetting(j))
              {
                filtered_vec.push_back(full_vec[j]);
                reduced_vector[count] = formic::real(full_vec[j]);
                count++;
              }

            hybridBLM_Input[i]     = filtered_vec;
            trimmed_old_updates[i] = reduced_vector;
          }

#if !defined(QMC_COMPLEX)
          EngineObj->setHybridBLM_Input(hybridBLM_Input);
#endif

          EngineObj->initialize(options_LMY_.nblocks, 0, options_LMY_.nkept, trimmed_old_updates, false);
          EngineObj->reset();
        }
        //If the Blocked LM is not part of a hybrid run, carry out the trimming of the old updates here
        else
        {
          app_log() << "Regular Blocked LM run. Need to filter old update vectors. " << std::endl;

          for (int i = 0; i < previous_update.size(); i++)
          {
            formic::ColVec<double> full_vec = previous_update[i];

            formic::ColVec<double> reduced_vector(new_num, 0.0);
            int count = 0;

            for (int j = 0; j < full_vec.size(); j++)
              if (EngineObj->getParameterSetting(j))
              {
                reduced_vector[count] = full_vec[j];
                count++;
              }

            trimmed_old_updates[i] = reduced_vector;
          }

          EngineObj->initialize(options_LMY_.nblocks, 0, options_LMY_.nkept, trimmed_old_updates, false);
          EngineObj->reset();
        }
      }
    }
    //If not filtering parameters and only storing samples, can proceed with the rest of the LM engine initialization
    else
    {
      EngineObj->initialize(options_LMY_.nblocks, 0, options_LMY_.nkept, previous_update, false);
      EngineObj->reset();
    }


    //This function call builds the matrices from the stored samples
    EngineObj->buildMatricesFromDerivatives();
  }


  // get dimension of the linear method matrices
  int N = numParams + 1;
  if (options_LMY_.filter_param)
    N = new_num + 1;

  // have the cost function prepare derivative vectors
  EngineObj->energy_target_compute();
  const RealType starting_cost = EngineObj->target_value();
  const RealType init_energy   = EngineObj->energy_mean();

  // print out the initial energy
  app_log() << std::endl
            << "*************************************************************************************************"
            << std::endl
            << "Solving the linear method equations on the initial sample with initial energy" << std::setw(20)
            << std::setprecision(12) << init_energy << std::endl
            << "*************************************************************************************************"
            << std::endl
            << std::endl;

  // prepare wavefunction update which does nothing if we do not use block lm
  EngineObj->wfn_update_prep();

  if (options_LMY_.block_lm)
  {
    if (!options_LMY_.store_samples)
    {
      optTarget->setneedGrads(true);

      int numOptParams = optTarget->getNumParams();

      // reset the engine object
      EngineObj->reset();

      // finish last sample
      finish();

      // take sample
      engine_start(EngineObj, *descentEngineObj, MinMethod);
    }
    else
    {
      EngineObj->clear_histories();
      EngineObj->reset();

      finish();

      if (options_LMY_.filter_param)
      {
        engine_start(EngineObj, *descentEngineObj, MinMethod);
        EngineObj->buildMatricesFromDerivatives();
      }
      else
      {
        engine_start(EngineObj, *descentEngineObj, MinMethod);
        app_log() << "Should be building matrices from stored samples" << std::endl;
        EngineObj->buildMatricesFromDerivatives();
      }
    }
  }

  //Need to wipe the stored samples after they are no longer needed and before the next iteration
  if (options_LMY_.store_samples)
  {
    EngineObj->clear_histories();
  }

  // say what we are doing
  app_log() << std::endl
            << "*********************************************************" << std::endl
            << "Solving the linear method equations on the initial sample" << std::endl
            << "*********************************************************" << std::endl
            << std::endl;

  // for each set of shifts, solve the linear method equations for the parameter update direction
  std::vector<std::vector<RealType>> parameterDirections;
#ifdef HAVE_LMY_ENGINE
  // call the engine to perform update
  EngineObj->wfn_update_compute();
#else
  solveShiftsWithoutLMYEngine(shifts_i, shifts_s, parameterDirections);
#endif

  // size update direction vector correctly
  parameterDirections.resize(shifts_i.size());
  for (int i = 0; i < shifts_i.size(); i++)
  {
    parameterDirections.at(i).assign(N, 0.0);
    if (true)
    {
      for (int j = 0; j < N; j++)
        parameterDirections.at(i).at(j) = std::real(EngineObj->wfn_update().at(i * N + j));
    }
    else
      parameterDirections.at(i).at(0) = 1.0;
  }

  //If paramters are being filtered need to expand the LM updates from the engine to the full parameter set.
  //There will be updates of 0 for parameters that were filtered out before derivative ratios were used by the engine.
  if (options_LMY_.filter_param)
  {
    std::vector<std::vector<RealType>> tmpParameterDirections;
    tmpParameterDirections.resize(shifts_i.size());

    for (int i = 0; i < shifts_i.size(); i++)
    {
      tmpParameterDirections.at(i).assign(numParams + 1, 0.0);
      int lm_update_idx = 0;
      for (int j = 0; j < numParams + 1; j++)
      {
        if (j == 0)
        {
          tmpParameterDirections.at(i).at(j) = parameterDirections.at(i).at(j);
          lm_update_idx++;
        }
        else if (EngineObj->getParameterSetting(j - 1) == true)
        {
          tmpParameterDirections.at(i).at(j) = parameterDirections.at(i).at(lm_update_idx);
          lm_update_idx++;
        }
      }
      parameterDirections.at(i) = tmpParameterDirections.at(i);
    }
  }

  //From this point, the comparison of the 3 diffferent shifts' updates should proceed as normal regardless of the sample storage or parameter filtration settings.

  // now that we are done with them, prevent further computation of derivative vectors
  optTarget->setneedGrads(false);

  // prepare vectors to hold the initial and current parameters
  std::vector<RealType> currParams(numParams, 0.0);

  // initialize the initial and current parameter vectors
  for (int i = 0; i < numParams; i++)
    currParams.at(i) = optTarget->Params(i);

  // create a vector telling which updates are within our constraints
  std::vector<bool> good_update(parameterDirections.size(), true);

  // compute the largest parameter change for each shift, and zero out updates that have too-large changes
  std::vector<RealType> max_change(parameterDirections.size(), 0.0);
  for (int k = 0; k < parameterDirections.size(); k++)
  {
    for (int i = 0; i < numParams; i++)
      max_change.at(k) =
          std::max(max_change.at(k), std::abs(parameterDirections.at(k).at(i + 1) / parameterDirections.at(k).at(0)));
    good_update.at(k) = (good_update.at(k) && max_change.at(k) <= options_LMY_.max_param_change);
  }

  // prepare to use the middle shift's update as the guiding function for a new sample
  for (int i = 0; i < numParams; i++)
    optTarget->Params(i) = currParams.at(i) + parameterDirections.at(central_index).at(i + 1);

  // say what we are doing
  app_log() << std::endl
            << "************************************************************" << std::endl
            << "Updating the guiding function with the middle shift's update" << std::endl
            << "************************************************************" << std::endl
            << std::endl;

  // generate the new sample on which we will compare the different shifts

  finish();
  app_log() << std::endl
            << "*************************************************************" << std::endl
            << "Generating a new sample based on the updated guiding function" << std::endl
            << "*************************************************************" << std::endl
            << std::endl;

  //Apparently the batched drivers are intended to be run only once, which
  //means that the origianl version of adaptive_three_shift will not work as
  //calling start or engine_start at this point will lead to vmcEngine being run again.
  //This will lead to a slight difference in behavior compared to the
  //legacy drivers as those could be run a second time to obtain samples based
  //on the wave function from the middle shift.
  //It is possible this difference may not make much difference in
  //practical optimization performance, but that is unexplored.


  // say what we are doing
  app_log() << std::endl
            << "******************************************************************" << std::endl
            << "Comparing different shifts' cost function values on updated sample" << std::endl
            << "******************************************************************" << std::endl
            << std::endl;

  // update the current parameters to those of the new guiding function
  for (int i = 0; i < numParams; i++)
    currParams.at(i) = optTarget->Params(i);

  // compute cost function for the initial parameters (by subtracting the middle shift's update back off)
  for (int i = 0; i < numParams; i++)
    optTarget->Params(i) = currParams.at(i) - parameterDirections.at(central_index).at(i + 1);
  optTarget->IsValid      = true;
  const RealType initCost = optTarget->LMYEngineCost(false, EngineObj);

  // compute the update directions for the smaller and larger shifts relative to that of the middle shift
  for (int i = 0; i < numParams; i++)
  {
    for (int j = 0; j < parameterDirections.size(); j++)
    {
      if (j != central_index)
        parameterDirections.at(j).at(i + 1) -= parameterDirections.at(central_index).at(i + 1);
    }
  }

  // prepare a vector to hold the cost function value for each different shift
  std::vector<RealType> costValues(options_LMY_.num_shifts, 0.0);

  // compute the cost function value for each shift and make sure the change is within our constraints
  for (int k = 0; k < parameterDirections.size(); k++)
  {
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currParams.at(i) + (k == central_index ? 0.0 : parameterDirections.at(k).at(i + 1));
    optTarget->IsValid = true;
    costValues.at(k)   = optTarget->LMYEngineCost(false, EngineObj);
    good_update.at(k)  = (good_update.at(k) &&
                         std::abs((initCost - costValues.at(k)) / initCost) < options_LMY_.max_relative_cost_change);
    if (!good_update.at(k))
      costValues.at(k) = std::abs(1.5 * initCost) + 1.0;
  }

  // find the best shift and the corresponding update direction
  const std::vector<RealType>* bestDirection = 0;
  int best_shift                             = -1;
  for (int k = 0;
       k < costValues.size() && std::abs((initCost - initCost) / initCost) < options_LMY_.max_relative_cost_change; k++)
    if (is_best_cost(k, costValues, shifts_i, initCost) && good_update.at(k))
    {
      best_shift    = k;
      bestDirection = &parameterDirections.at(k);
    }

  // print the results for each shift
  app_log() << std::endl;
  print_cost_summary_header();
  print_cost_summary(0.0, 0.0, 0.0, initCost, -1, best_shift, true);
  for (int k = 0; k < good_update.size(); k++)
    print_cost_summary(shifts_i.at(k), shifts_s.at(k), max_change.at(k), costValues.at(k), k, best_shift,
                       good_update.at(k));

  // if any of the shifts produced a good update, apply the best such update and remember those shifts for next time
  if (bestDirection)
  {
    bestShift_i = shifts_i.at(best_shift);
    bestShift_s = shifts_s.at(best_shift);
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currParams.at(i) + (best_shift == central_index ? 0.0 : bestDirection->at(i + 1));
    app_log() << std::endl
              << "*****************************************************************************" << std::endl
              << "Applying the update for shift_i = " << std::scientific << std::right << std::setw(12)
              << std::setprecision(4) << bestShift_i << "     and shift_s = " << std::scientific << std::right
              << std::setw(12) << std::setprecision(4) << bestShift_s << std::endl
              << "*****************************************************************************" << std::endl
              << std::endl;

    // otherwise revert to the old parameters and set the next shift to be larger
  }
  else
  {
    bestShift_i *= 10.0;
    bestShift_s *= 10.0;
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currParams.at(i) - parameterDirections.at(central_index).at(i + 1);
    app_log() << std::endl
              << "***********************************************************" << std::endl
              << "Reverting to old parameters and increasing shift magnitudes" << std::endl
              << "***********************************************************" << std::endl
              << std::endl;
  }

  // save the update for future linear method iterations
  if (options_LMY_.block_lm && bestDirection)
  {
    // save the difference between the updated and old variables
    formic::ColVec<RealType> update_dirs(numParams, 0.0);
    for (int i = 0; i < numParams; i++)
      // take the real part since blocked LM currently does not support complex parameter optimization
      update_dirs.at(i) = std::real(bestDirection->at(i + 1) + parameterDirections.at(central_index).at(i + 1));
    previous_update.insert(previous_update.begin(), update_dirs);

    // eliminate the oldest saved update if necessary
    while (previous_update.size() > options_LMY_.nolds)
      previous_update.pop_back();
  }

  // return the cost function grads flag to what it was
  optTarget->setneedGrads(saved_grads_flag);

  // perform some finishing touches for this linear method iteration
  finish();

  // set the number samples to be initial one
  optTarget->setNumSamples(init_num_samp);

  //app_log() << "block first second third end " << options_LMY_.block_first << options_LMY_.block_second << options_LMY_.block_third << endl;
  // return whether the cost function's report counter is positive
  return (optTarget->getReportCounter() > 0);
}
#endif

bool QMCFixedSampleLinearOptimizeBatched::one_shift_run()
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

  // say what we are doing
  app_log() << std::endl
            << "*****************************************" << std::endl
            << "Building overlap and Hamiltonian matrices" << std::endl
            << "*****************************************" << std::endl;

  // allocate the matrices we will need
  Matrix<RealType> ovlMat(N, N);
  ovlMat = 0.0;
  Matrix<RealType> hamMat(N, N);
  hamMat = 0.0;
  Matrix<RealType> invMat(N, N);
  invMat = 0.0;
  Matrix<RealType> prdMat(N, N);
  prdMat = 0.0;

  // build the overlap and hamiltonian matrices
  optTarget->fillOverlapHamiltonianMatrices(hamMat, ovlMat);
  invMat.copy(ovlMat);

  if (do_output_matrices_csv_)
  {
    output_overlap_.output(ovlMat);
    output_hamiltonian_.output(hamMat);
  }

  hdf_archive hout;
  if (do_output_matrices_hdf_)
  {
    std::string newh5 = get_root_name() + ".linear_matrices.h5";
    hout.create(newh5, H5F_ACC_TRUNC);
    hout.write(ovlMat, "overlap");
    hout.write(hamMat, "Hamiltonian");
    hout.write(bestShift_i, "bestShift_i");
    hout.write(bestShift_s, "bestShift_s");
  }

  // apply the identity shift
  for (int i = 1; i < N; i++)
  {
    hamMat(i, i) += bestShift_i;
    if (invMat(i, i) == 0)
      invMat(i, i) = bestShift_i * bestShift_s;
  }

  // compute the inverse of the overlap matrix
  {
    ScopedTimer local(involvmat_timer_);
    invert_matrix(invMat, false);
  }

  // apply the overlap shift
  for (int i = 1; i < N; i++)
    for (int j = 1; j < N; j++)
      hamMat(i, j) += bestShift_s * ovlMat(i, j);

  // multiply the shifted hamiltonian matrix by the inverse of the overlap matrix
  qmcplusplus::MatrixOperators::product(invMat, hamMat, prdMat);

  // transpose the result (why?)
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      std::swap(prdMat(i, j), prdMat(j, i));

  // compute the lowest eigenvalue of the product matrix and the corresponding eigenvector
  RealType lowestEV = 0.;
  {
    ScopedTimer local(eigenvalue_timer_);
    lowestEV = getLowestEigenvector(prdMat, parameterDirections);
  }

  // compute the scaling constant to apply to the update
  objFuncWrapper_.Lambda = getNonLinearRescale(parameterDirections, ovlMat, *optTarget);

  if (do_output_matrices_hdf_)
  {
    hout.write(lowestEV, "lowest_eigenvalue");
    hout.write(parameterDirections, "scaled_eigenvector");
    hout.write(objFuncWrapper_.Lambda, "non_linear_rescale");
    hout.close();
  }

  // scale the update by the scaling constant
  for (int i = 0; i < numParams; i++)
    parameterDirections.at(i + 1) *= objFuncWrapper_.Lambda;

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
  optTarget->IsValid     = true;
  const RealType newCost = optTarget->Cost(false);

  app_log() << std::endl
            << "******************************************************************************" << std::endl
            << "Init Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << initCost
            << "    New Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4) << newCost
            << "  Delta Cost = " << std::scientific << std::right << std::setw(12) << std::setprecision(4)
            << newCost - initCost << std::endl
            << "******************************************************************************" << std::endl;

  if (!optTarget->IsValid || qmcplusplus::isnan(newCost))
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

  // return whether the cost function's report counter is positive
  return (optTarget->getReportCounter() > 0);
}

#ifdef HAVE_LMY_ENGINE
//Function for optimizing using gradient descent
bool QMCFixedSampleLinearOptimizeBatched::descent_run()
{
  //Compute Lagrangian derivatives needed for parameter updates with engine_checkConfigurations, which is called inside engine_start
  engine_start(EngineObj, *descentEngineObj, MinMethod);

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

  //If descent is being run as part of a hybrid optimization, need to check if a vector of
  //parameter differences should be stored.
  if (options_LMY_.doHybrid)
  {
    int store_num = descentEngineObj->retrieveStoreFrequency();
    bool store    = hybridEngineObj->queryStore(store_num, OptimizerType::DESCENT);
    if (store)
    {
      descentEngineObj->storeVectors(results);
    }
  }

  finish();
  return (optTarget->getReportCounter() > 0);
}
#endif


//Function for controlling the alternation between sections of descent optimization and BLM optimization.
#ifdef HAVE_LMY_ENGINE
bool QMCFixedSampleLinearOptimizeBatched::hybrid_run()
{
  app_log() << "This is methodName: " << MinMethod << std::endl;

  //Either the adaptive BLM or descent optimization is run

  //Ensure LM engine knows it is being used as part of a hybrid run
  EngineObj->setOnHybrid(true);

  if (options_LMY_.current_optimizer_type == OptimizerType::ADAPTIVE)
  {
    //If the optimization just switched to using the BLM, need to transfer a set
    //of vectors to the BLM engine.
    if (options_LMY_.previous_optimizer_type == OptimizerType::DESCENT)
    {
      descentEngineObj->resetStorageCount();
      std::vector<std::vector<ValueType>> hybridBLM_Input = descentEngineObj->retrieveHybridBLM_Input();
#if !defined(QMC_COMPLEX)
      //FIXME once complex is fixed in BLM engine
      EngineObj->setHybridBLM_Input(hybridBLM_Input);
#endif
    }
    adaptive_three_shift_run();
  }

  if (options_LMY_.current_optimizer_type == OptimizerType::DESCENT)
    descent_run();

  app_log() << "Finished a hybrid step" << std::endl;
  return (optTarget->getReportCounter() > 0);
}
#endif

} // namespace qmcplusplus
