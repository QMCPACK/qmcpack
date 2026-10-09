//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCFixedSampleLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "RandomNumberControl.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCDrivers/WFOpt/GradientTest.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include "Message/UniformCommunicateError.h"
#include "Numerics/DeterminantOperators.h"
#include "LinearMethod.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <stdexcept>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{
using MatrixOperators::product;


QMCFixedSampleLinearOptimize::QMCFixedSampleLinearOptimize(const ProjectData& project_data,
                                                           MCWalkerConfiguration& w,
                                                           TrialWaveFunction& psi,
                                                           QMCHamiltonian& h,
                                                           Communicate& comm)
    : QMCDriver(project_data, w, psi, h, comm, "QMCFixedSampleLinearOptimize"),
      nstabilizers(3),
      stabilizerScale(2.0),
      bigChange(50),
      exp0(-16),
      stepsize(0.25),
      StabilizerMethod("best"),
      bestShift_i(-1.0),
      bestShift_s(-1.0),
      shift_i_input(0.01),
      shift_s_input(1.00),
      accept_history(3),
      shift_s_base(4.0),
      MinMethod("OneShiftOnly"),
      current_optimizer_type_(OptimizerType::NONE),
      do_output_matrices_(false),
      output_matrices_initialized_(false),
      freeze_parameters_(false),
      Max_iterations(1),
      wfNode(NULL),
      param_tol(1e-4),
      generate_samples_timer_(createGlobalTimer("QMCLinearOptimize::generateSamples", timer_level_medium)),
      initialize_timer_(createGlobalTimer("QMCLinearOptimize::Initialize", timer_level_medium)),
      eigenvalue_timer_(createGlobalTimer("QMCLinearOptimize::EigenvalueSolve", timer_level_medium)),
      involvmat_timer_(createGlobalTimer("QMCLinearOptimize::invertOverlapMat", timer_level_medium)),
      line_min_timer_(createGlobalTimer("QMCLinearOptimize::Line_Minimization", timer_level_medium)),
      cost_function_timer_(createGlobalTimer("QMCLinearOptimize::CostFunction", timer_level_medium))
{
  IsQMCDriver = false;
  //set the optimization flag
  qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  RootName = "pot";
  m_param.add(Max_iterations, "max_its");
  m_param.add(nstabilizers, "nstabilizers");
  m_param.add(stabilizerScale, "stabilizerscale");
  m_param.add(bigChange, "bigchange");
  m_param.add(MinMethod, "MinMethod");
  m_param.add(exp0, "exp0");
  m_param.add(shift_i_input, "shift_i");
  m_param.add(shift_s_input, "shift_s");
  m_param.add(param_tol, "alloweddifference");
}

QMCFixedSampleLinearOptimize::~QMCFixedSampleLinearOptimize() = default;

void QMCFixedSampleLinearOptimize::test_run()
{
  // generate samples and compute weights, local energies, and derivative vectors
  start();

  testEngineObj->run(*optTarget, get_root_name());

  finish();
}

void QMCFixedSampleLinearOptimize::run()
{
  if (do_output_matrices_ && !output_matrices_initialized_)
  {
    size_t numParams = optTarget->getNumParams();
    size_t N         = numParams + 1;
    output_overlap_.init_file(get_root_name(), "ovl", N);
    output_hamiltonian_.init_file(get_root_name(), "ham", N);
    output_matrices_initialized_ = true;
  }

  if (doGradientTest)
  {
    app_log() << "Doing gradient test run" << std::endl;
    test_run();
    return;
  }


  if (current_optimizer_type_ == OptimizerType::ONESHIFTONLY)
  {
    one_shift_run();
    return;
  }

  start();
  bool Valid(true);
  int Total_iterations(0);
  //size of matrix
  size_t numParams = optTarget->getNumParams();
  size_t N         = numParams + 1;
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
    validFuncVal          = optTarget->isEffectiveWeightValid(effective_weight);
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
      {
        ScopedTimer local(involvmat_timer_);
        invert_matrix(RightT, false);
      }
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
        Lambda = LinearMethod::getNonLinearRescale(currentParameterDirections, S, *optTarget);
      }
      //       biggest gradient in the parameter direction vector
      RealType bigVec(0);
      for (int i = 0; i < numParams; i++)
        bigVec = std::max(bigVec, std::abs(currentParameterDirections[i + 1]));
      //       this can be overwritten during the line minimization
      RealType evaluated_cost(startCost);
      if (MinMethod == "rescale")
      {
        if (std::abs(Lambda * bigVec) > bigChange)
        {
          goodStep = false;
          app_log() << "  Failed Step. Magnitude of largest parameter change: " << std::abs(Lambda * bigVec)
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
          optTarget->Params(i) = currentParameters[i] + Lambda * currentParameterDirections[i + 1];
      }
      else
      {
        for (int i = 0; i < numParams; i++)
          optparam[i] = currentParameters[i];
        for (int i = 0; i < numParams; i++)
          optdir[i] = currentParameterDirections[i + 1];
        TOL              = param_tol / bigVec;
        AbsFuncTol       = true;
        largeQuarticStep = bigChange / bigVec;
        LambdaMax        = 0.5 * Lambda;
        line_min_timer_.start();
        if (MinMethod == "quartic")
        {
          int npts(7);
          quadstep         = stepsize * Lambda;
          largeQuarticStep = bigChange / bigVec;
          Valid            = lineoptimization3(costfunc_evaluator, npts, evaluated_cost);
        }
        else
          Valid = lineoptimization2(costfunc_evaluator);
        line_min_timer_.stop();
        RealType biggestParameterChange = bigVec * std::abs(Lambda);
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
            optTarget->Params(i) = optparam[i] + Lambda * optdir[i];
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
      //APP_ABORT("QMCFixedSampleLinearOptimize::run TOO MANY FAILURES");
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
bool QMCFixedSampleLinearOptimize::put(xmlNodePtr q)
{
  std::string vmcMove("pbyp");
  std::string ReportToH5("no");
  std::string OutputMatrices("no");
  std::string FreezeParameters("no");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(vmcMove, "move");
  oAttrib.add(ReportToH5, "hdf5");

  m_param.add(OutputMatrices, "output_matrices_csv");
  m_param.add(FreezeParameters, "freeze_parameters");

  oAttrib.put(q);
  m_param.put(q);

  do_output_matrices_ = (OutputMatrices != "no");
  freeze_parameters_  = (FreezeParameters != "no");

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

  return processOptXML(q, vmcMove, ReportToH5 == "yes");
}

bool QMCFixedSampleLinearOptimize::processOptXML(xmlNodePtr opt_xml, const std::string& vmcMove, bool reportH5)
{
  m_param.put(opt_xml);

  auto iter = OptimizerNames.find(MinMethod);
  if (iter == OptimizerNames.end())
    throw std::runtime_error("Unknown MinMethod!\n");
  current_optimizer_type_ = OptimizerNames.at(MinMethod);
  if (current_optimizer_type_ == OptimizerType::DESCENT)
    throw std::runtime_error("MinMethod=descent is supported only by the batched optimization driver.\n");

  // check shift sanity
  if (shift_i_input <= 0.0)
    throw std::runtime_error("shift_i must be positive in QMCFixedSampleLinearOptimize::put");
  if (shift_s_input <= 0.0)
    throw std::runtime_error("shift_s must be positive in QMCFixedSampleLinearOptimize::put");

  // if this is the first time this function has been called, set the initial shifts
  if (current_optimizer_type_ == OptimizerType::ONESHIFTONLY)
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
  // no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers = W.getActiveWalkers();

  // Destroy old object to stop timer to correctly order timer with object lifetime scope
  vmcEngine.reset(nullptr);
  vmcEngine = std::make_unique<VMC>(project_data_, W, Psi, H, RandomNumberControl::Children, myComm, false);
  vmcEngine->setUpdateMode(vmcMove[0] == 'p');


  vmcEngine->setStatus(RootName, h5FileRoot, AppendRun);
  vmcEngine->process(qsave);

  bool success = true;
  //allways reset optTarget
  optTarget = std::make_unique<QMCCostFunction>(W, Psi, H, myComm);
  optTarget->setStream(&app_log());
  if (reportH5)
    optTarget->reportH5 = true;
  success = optTarget->put(qsave);

  return success;
}
void QMCFixedSampleLinearOptimize::one_shift_run()
{
  // ensure the cost function is set to compute derivative vectors
  optTarget->setneedGrads(true);

  // generate samples and compute weights, local energies, and derivative vectors
  start();

  // get number of optimizable parameters
  size_t numParams = optTarget->getNumParams();

  // get dimension of the linear method matrices
  size_t N = numParams + 1;

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

  if (do_output_matrices_)
  {
    output_overlap_.output(ovlMat);
    output_hamiltonian_.output(hamMat);
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
  {
    ScopedTimer local(eigenvalue_timer_);
    LinearMethod::getLowestEigenvector(prdMat, parameterDirections);
  }

  // compute the scaling constant to apply to the update
  auto lambda = LinearMethod::getNonLinearRescale(parameterDirections, ovlMat, *optTarget);

  // scale the update by the scaling constant
  for (int i = 0; i < numParams; i++)
    parameterDirections.at(i + 1) *= lambda;

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

void QMCFixedSampleLinearOptimize::start()
{
  {
    //generate samples
    ScopedTimer local(generate_samples_timer_);
    generateSamples();
    //store active number of walkers
    NumOfVMCWalkers = W.getActiveWalkers();
  }

  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  //reset the rootname
  optTarget->setRootName(RootName);
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "   Reading configurations from h5FileRoot " << h5FileRoot << std::endl;
  {
    //get configuration from the previous run
    ScopedTimer local(initialize_timer_);
    Timer t2;
    optTarget->getConfigurations(h5FileRoot);
    optTarget->setRng(vmcEngine->getRngRefs());

    // Compute wfn parameter derivatives
    NullEngineHandle handle;
    optTarget->checkConfigurations(handle);

    // check recomputed variance against VMC
    auto sigma2_vmc   = vmcEngine->getBranchEngine()->vParam[SimpleFixedNodeBranch::SBVP::SIGMA2];
    auto sigma2_check = optTarget->getVariance();
    if (optTarget->getNumSamples() > 1 && (sigma2_check > 2.0 * sigma2_vmc || sigma2_check < 0.5 * sigma2_vmc))
      throw std::runtime_error(
          "Safeguard failure: checkConfigurations variance out of [0.5, 2.0] * reference! Please report this bug.\n");
    app_log() << "  Execution time = " << std::setprecision(4) << t2.elapsed() << std::endl;
  }
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << R"(<opt stage="main" walkers=")" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}


void QMCFixedSampleLinearOptimize::finish()
{
  MyCounter++;
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;

  if (optTarget->reportH5)
    optTarget->reportParametersH5();
  optTarget->reportParameters();


  int nw_removed = W.getActiveWalkers() - NumOfVMCWalkers;
  app_log() << "   Restore the number of walkers to " << NumOfVMCWalkers << ", removing " << nw_removed << " walkers."
            << std::endl;
  if (nw_removed > 0)
    W.destroyWalkers(nw_removed);
  else
    W.createWalkers(-nw_removed);
  app_log() << "</opt>" << std::endl;
  app_log() << "</optimization-report>" << std::endl;
}

void QMCFixedSampleLinearOptimize::generateSamples()
{
  app_log() << "<optimization-report>" << std::endl;
  vmcEngine->qmc_driver_mode.set(QMC_WARMUP, 1);
  //  vmcEngine->run();
  //  vmcEngine->setValue("blocks",nBlocks);
  //  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  //  app_log() << "</vmc>" << std::endl;
  //}
  //     if (W.getActiveWalkers()>NumOfVMCWalkers)
  //     {
  //         W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
  //         app_log() << "  QMCFixedSampleLinearOptimize::generateSamples removed walkers." << std::endl;
  //         app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << std::endl;
  //     }
  vmcEngine->qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  vmcEngine->qmc_driver_mode.set(QMC_WARMUP, 0);
  //vmcEngine->setValue("recordWalkers",1);//set record
  vmcEngine->setValue("current", 0); //reset CurrentStep
  app_log() << R"(<vmc stage="main" blocks=")" << nBlocks << "\">" << std::endl;
  t1.restart();
  //     W.reset();
  branchEngine->flush(0);
  branchEngine->reset();
  vmcEngine->run();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "</vmc>" << std::endl;
  h5FileRoot = RootName;
}

} // namespace qmcplusplus
