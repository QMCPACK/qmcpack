//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"

#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"

//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionCUDA.h"
#endif
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{
QMCLinearOptimize::QMCLinearOptimize(MCWalkerConfiguration& w,
                                     TrialWaveFunction& psi,
                                     QMCHamiltonian& h,
                                     Communicate* comm,
                                     const std::string& QMC_driver_type)
    : QMCDriver(w, psi, h, comm, QMC_driver_type), wfNode(NULL), param_tol(1e-4)
{
  IsQMCDriver = false;
  //     //set the optimization flag
  qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  m_param.add(param_tol, "alloweddifference");
  //Set parameters for line minimization:
}

void QMCLinearOptimize::start()
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
    optTarget->checkConfigurations();
    // check recomputed variance against VMC
    auto sigma2_vmc   = vmcEngine->getBranchEngine()->vParam[SimpleFixedNodeBranch::SBVP::SIGMA2];
    auto sigma2_check = optTarget->getVariance();
    if (sigma2_check > 2.0 * sigma2_vmc || sigma2_check < 0.5 * sigma2_vmc)
      throw std::runtime_error(
          "Safeguard failure: checkConfigurations variance out of [0.5, 2.0] * reference! Please report this bug.\n");
    app_log() << "  Execution time = " << std::setprecision(4) << t2.elapsed() << std::endl;
  }
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << "<opt stage=\"main\" walkers=\"" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}

#ifdef HAVE_LMY_ENGINE
void QMCLinearOptimize::engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                                     DescentEngine& descentEngineObj,
                                     std::string MinMethod)
{
  app_log() << "entering engine_start function" << std::endl;

  {
    //generate samples
    ScopedTimer local(generate_samples_timer_);
    generateSamples();
    //store active number of walkers
    NumOfVMCWalkers = W.getActiveWalkers();
  }

  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;

  // reset the root name
  optTarget->setRootName(RootName);
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "     Reading configurations from h5FileRoot " << h5FileRoot << std::endl;
  {
    // get configuration from the previous run
    ScopedTimer local(initialize_timer_);
    Timer t2;
    optTarget->getConfigurations(h5FileRoot);
    optTarget->setRng(vmcEngine->getRngRefs());
    optTarget->engine_checkConfigurations(EngineObj, descentEngineObj,
                                          MinMethod); // computes derivative ratios and pass into engine
    app_log() << "  Execution time = " << std::setprecision(4) << t2.elapsed() << std::endl;
  }
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << "<opt stage=\"main\" walkers=\"" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}
#endif

void QMCLinearOptimize::finish()
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

void QMCLinearOptimize::generateSamples()
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
  //         app_log() << "  QMCLinearOptimize::generateSamples removed walkers." << std::endl;
  //         app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << std::endl;
  //     }
  vmcEngine->qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  vmcEngine->qmc_driver_mode.set(QMC_WARMUP, 0);
  //vmcEngine->setValue("recordWalkers",1);//set record
  vmcEngine->setValue("current", 0); //reset CurrentStep
  app_log() << "<vmc stage=\"main\" blocks=\"" << nBlocks << "\">" << std::endl;
  t1.restart();
  //     W.reset();
  branchEngine->flush(0);
  branchEngine->reset();
  vmcEngine->run();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "</vmc>" << std::endl;
  h5FileRoot = RootName;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool QMCLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("no");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU, "gpu");
  oAttrib.add(vmcMove, "move");
  oAttrib.put(q);
  xmlNodePtr optNode = q;
  xmlNodePtr cur     = optNode->children;
  int pid            = OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    cur = cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers = W.getActiveWalkers();
  bool success    = true;
  //allways reset optTarget
#if defined(QMC_CUDA)
  if (useGPU == "yes")
    optTarget = std::make_unique<QMCCostFunctionCUDA>(W, Psi, H, myComm);
  else
#endif
    optTarget = std::make_unique<QMCCostFunction>(W, Psi, H, myComm);
  optTarget->setStream(&app_log());
  success = optTarget->put(q);

  //create VMC engine
  if (vmcEngine == 0)
  {
#if defined(QMC_CUDA)
    if (useGPU == "yes")
      vmcEngine = std::make_unique<VMCcuda>(W, Psi, H, myComm, false);
    else
#endif
      vmcEngine = std::make_unique<VMC>(W, Psi, H, myComm, false);
    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
  }

  vmcEngine->setStatus(RootName, h5FileRoot, AppendRun);
  vmcEngine->process(optNode);
  return success;
}

} // namespace qmcplusplus
