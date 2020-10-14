//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCDriverFactory.cpp
 * @brief Implments QMCMain operators.
 */
#include <queue>


#include "QMCDrivers/MCPopulation.h"
#include "Utilities/qmc_common.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/VMC/VMCFactoryNew.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/DMC/DMCFactoryNew.h"
#include "QMCDrivers/RMC/RMCFactory.h"
#include "QMCDrivers/WFOpt/QMCOptimize.h"
#include "QMCDrivers/WFOpt/QMCOptimizeBatched.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"
#include "QMCDrivers/WFOpt/QMCCorrelatedSamplingLinearOptimize.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCDrivers/WFOpt/QMCWFOptFactoryNew.h"

namespace qmcplusplus
{
/** Read the xml specifify the driver for this QMC section
 *
 *  Copy elision should result in just a move of the
 *  DriverAssemblyState
 *
 *  Most (all) of this should be done by calling QMCDriverInput::readXML
 *  At some point in driver refactoring this should go there and
 *  QMCDriverInput created before the giant switch
 */
QMCDriverFactory::DriverAssemblyState QMCDriverFactory::readSection(int curSeries, xmlNodePtr cur)
{
  DriverAssemblyState das;
  std::string curName((const char*)cur->name);
  std::string update_mode("pbyp");
  std::string qmc_mode("invalid");
  std::string multi_tag("no");
  std::string warp_tag("no");
  std::string append_tag("no");
#if defined(QMC_CUDA)
  std::string gpu_tag("yes");
#else
  std::string gpu_tag("no");
#endif
  OhmmsAttributeSet aAttrib;
  aAttrib.add(qmc_mode, "method");
  aAttrib.add(update_mode, "move");
  aAttrib.add(multi_tag, "multiple");
  aAttrib.add(warp_tag, "warp");
  aAttrib.add(append_tag, "append");
  aAttrib.add(gpu_tag, "gpu");
  aAttrib.add(das.traces_tag, "trace");
  aAttrib.put(cur);
  das.append_run                 = (append_tag == "yes");
  das.what_to_do[SPACEWARP_MODE] = (warp_tag == "yes");
  das.what_to_do[MULTIPLE_MODE]  = (multi_tag == "yes");
  das.what_to_do[UPDATE_MODE]    = (update_mode == "pbyp");
#if defined(QMC_CUDA)
  das.what_to_do[GPU_MODE] = (gpu_tag == "yes");
#endif
  infoSummary.flush();
  infoLog.flush();

  std::string wf_test_name("wftest");

  if (curName != "qmc")
    qmc_mode = curName;
  int nchars = qmc_mode.size();
  if (qmc_mode.find("linear_batch") < nchars)
  {
    das.new_run_type = QMCRunType::LINEAR_OPTIMIZE_BATCH;
  }
  else if (qmc_mode.find("linear") < nchars)
  {
    if (qmc_mode.find("cslinear") < nchars)
      das.new_run_type = QMCRunType::CS_LINEAR_OPTIMIZE;
    else
      das.new_run_type = QMCRunType::LINEAR_OPTIMIZE;
  }
  else if (qmc_mode.find("opt_batch") < nchars)
  {
    das.new_run_type = QMCRunType::OPTIMIZE_BATCH;
  }
  else if (qmc_mode.find("opt") < nchars)
  {
    das.new_run_type = QMCRunType::OPTIMIZE;
  }
  else
  {
    if (qmc_mode.find("ptcl") < nchars)
      das.what_to_do[UPDATE_MODE] = 1;
    if (qmc_mode.find("mul") < nchars)
      das.what_to_do[MULTIPLE_MODE] = 1;
    if (qmc_mode.find("warp") < nchars)
      das.what_to_do[SPACEWARP_MODE] = 1;
    //       if (qmc_mode.find("rmcPbyP")<nchars)
    //       {
    //         das.new_run_type=RMC_PBYP_RUN;
    //       }
    //       else
    if (qmc_mode.find("rmc") < nchars)
    {
      das.new_run_type = QMCRunType::RMC;
    }
    else if (qmc_mode.find("vmc_batch") < nchars) // order matters here
    {
      das.new_run_type = QMCRunType::VMC_BATCH;
    }
    else if (qmc_mode.find("vmc") < nchars)
    {
      das.new_run_type = QMCRunType::VMC;
    }
    else if (qmc_mode.find("dmc_batch") < nchars) // order matters here
    {
      das.new_run_type = QMCRunType::DMC_BATCH;
    }
    else if (qmc_mode.find("dmc") < nchars)
    {
      das.new_run_type = QMCRunType::DMC;
    }
    else if (qmc_mode == wf_test_name)
    {
      das.new_run_type = QMCRunType::WF_TEST;
    }
    else
    {
      app_log() << "Unknown qmc method: " << qmc_mode << std::endl;
    }
  }
  return das;
}

std::unique_ptr<QMCDriverInterface> QMCDriverFactory::newQMCDriver(std::unique_ptr<QMCDriverInterface> last_driver,
                                                                   int curSeries,
                                                                   xmlNodePtr cur,
                                                                   QMCDriverFactory::DriverAssemblyState& das,
                                                                   MCWalkerConfiguration& qmc_system,
                                                                   ParticleSetPool& particle_pool,
                                                                   WaveFunctionPool& wavefunction_pool,
                                                                   HamiltonianPool& hamiltonian_pool,
                                                                   MCPopulation& population,
                                                                   Communicate* comm)
{
  //create a driver
  std::unique_ptr<QMCDriverInterface> new_driver =
      createQMCDriver(cur, das, qmc_system, particle_pool, wavefunction_pool, hamiltonian_pool, population, comm);
  //initialize QMCDriver::myComm
  // branchEngine has to be transferred to a new QMCDriver
  // but we also have to deal with it splitting into new and legacy classes

  if (last_driver && last_driver->getRunType() == QMCRunType::DUMMY)
  {
    throw std::runtime_error("QMCDriverFactory::newQMCDriver\n Other qmc sections cannot come after <qmc method=\"test\">.\n");
  }

  switch (das.new_run_type)
  {
  case QMCRunType::DUMMY:
  case QMCRunType::VMC:
  case QMCRunType::CSVMC:
  case QMCRunType::DMC:
  case QMCRunType::RMC:
  case QMCRunType::OPTIMIZE:
  case QMCRunType::VMC_OPT:
  case QMCRunType::LINEAR_OPTIMIZE:
  case QMCRunType::CS_LINEAR_OPTIMIZE:
  case QMCRunType::WF_TEST:
  {
    QMCDriver::BranchEngineType* branchEngine = nullptr;
    if (last_driver)
    {
      branchEngine = last_driver->getBranchEngine();
      branchEngine->resetRun(cur);
    }
    if (branchEngine)
      new_driver->setBranchEngine(branchEngine);
  }
  break;
  case QMCRunType::VMC_BATCH:
  case QMCRunType::DMC_BATCH:
  case QMCRunType::OPTIMIZE_BATCH:
  case QMCRunType::LINEAR_OPTIMIZE_BATCH:
  {
    SFNBranch* branchEngine = nullptr;
    if (last_driver)
    {
      //Checking the QMCRunType gives us the type information
      QMCDriverNew& old_unified_driver = static_cast<QMCDriverNew&>(*last_driver);

      branchEngine = old_unified_driver.getNewBranchEngine();
      // Someone helpful added a resetRun call here, a call from legacy CUDA that if it did help I need to know why.
      // It did seem to fix an issue with legacy which is interesting....
    }
    if (branchEngine)
    {
      QMCDriverNew& new_unified_driver = static_cast<QMCDriverNew&>(*new_driver);
      new_unified_driver.setNewBranchEngine(branchEngine);
    }
  }
  break;
  }
  //initialize to 0

  infoSummary.flush();
  infoLog.flush();
  //add trace information
  bool allow_traces = das.traces_tag == "yes" ||
      (das.traces_tag == "none" && (das.new_run_type == QMCRunType::VMC || das.new_run_type == QMCRunType::DMC));
  new_driver->requestTraces(allow_traces);
  return new_driver;
}

std::unique_ptr<QMCDriverInterface> QMCDriverFactory::createQMCDriver(xmlNodePtr cur,
                                                                      DriverAssemblyState& das,
                                                                      MCWalkerConfiguration& qmc_system,
                                                                      ParticleSetPool& particle_pool,
                                                                      WaveFunctionPool& wavefunction_pool,
                                                                      HamiltonianPool& hamiltonian_pool,
                                                                      MCPopulation& population,
                                                                      Communicate* comm)
{
  ///////////////////////////////////////////////
  // get primaryPsi and primaryH
  ///////////////////////////////////////////////
  TrialWaveFunction* primaryPsi = 0;
  QMCHamiltonian* primaryH      = 0;
  std::queue<TrialWaveFunction*> targetPsi; //FIFO
  std::queue<QMCHamiltonian*> targetH;      //FIFO
  xmlNodePtr tcur = cur->children;
  std::unique_ptr<QMCDriverInterface> new_driver;
  while (tcur != NULL)
  {
    if (xmlStrEqual(tcur->name, (const xmlChar*)"qmcsystem"))
    {
      const XMLAttrString wf_name(tcur, "wavefunction");
      if (!wf_name.empty())
      {
        targetPsi.push(wavefunction_pool.getWaveFunction(wf_name));
      }
      else
      {
        app_warning() << " qmcsystem does not have wavefunction. Assign 0" << std::endl;
        targetPsi.push(0);
      }
      const XMLAttrString ham_name(tcur, "hamiltonian");
      if (!ham_name.empty())
      {
        targetH.push(hamiltonian_pool.getHamiltonian(ham_name));
      }
      else
      {
        app_warning() << " qmcsystem does not have hamiltonian. Assign 0" << std::endl;
        targetH.push(0);
      }
    }
    tcur = tcur->next;
  }
  //mark the first targetPsi and targetH as the primaries
  if (targetH.empty())
  {
    primaryPsi = wavefunction_pool.getPrimary();
    primaryH   = hamiltonian_pool.getPrimary();
  }
  else
  {
    primaryPsi = targetPsi.front();
    // targetPsi.pop();
    primaryH = targetH.front();
    //  targetH.pop();
  }
  //set primaryH->Primary
  primaryH->setPrimary(true);

  ////flux is evaluated only with single-configuration VMC
  //if(curRunType ==QMCRunType::VMC && !curQmcModeBits[MULTIPLE_MODE])
  //{
  //  OperatorBase* flux=primaryH->getHamiltonian("Flux");
  //  if(flux == 0) primaryH->addOperator(new ConservedEnergy,"Flux");
  //}
  //else
  //{
  //  primaryH->remove("Flux");
  //}
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  if (das.new_run_type == QMCRunType::VMC || das.new_run_type == QMCRunType::CSVMC)
  {
    //VMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
    VMCFactory fac(das.what_to_do.to_ulong(), cur);
    new_driver.reset(
        fac.create(qmc_system, *primaryPsi, *primaryH, comm));
    //TESTING CLONE
    //TrialWaveFunction* psiclone=primaryPsi->makeClone(qmc_system);
    //qmcDriver = fac.create(qmc_system,*psiclone,*primaryH,particle_pool,hamiltonian_pool);
  }
  else if (das.new_run_type == QMCRunType::VMC_BATCH)
  {
    VMCFactoryNew fac(cur, das.what_to_do[UPDATE_MODE], qmc_common.qmc_counter);
    new_driver.reset(
        fac.create(population, *primaryPsi, *primaryH, qmc_system.getSampleStack(), comm));
  }
  else if (das.new_run_type == QMCRunType::DMC)
  {
    DMCFactory fac(das.what_to_do[UPDATE_MODE], das.what_to_do[GPU_MODE], cur);
    new_driver.reset(fac.create(qmc_system, *primaryPsi, *primaryH, comm));
  }
  else if (das.new_run_type == QMCRunType::DMC_BATCH)
  {
    DMCFactoryNew fac(cur, das.what_to_do[UPDATE_MODE], qmc_common.qmc_counter);
    new_driver.reset(fac.create(population, *primaryPsi, *primaryH, comm));
  }
  else if (das.new_run_type == QMCRunType::RMC)
  {
    RMCFactory fac(das.what_to_do[UPDATE_MODE], cur);
    new_driver.reset(
        fac.create(qmc_system, *primaryPsi, *primaryH, comm));
  }
  else if (das.new_run_type == QMCRunType::OPTIMIZE)
  {
    QMCOptimize* opt = new QMCOptimize(qmc_system, *primaryPsi, *primaryH, comm);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(qmc_system,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
    new_driver.reset(opt);
  }
  else if (das.new_run_type == QMCRunType::OPTIMIZE_BATCH)
  {
    QMCOptimizeBatched* opt =
        QMCWFOptFactoryNew(cur, qmc_common.qmc_counter, qmc_system, *primaryPsi, *primaryH, population, qmc_system.getSampleStack(), comm);
    opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
    new_driver.reset(opt);
  }
  else if (das.new_run_type == QMCRunType::LINEAR_OPTIMIZE)
  {
#ifdef MIXED_PRECISION
    APP_ABORT("QMCDriverFactory::createQMCDriver : method=\"linear\" is not safe with CPU mixed precision. Please use "
              "full precision build instead.");
#endif
    QMCFixedSampleLinearOptimize* opt =
        new QMCFixedSampleLinearOptimize(qmc_system, *primaryPsi, *primaryH, comm);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(qmc_system,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
    new_driver.reset(opt);
  }
  else if (das.new_run_type == QMCRunType::LINEAR_OPTIMIZE_BATCH)
  {
#ifdef MIXED_PRECISION
    APP_ABORT(
        "QMCDriverFactory::createQMCDriver : method=\"linear_batch\" is not safe with CPU mixed precision. Please use "
        "full precision build instead.");
#endif
    QMCFixedSampleLinearOptimizeBatched* opt =
        QMCWFOptLinearFactoryNew(cur, qmc_common.qmc_counter, qmc_system, *primaryPsi, *primaryH, population, qmc_system.getSampleStack(), comm);
    opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
    new_driver.reset(opt);
  }
  else if (das.new_run_type == QMCRunType::CS_LINEAR_OPTIMIZE)
  {
#if defined(QMC_CUDA)
    app_log() << "cslinear is not supported. Switch to linear method. " << std::endl;
    QMCFixedSampleLinearOptimize* opt =
        new QMCFixedSampleLinearOptimize(qmc_system, *primaryPsi, *primaryH, comm);
#else
    QMCCorrelatedSamplingLinearOptimize* opt =
        new QMCCorrelatedSamplingLinearOptimize(qmc_system, *primaryPsi, *primaryH, comm);
#endif
    opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
    new_driver.reset(opt);
  }
  else if (das.new_run_type == QMCRunType::WF_TEST)
  {
    app_log() << "Testing wavefunctions." << std::endl;
    QMCDriverInterface* temp_ptr =
        new WaveFunctionTester(qmc_system, *primaryPsi, *primaryH, particle_pool, comm);
    new_driver.reset(temp_ptr);
  }
  else
  {
    APP_ABORT("Unhandled run type: " << static_cast<int>(das.new_run_type));
  }
  if (das.what_to_do[MULTIPLE_MODE])
  {
    while (targetH.size())
    {
      new_driver->add_H_and_Psi(targetH.front(), targetPsi.front());
      targetH.pop();
      targetPsi.pop();
    }
  }
  return new_driver;
}
} // namespace qmcplusplus
