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
#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/RMC/RMCFactory.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/QMCCorrelatedSamplingLinearOptimize.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include <queue>
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
///initialize the static data member
//ParticleSetPool* QMCDriverFactory::ptclPool = new ParticleSetPool;
QMCDriverFactory::QMCDriverFactory(Communicate* c) : MPIObjectBase(c), qmcSystem(0), qmcDriver(0), curRunType(QMCDriverFactory::QMCRunType::DUMMY)
{
  ////create ParticleSetPool
  ptclPool = new ParticleSetPool(myComm);
  //create WaveFunctionPool
  psiPool = new WaveFunctionPool(myComm);
  psiPool->setParticleSetPool(ptclPool);
  //create HamiltonianPool
  hamPool = new HamiltonianPool(myComm);
  hamPool->setParticleSetPool(ptclPool);
  hamPool->setWaveFunctionPool(psiPool);
}

QMCDriverFactory::~QMCDriverFactory()
{
  delete hamPool;
  delete psiPool;
  delete ptclPool;
}

void QMCDriverFactory::putCommunicator(xmlNodePtr cur)
{
  //BROKEN: myComm is ALWAYS initialized by the constructor
  if (myComm)
    return;
  ParameterSet params;
  int nparts = 1;
  params.add(nparts, "groups", "int");
  params.add(nparts, "twistAngles", "int");
  params.put(cur);
  if (nparts > 1)
  {
    app_log() << "  Communiator groups = " << nparts << std::endl;
    myComm = new Communicate(*OHMMS::Controller, nparts);
  }
}

/** Read the xml specifify the driver for this QMC section
 *
 *  Copy elision should result in just a move of the
 *  DriverAssemblyState
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
  das.append_run          = (append_tag == "yes");
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
  if (qmc_mode.find("linear") < nchars)
  {
    if (qmc_mode.find("cslinear") < nchars)
      das.new_run_type =QMCRunType::CS_LINEAR_OPTIMIZE;
    else
      das.new_run_type =QMCRunType::LINEAR_OPTIMIZE;
  }
  else if (qmc_mode.find("opt") < nchars)
  {
    das.new_run_type =QMCRunType::OPTIMIZE;
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
      das.new_run_type =QMCRunType::RMC;
    }
    else if (qmc_mode.find("vmc_batch") < nchars) // order matters here
    {
      das.new_run_type = QMCRunType::VMC_BATCH;
    }
    else if (qmc_mode.find("vmc") < nchars)
    {
      das.new_run_type =QMCRunType::VMC;
    }
    else if (qmc_mode.find("dmc") < nchars)
    {
      das.new_run_type =QMCRunType::DMC;
    }
    else if (qmc_mode == wf_test_name)
    {
      das.new_run_type =QMCRunType::WF_TEST;
    }
    else
    {
      app_log() << "Unknown qmc method: " << qmc_mode << std::endl;
    }
  }
  return das;
}

bool QMCDriverFactory::setQMCDriver(int curSeries, xmlNodePtr cur, QMCDriverFactory::DriverAssemblyState& das)
{
  unsigned long newQmcMode = das.what_to_do.to_ulong();
  //initialize to 0
  QMCDriver::BranchEngineType* branchEngine = nullptr;
  if (qmcDriver)
  {
    if (curRunType == QMCRunType::DUMMY)
    {
      APP_ABORT("QMCDriverFactory::setQMCDriver\n Other qmc sections cannot come after <qmc method=\"test\">.\n");
    }

    branchEngine = qmcDriver->getBranchEngine(); 
    delete qmcDriver;
    //set to 0 so that a new driver is created
    qmcDriver = 0;
    //if the current qmc method is different from the previous one, append_run is set to false
    
    if(das.new_run_type == curRunType && newQmcMode == curQmcMode)
      das.append_run = true;
    else
      das.append_run = false;
  }
  if (curSeries == 0)
    das.append_run = false;

  //need to create a qmcDriver
  curRunType     = das.new_run_type;
  curQmcMode     = newQmcMode;
  curQmcModeBits = das.what_to_do;
  //create a driver
  createQMCDriver(cur);
  //initialize QMCDriver::myComm
  //branchEngine has to be transferred to a new QMCDriver
  if (branchEngine)
    qmcDriver->setBranchEngine(branchEngine);
  infoSummary.flush();
  infoLog.flush();
  //add trace information
  bool allow_traces = das.traces_tag == "yes" || (das.traces_tag == "none" && (das.new_run_type ==QMCRunType::VMC || das.new_run_type ==QMCRunType::DMC));
  qmcDriver->requestTraces(allow_traces);
  return das.append_run;
}

void QMCDriverFactory::createQMCDriver(xmlNodePtr cur)
{
  ///////////////////////////////////////////////
  // get primaryPsi and primaryH
  ///////////////////////////////////////////////
  TrialWaveFunction* primaryPsi = 0;
  QMCHamiltonian* primaryH      = 0;
  std::queue<TrialWaveFunction*> targetPsi; //FIFO
  std::queue<QMCHamiltonian*> targetH;      //FIFO
  xmlNodePtr tcur = cur->children;
  while (tcur != NULL)
  {
    if (xmlStrEqual(tcur->name, (const xmlChar*)"qmcsystem"))
    {
      const xmlChar* t = xmlGetProp(tcur, (const xmlChar*)"wavefunction");
      if (t != NULL)
      {
        targetPsi.push(psiPool->getWaveFunction((const char*)t));
      }
      else
      {
        app_warning() << " qmcsystem does not have wavefunction. Assign 0" << std::endl;
        targetPsi.push(0);
      }
      t = xmlGetProp(tcur, (const xmlChar*)"hamiltonian");
      if (t != NULL)
      {
        targetH.push(hamPool->getHamiltonian((const char*)t));
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
    primaryPsi = psiPool->getPrimary();
    primaryH   = hamPool->getPrimary();
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
  //  QMCHamiltonianBase* flux=primaryH->getHamiltonian("Flux");
  //  if(flux == 0) primaryH->addOperator(new ConservedEnergy,"Flux");
  //}
  //else
  //{
  //  primaryH->remove("Flux");
  //}
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  if (curRunType ==QMCRunType::VMC || curRunType ==QMCRunType::CSVMC)
  {
    //VMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
    VMCFactory fac(curQmcModeBits.to_ulong(), cur);
    qmcDriver = fac.create(*qmcSystem, *primaryPsi, *primaryH, *ptclPool, *hamPool, *psiPool, myComm);
    //TESTING CLONE
    //TrialWaveFunction* psiclone=primaryPsi->makeClone(*qmcSystem);
    //qmcDriver = fac.create(*qmcSystem,*psiclone,*primaryH,*ptclPool,*hamPool);
  }
  else if (curRunType == QMCRunType::VMC_BATCH)
  {
    APP_ABORT("VMCBatch driver not yet supported");
    // DFCreator<VMC_BATCH> dfc(curQmcModeBits.to_ulong(), cur);
    // qmcDriver = dfc(*qmcSystem, *primaryPsi, *primaryH, *ptclPool, *hamPool, *psiPool, myComm);

  }
  else if (curRunType ==QMCRunType::DMC)
  {
    DMCFactory fac(curQmcModeBits[UPDATE_MODE], curQmcModeBits[GPU_MODE], cur);
    qmcDriver = fac.create(*qmcSystem, *primaryPsi, *primaryH, *hamPool, *psiPool, myComm);
  }
  else if (curRunType ==QMCRunType::RMC)
  {
    RMCFactory fac(curQmcModeBits[UPDATE_MODE], cur);
    qmcDriver = fac.create(*qmcSystem, *primaryPsi, *primaryH, *ptclPool, *hamPool, *psiPool, myComm);
  }
  else if (curRunType ==QMCRunType::OPTIMIZE)
  {
    QMCOptimize* opt = new QMCOptimize(*qmcSystem, *primaryPsi, *primaryH, *hamPool, *psiPool, myComm);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(*qmcSystem,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver = opt;
  }
  else if (curRunType ==QMCRunType::LINEAR_OPTIMIZE)
  {
#ifdef MIXED_PRECISION
    APP_ABORT("QMCDriverFactory::createQMCDriver : method=\"linear\" is not safe with CPU mixed precision. Please use "
              "full precision build instead.");
#endif
    QMCFixedSampleLinearOptimize* opt =
        new QMCFixedSampleLinearOptimize(*qmcSystem, *primaryPsi, *primaryH, *hamPool, *psiPool, myComm);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(*qmcSystem,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver = opt;
  }
  else if (curRunType ==QMCRunType::CS_LINEAR_OPTIMIZE)
  {
#if defined(QMC_CUDA)
    app_log() << "cslinear is not supported. Switch to linear method. " << std::endl;
    QMCFixedSampleLinearOptimize* opt =
        new QMCFixedSampleLinearOptimize(*qmcSystem, *primaryPsi, *primaryH, *hamPool, *psiPool, myComm);
#else
    QMCCorrelatedSamplingLinearOptimize* opt =
        new QMCCorrelatedSamplingLinearOptimize(*qmcSystem, *primaryPsi, *primaryH, *hamPool, *psiPool, myComm);
#endif
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver = opt;
  }
  else if (curRunType ==QMCRunType::WF_TEST)
  {
    app_log() << "Testing wavefunctions." << std::endl;
    qmcDriver = new WaveFunctionTester(*qmcSystem, *primaryPsi, *primaryH, *ptclPool, *psiPool, myComm);
  }
  else
  {
    APP_ABORT("Unhandled run type: " << static_cast<int>(curRunType));
  }
  if (curQmcModeBits[MULTIPLE_MODE])
  {
    while (targetH.size())
    {
      qmcDriver->add_H_and_Psi(targetH.front(), targetPsi.front());
      targetH.pop();
      targetPsi.pop();
    }
  }
}
} // namespace qmcplusplus
