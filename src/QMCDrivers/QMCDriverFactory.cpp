//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCDriverFactory.cpp
 * @brief Implments QMCMain operators.
 */

#include "QMCDriverFactory.h"
#include <queue>
#include "QMCDrivers/MCPopulation.h"
#include "Utilities/qmc_common.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "PsiHamNamePairReader.h"
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/RMC/RMCFactory.h"
#include "VMC/VMCDriverInput.h"
#include "VMC/VMCBatched.h"
#include "DMC/DMCDriverInput.h"
#include "DMC/DMCBatched.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/WFOpt/QMCFixedSampleLinearOptimizeBatched.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Estimators/EstimatorInputDelegates.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Message/UniformCommunicateError.h"
#include "RandomNumberControl.h"

namespace qmcplusplus
{
QMCDriverFactory::QMCDriverFactory(const ProjectData& project_data) : project_data_(project_data) {}

/** Read the xml defining the driver for this QMC section
 *
 *  Copy elision should result in just a move of the
 *  DriverAssemblyState
 *
 *  Most (all) of this should be done by calling QMCDriverInput::readXML
 *  At some point in driver refactoring this should go there and
 *  QMCDriverInput created before the giant switch
 */
QMCDriverFactory::DriverAssemblyState QMCDriverFactory::readSection(xmlNodePtr cur) const
{
  DriverAssemblyState das;
  std::string curName(castXMLCharToChar(cur->name));
  std::string update_mode("pbyp");
  std::string qmc_mode;
  std::string multi_tag("no");
  std::string warp_tag("no");
  std::string append_tag("no");
  std::string profiling_tag("no");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(qmc_mode, "method",
              {"", "vmc", "vmc_batch", "dmc", "dmc_batch", "csvmc", "rmc", "linear", "linear_batch", "wftest"});
  aAttrib.add(update_mode, "move");
  aAttrib.add(multi_tag, "multiple");
  aAttrib.add(warp_tag, "warp");
  aAttrib.add(append_tag, "append");
  aAttrib.add(profiling_tag, "profiling");
  aAttrib.add(das.traces_tag, "trace");
  aAttrib.add(das.walkerlogs_tag, "walkerlog");
  aAttrib.put(cur);
  das.append_run                 = (append_tag == "yes");
  das.enable_profiling           = (profiling_tag == "yes");
  das.what_to_do[SPACEWARP_MODE] = (warp_tag == "yes");
  das.what_to_do[MULTIPLE_MODE]  = (multi_tag == "yes");
  das.what_to_do[UPDATE_MODE]    = (update_mode == "pbyp");
  infoSummary.flush();
  infoLog.flush();

  // Really by position if you don't write <qmc name="xxx",,,>
  // you can write <vmc ...> so hacky
  if (curName != "qmc")
    qmc_mode = curName;

  const int nchars = qmc_mode.size();

  using DV = ProjectData::DriverVersion;
  switch (project_data_.getDriverVersion())
  {
  case DV::BATCH:
    if (qmc_mode.find("vmc") < nchars) // order matters here
      das.new_run_type = QMCRunType::VMC_BATCH;
    else if (qmc_mode.find("dmc") < nchars) // order matters here
      das.new_run_type = QMCRunType::DMC_BATCH;
    else if (qmc_mode.find("linear") < nchars)
      das.new_run_type = QMCRunType::LINEAR_OPTIMIZE_BATCH;
    else
      throw UniformCommunicateError("QMC mode unknown. Valid modes for batched drivers are : vmc, dmc, linear.");
    break;
  // Begin to separate driver version = batch input reading from the legacy input parsing
  case DV::LEGACY:
    if (qmc_mode.find("linear_batch") < nchars) // order matters here
      das.new_run_type = QMCRunType::LINEAR_OPTIMIZE_BATCH;
    else if (qmc_mode.find("linear") < nchars)
      das.new_run_type = QMCRunType::LINEAR_OPTIMIZE;
    else if (qmc_mode.find("rmc") < nchars)
      das.new_run_type = QMCRunType::RMC;
    else if (qmc_mode.find("csvmc") < nchars)
      das.new_run_type = QMCRunType::CSVMC;
    else if (qmc_mode.find("vmc_batch") < nchars) // order matters here
      das.new_run_type = QMCRunType::VMC_BATCH;
    else if (qmc_mode.find("vmc") < nchars)
      das.new_run_type = QMCRunType::VMC;
    else if (qmc_mode.find("dmc_batch") < nchars) // order matters here
      das.new_run_type = QMCRunType::DMC_BATCH;
    else if (qmc_mode.find("dmc") < nchars)
      das.new_run_type = QMCRunType::DMC;
    else if (qmc_mode == "wftest")
      das.new_run_type = QMCRunType::WF_TEST;
    else
      throw std::runtime_error("qmc method cannot be empty!");
  }
  return das;
}

std::unique_ptr<QMCDriverInterface> QMCDriverFactory::createQMCDriver(xmlNodePtr cur,
                                                                      DriverAssemblyState& das,
                                                                      const std::optional<EstimatorManagerInput>& emi,
                                                                      MCWalkerConfiguration& qmc_system,
                                                                      ParticleSetPool& particle_pool,
                                                                      WaveFunctionPool& wavefunction_pool,
                                                                      HamiltonianPool& hamiltonian_pool,
                                                                      Communicate* comm) const
{
  std::unique_ptr<QMCDriverInterface> new_driver;

  auto getPsi = [&wavefunction_pool](const std::string& name) -> TrialWaveFunction& {
    if (auto psi_optional = wavefunction_pool.getWaveFunction(name); psi_optional)
      return *psi_optional;
    else if (name.empty())
      throw UniformCommunicateError("Failed to find a wavefunction! Please specify the name of wavefunction using a "
                                    "qmcsystem node in the driver input.");
    else
      throw UniformCommunicateError("Failed to find the wavefunction named \"" + name + "\"!");
  };

  if (das.new_run_type == QMCRunType::CSVMC)
  { // CSVMC requires multiple pairs of Psi and Ham
    std::vector<TrialWaveFunction*> multi_psi;
    std::vector<QMCHamiltonian*> multi_ham;
    for (const auto& name_pair : PsiHamNamePairReader::readMultiplePairs(cur))
    {
      multi_psi.emplace_back(&getPsi(name_pair.first));
      multi_ham.emplace_back(hamiltonian_pool.getHamiltonian(name_pair.second));
    }
    new_driver = std::make_unique<CSVMC>(project_data_, qmc_system, std::move(multi_psi), std::move(multi_ham), comm);
    new_driver->setUpdateMode(das.what_to_do[UPDATE_MODE]);
  }
  else
  { // other drivers require at most one pair of Psi and Ham
    auto one_pair = PsiHamNamePairReader::readOnePair(cur);
    // get primaryPsi
    auto& primaryPsi = getPsi(one_pair ? one_pair->first : "");
    // get primaryH
    QMCHamiltonian* primaryH = hamiltonian_pool.getPrimary();

    auto makeEstimatorManager =
        [&](const std::optional<EstimatorManagerInput>& global_emi,
            const std::optional<EstimatorManagerInput>& driver_emi) -> UPtr<EstimatorManagerNew> {
      // This is done so that the application level input structures reflect the actual input to the code.
      // While the actual simulation objects still take singular input structures at construction.
      auto makeEstimatorManagerInput = [](auto& global_emi, auto& local_emi) -> EstimatorManagerInput {
        if (global_emi.has_value() && local_emi.has_value())
          return {global_emi.value(), local_emi.value()};
        else if (global_emi.has_value())
          return {global_emi.value()};
        else if (local_emi.has_value())
          return {local_emi.value()};
        else
          return {};
      };

      auto estimator_manager = std::make_unique<EstimatorManagerNew>(*primaryH, comm);
      estimator_manager->constructEstimators(makeEstimatorManagerInput(global_emi, driver_emi), qmc_system, primaryPsi,
                                             *primaryH, particle_pool.getPool());
      return estimator_manager;
    };

    if (das.new_run_type == QMCRunType::VMC)
    {
      new_driver = std::make_unique<VMC>(project_data_, qmc_system, primaryPsi, *primaryH,
                                         RandomNumberControl::getChildren(), comm, das.enable_profiling);
      new_driver->setUpdateMode(das.what_to_do[UPDATE_MODE]);
    }
    else if (das.new_run_type == QMCRunType::VMC_BATCH)
    {
      if (!das.what_to_do[UPDATE_MODE])
        throw UniformCommunicateError("Batched driver only supports particle-by-particle moves.");

      app_summary() << "\n========================================"
                       "\n  Reading VMC driver XML input section"
                       "\n========================================"
                    << std::endl;

      QMCDriverInput qmcdriver_input;
      VMCDriverInput vmcdriver_input;
      try
      {
        qmcdriver_input.readXML(cur);
        vmcdriver_input.readXML(cur);
      }
      catch (const std::exception& e)
      {
        throw UniformCommunicateError(e.what());
      }

      // I don't like that QMCDriverFactory is unpacking the driver input here, ideally only the driver should need to
      // depend on the content and implementation of the input.  This seems to be to be a bigger deal that passing the
      // known at this level PSPool down.
      new_driver =
          std::make_unique<VMCBatched>(project_data_, std::move(qmcdriver_input),
                                       makeEstimatorManager(emi, qmcdriver_input.get_estimator_manager_input()),
                                       std::move(vmcdriver_input), qmc_system,
                                       MCPopulation(comm->size(), comm->rank(), &qmc_system, &primaryPsi, primaryH),
                                       RandomNumberControl::getChildrenRefs(), qmc_system.getSampleStack(), comm);

      new_driver->setUpdateMode(1);
    }
    else if (das.new_run_type == QMCRunType::DMC)
    {
      DMCFactory fac(das.what_to_do[UPDATE_MODE], das.what_to_do[GPU_MODE], cur);
      new_driver = fac.create(project_data_, qmc_system, primaryPsi, *primaryH, comm, das.enable_profiling);
    }
    else if (das.new_run_type == QMCRunType::DMC_BATCH)
    {
      app_summary() << "\n========================================"
                       "\n  Reading DMC driver XML input section"
                       "\n========================================"
                    << std::endl;

      QMCDriverInput qmcdriver_input;
      DMCDriverInput dmcdriver_input;
      try
      {
        qmcdriver_input.readXML(cur);
        dmcdriver_input.readXML(cur);
      }
      catch (const std::exception& e)
      {
        throw UniformCommunicateError(e.what());
      }

      new_driver =
          std::make_unique<DMCBatched>(project_data_, std::move(qmcdriver_input),
                                       makeEstimatorManager(emi, qmcdriver_input.get_estimator_manager_input()),
                                       std::move(dmcdriver_input), qmc_system,
                                       MCPopulation(comm->size(), comm->rank(), &qmc_system, &primaryPsi, primaryH),
                                       RandomNumberControl::getChildrenRefs(), comm);
    }
    else if (das.new_run_type == QMCRunType::RMC)
    {
      RMCFactory fac(das.what_to_do[UPDATE_MODE], cur);
      new_driver = fac.create(project_data_, qmc_system, primaryPsi, *primaryH, comm);
    }
    else if (das.new_run_type == QMCRunType::LINEAR_OPTIMIZE)
    {
#ifdef MIXED_PRECISION
      APP_ABORT(
          "QMCDriverFactory::createQMCDriver : method=\"linear\" is not safe with CPU mixed precision. Please use "
          "full precision build instead.");
#endif
      QMCFixedSampleLinearOptimize* opt =
          new QMCFixedSampleLinearOptimize(project_data_, qmc_system, primaryPsi, *primaryH, comm);
      //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(qmc_system,primaryPsi,*primaryH );
      opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
      new_driver.reset(opt);
    }
    else if (das.new_run_type == QMCRunType::LINEAR_OPTIMIZE_BATCH)
    {
#ifdef MIXED_PRECISION
      APP_ABORT("QMCDriverFactory::createQMCDriver : method=\"linear_batch\" is not safe with CPU mixed precision. "
                "Please use full precision build instead.");
#endif
      app_summary() << "\n========================================"
                       "\n  Reading WFOpt driver XML input section"
                       "\n========================================"
                    << std::endl;

      QMCDriverInput qmcdriver_input;
      VMCDriverInput vmcdriver_input;
      try
      {
        qmcdriver_input.readXML(cur);
        vmcdriver_input.readXML(cur);
      }
      catch (const std::exception& e)
      {
        throw UniformCommunicateError(e.what());
      }

      if (qmcdriver_input.get_estimator_manager_input())
        app_warning() << "The batched wavefunction optimization driver ignores local estimator input.";
      if (emi)
        app_warning() << "The batched wavefunction optimization driver ignores global estimator input.";

      auto opt = std::make_unique<QMCFixedSampleLinearOptimizeBatched>(project_data_, std::move(qmcdriver_input),
                                                                       std::move(vmcdriver_input), qmc_system,
                                                                       MCPopulation(comm->size(), comm->rank(),
                                                                                    &qmc_system, &primaryPsi, primaryH),
                                                                       RandomNumberControl::getChildrenRefs(),
                                                                       qmc_system.getSampleStack(), comm);
      opt->setWaveFunctionNode(wavefunction_pool.getWaveFunctionNode("psi0"));
      new_driver = std::move(opt);
    }
    else if (das.new_run_type == QMCRunType::WF_TEST)
    {
      app_log() << "Testing wavefunctions." << std::endl;
      QMCDriverInterface* temp_ptr =
          new WaveFunctionTester(project_data_, qmc_system, primaryPsi, *primaryH, particle_pool, comm);
      new_driver.reset(temp_ptr);
    }
    else
    {
      APP_ABORT("Unhandled run type: " << static_cast<int>(das.new_run_type));
    }
  }

  infoSummary.flush();
  infoLog.flush();
  //add trace information
  bool allow_traces = das.traces_tag == "yes" ||
      (das.traces_tag == "none" && (das.new_run_type == QMCRunType::VMC || das.new_run_type == QMCRunType::DMC));
  new_driver->requestTraces(allow_traces);

  //add trace information
  bool allow_walker_logs = das.walkerlogs_tag == "yes" ||
      (das.walkerlogs_tag == "none" &&
       (das.new_run_type == QMCRunType::VMC || das.new_run_type == QMCRunType::DMC ||
        das.new_run_type == QMCRunType::VMC_BATCH || das.new_run_type == QMCRunType::DMC_BATCH));
  new_driver->requestWalkerLogs(allow_walker_logs);

  return new_driver;
}

} // namespace qmcplusplus
