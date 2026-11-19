//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include "Configuration.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "Platforms/Host/OutputManager.h"
#include "DeviceManager.h"
#include "QMCApp/CommandLineOptions.h"
#include "QMCApp/setupQMCState.h"
#include "QMCApp/QMCMain.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/TimerManager.h"
#include "Utilities/qmc_common.h"
#include "Concurrency/OpenMP.h"

void output_hardware_info(Communicate& comm, Libxml2Document& doc, xmlNodePtr root);

/** @file qmcapp.cpp
 *@brief a main function for QMC simulation.
 *
 * @ingroup qmcapp
 * @brief main function for qmcapp executable.
 *
 *Actual works are done by QMCAppBase and its derived classe.
 *For other simulations, one can derive a class from QMCApps, similarly to MolecuApps.
 */
int main(int argc, char** argv)
{
  using namespace qmcplusplus;
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv, boost::mpi3::thread_level::funneled);
  OHMMS::Controller = std::make_unique<Communicate>(env.world());
  if (const auto [min_ranks, max_ranks] = OHMMS::Controller->getMinMaxNodeRankCount();
      min_ranks != max_ranks && OHMMS::Controller->rank() == 0)
    std::cerr << "Warning: Unequal number of MPI ranks per node (min: " << min_ranks << ", max: " << max_ranks << ")."
              << std::endl;
#endif
  try
  {
    CommandLineOptions options;
    options.parse(argc, argv);
    if (options.should_normal_stop())
    {
      if (OHMMS::Controller->rank() == 0)
        options.print_version_and_usage(std::cout);
      OHMMS::Controller.reset();
      return 0;
    }

    if (options.get_debug())
      ReportEngine::enableOutput();

    // Default setting is 'timer_level_coarse'
    if (options.get_enable_timers())
    {
#ifndef ENABLE_TIMERS
      if (OHMMS::Controller->rank() == 0)
        std::cerr
            << "The '-enable-timers' command line option will have no effect. This executable was built without "
               "ENABLE_TIMER set."
            << std::endl;
#endif
      if (!options.get_timer_level().empty())
        getGlobalTimerManager().set_timer_threshold(options.get_timer_level());
    }

    if (options.get_verbosity().has_value())
      outputManager.setVerbosity(*options.get_verbosity());

    setupQMCState(options, *OHMMS::Controller, qmc_common);

    Communicate node_comm{OHMMS::Controller->createNodeComm()};
    // assign accelerators within a physical node
    DeviceManager::initializeGlobalDeviceManager(node_comm.rank(), node_comm.size());

    const auto& inputs = options.get_inputs();

    std::unique_ptr<Communicate> comm_job;
    if (inputs.size() > 1 || options.get_vnode_ranks().has_value())
    {
      const int vnode_ranks = options.get_vnode_ranks().value_or(OHMMS::Controller->getVirtualNodeRankCount());
      comm_job              = std::make_unique<Communicate>(OHMMS::Controller->split(inputs.size(), vnode_ranks));
      if (inputs.size() > 1)
        qmc_common.setupFilenameInfixForEnsembleJobs(comm_job->getGroupID(), comm_job->getNumGroups());
    }

    Communicate& qmcComm          = comm_job ? *comm_job : *OHMMS::Controller;
    const std::string& input_name = inputs.size() > 1 ? inputs[qmcComm.getGroupID()] : inputs[0];

    if (qmcComm.rank() != 0)
      outputManager.shutOff(); // might need to redirect debug stream to a file per rank if debugging is enabled
    else if(inputs.size() > 1)
    {
      const std::string fn =
          input_name.substr(0, input_name.size() - 4) + qmc_common.get_filename_infix_for_ensemble_jobs() + ".qmc";
      infoSummary.redirectToFile(fn);
      infoLog.redirectToSameStream(infoSummary);
      infoError.redirectToSameStream(infoSummary);
    }

    app_log() << "  Input file(s): ";
    for (int k = 0; k < inputs.size(); ++k)
      app_log() << inputs[k] << " ";
    app_log() << std::endl;

    auto qmc = std::make_unique<QMCMain>(qmcComm);

    try
    {
      if (!qmc->parse(input_name))
        throw UniformCommunicateError("main(). Input invalid.");

      if (!qmc->execute())
        throw UniformCommunicateError("main(). QMC execution failed.");
    }
    catch (const UniformCommunicateError& ue)
    {
      std::ostringstream msg;
      if (inputs.size() > 1)
        msg << "Ensemble job " << qmcComm.getGroupID() << ", input file " << input_name << ": ";
      else
        msg << "Input file " << input_name << ": ";
      msg << ue.what();
      qmcComm.barrier_and_abort(msg.str());
    }

    Libxml2Document timingDoc;
    timingDoc.newDoc("resources");
    output_hardware_info(qmcComm, timingDoc, timingDoc.getRoot());
    getGlobalTimerManager().output_timing(qmcComm, timingDoc, timingDoc.getRoot());
    qmc->getParticlePool().output_particleset_info(timingDoc, timingDoc.getRoot());
    if (OHMMS::Controller->rank() == 0)
    {
      timingDoc.dump(qmc->getTitle() + ".info.xml");
    }
    getGlobalTimerManager().print(qmcComm);

    qmc.reset();
  }
  catch (const UniformCommunicateError& ue)
  {
    OHMMS::Controller->barrier_and_abort("Global parallel scope: " + std::string(ue.what()));
  }
  catch (const std::exception& e)
  {
    APP_ABORT("Global rank " + std::to_string(OHMMS::Controller->rank()) + " Unhandled Exception: " + e.what());
  }
  catch (...)
  {
    APP_ABORT("Global rank " + std::to_string(OHMMS::Controller->rank()) +
              " Unhandled Exception (not derived from std::exception)");
  }

  OHMMS::Controller->barrier(); // Barrier needed to ensure all groups in ensemble runs complete successfully

  if (OHMMS::Controller->rank() == 0)
    std::cout << std::endl << "QMCPACK execution completed successfully" << std::endl;

  // OHMMS::Controller holds duplicated (linked) mpi3::environment and must be freed first.
  OHMMS::Controller.reset();
  return 0;
}

void output_hardware_info(Communicate& comm, Libxml2Document& doc, xmlNodePtr root)
{
  xmlNodePtr hardware = doc.addChild(root, "hardware");

  bool using_mpi = false;
#ifdef HAVE_MPI
  using_mpi = true;
  doc.addChild(hardware, "mpi_size", comm.size());
#endif
  doc.addChild(hardware, "mpi", using_mpi);

  bool using_openmp = false;
#ifdef _OPENMP
  using_openmp = true;
  doc.addChild(hardware, "openmp_threads", omp_get_max_threads());
#endif
  doc.addChild(hardware, "openmp", using_openmp);
}
