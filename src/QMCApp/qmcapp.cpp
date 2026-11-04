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

#include <stdexcept>
#include <memory>
#include <optional>
#include <iostream>
#include "Configuration.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/ProgressReportEngine.h"
#include "Platforms/Host/OutputManager.h"
#include "OhmmsData/FileUtility.h"
#include "Host/sysutil.h"
#include "ProjectData.h"
#include "DeviceManager.h"
#include "QMCApp/QMCMain.h"
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
    Communicate node_comm{OHMMS::Controller->createNodeComm()};
    // assign accelerators within a physical node
    DeviceManager::initializeGlobalDeviceManager(node_comm.rank(), node_comm.size());

    //qmc_common  and MPI is initialized
    qmcplusplus::qmc_common.initialize(argc, argv);
    std::vector<std::string> fgroup_names_cmd, fgroup_names_txt;
    std::optional<int> vnode_ranks;

    int i = 1;
    while (i < argc)
    {
      std::string c(argv[i]);
      if (c[0] == '-')
      {
        if (c == "-debug")
          ReportEngine::enableOutput();

        // Default setting is 'timer_level_coarse'
        if (c.find("-enable-timers") < c.size())
        {
#ifndef ENABLE_TIMERS
          std::cerr
              << "The '-enable-timers' command line option will have no effect. This executable was built without "
                 "ENABLE_TIMER set."
              << std::endl;
#endif
          int pos = c.find("=");
          if (pos != std::string::npos)
          {
            std::string timer_level = c.substr(pos + 1);
            getGlobalTimerManager().set_timer_threshold(timer_level);
          }
        }
        if (c.find("-verbosity") < c.size())
        {
          int pos = c.find("=");
          if (pos != std::string::npos)
          {
            std::string verbose_level = c.substr(pos + 1);
            if (verbose_level == "low")
            {
              outputManager.setVerbosity(Verbosity::LOW);
            }
            else if (verbose_level == "high")
            {
              outputManager.setVerbosity(Verbosity::HIGH);
            }
            else if (verbose_level == "debug")
            {
              outputManager.setVerbosity(Verbosity::DEBUG);
            }
            else
            {
              throw UniformCommunicateError("Unknown verbosity level: " + verbose_level);
            }
          }
        }
        if (c.find("vnode-ranks") < c.size())
        {
          int pos = c.find("=");
          try
          {
            if (pos != std::string::npos)
              vnode_ranks = std::stoi(c.substr(pos + 1));
            else if (i + 1 < argc)
              vnode_ranks = std::stoi(argv[++i]);
            else
              throw UniformCommunicateError("Missing argument for --vnode-ranks");
          }
          catch (const std::invalid_argument&)
          {
            throw UniformCommunicateError("Invalid argument for --vnode-ranks");
          }
          catch (const std::out_of_range&)
          {
            throw UniformCommunicateError("Argument out of range for --vnode-ranks");
          }
          if (*vnode_ranks <= 0)
            throw UniformCommunicateError("--vnode-ranks must be a positive integer.");
        }
      }
      else
      {
        if (c.find(".xml") == c.size() - 4)
          fgroup_names_cmd.push_back(argv[i]);
        else
        {
          std::ifstream fin(argv[i], std::ifstream::in);
          bool valid = !fin.fail();
          while (valid)
          {
            std::vector<std::string> words;
            getwords(words, fin);
            if (words.size())
            {
              if (words[0].find(".xml") == words[0].size() - 4)
                fgroup_names_txt.push_back(words[0]);
            }
            else
              valid = false;
          }
        }
      }
      ++i;
    }

    std::vector<std::string> inputs(fgroup_names_cmd.size() + fgroup_names_txt.size());
    copy(fgroup_names_txt.begin(), fgroup_names_txt.end(), inputs.begin());
    i = fgroup_names_txt.size();
    for (int k = 0; k < fgroup_names_cmd.size(); ++k)
      inputs[i++] = fgroup_names_cmd[k];
    if (inputs.empty())
    {
      std::string msg = "No valid input file is given.\n"
                        "Usage: qmcpack [options] <input-files.xml>\n"
                        "Ensemble runs may be initialized by specifying either multiple .xml files or text files "
                        "containing lists of .xml input files.";
      throw UniformCommunicateError(msg);
    }

    // set virtual node rank count
    if (const int max_vnode_ranks = OHMMS::Controller->getUniformNodeRankCount(); vnode_ranks.has_value())
    {
      if (max_vnode_ranks % *vnode_ranks != 0)
      {
        std::ostringstream msg;
        msg << "main(). --vnode-ranks (" << *vnode_ranks
            << ") can only be a factor of the uniform node rank count (" << max_vnode_ranks << ").";
        throw UniformCommunicateError(msg.str());
      }
      qmc_common.set_virtual_node_rank_count(*vnode_ranks);
    }
    else
      qmc_common.set_virtual_node_rank_count(max_vnode_ranks);

    if (inputs.size() > 1)
    {
      const int total_blocks = OHMMS::Controller->size() / qmc_common.get_virtual_node_rank_count();
      if (total_blocks < inputs.size())
      {
        std::ostringstream msg;
        msg << "main(). Current " << OHMMS::Controller->size() << " MPI ranks (" << total_blocks << " blocks of size "
            << qmc_common.get_virtual_node_rank_count() << ") cannot accommodate all the " << inputs.size()
            << " individual calculations in the ensemble. "
            << "Increase the number of MPI ranks or reduce the number of calculations.";
        throw UniformCommunicateError(msg.str());
      }
    }

    std::unique_ptr<Communicate> comm_job;
    if (inputs.size() > 1)
    {
      comm_job =
          std::make_unique<Communicate>(OHMMS::Controller->split(inputs.size(), qmc_common.get_virtual_node_rank_count()));
      qmc_common.setupFilenameInfixForEnsembleJobs(comm_job->getGroupID(), comm_job->getNumGroups());
    }

    Communicate& qmcComm          = inputs.size() > 1 ? *comm_job : *OHMMS::Controller;
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
