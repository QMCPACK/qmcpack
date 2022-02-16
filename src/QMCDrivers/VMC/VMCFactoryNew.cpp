//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "VMCFactoryNew.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "Concurrency/Info.hpp"

namespace qmcplusplus
{
QMCDriverInterface* VMCFactoryNew::create(const ProjectData& project_data,
                                          MCPopulation&& pop,
                                          SampleStack& samples,
                                          Communicate* comm)
{
#if defined(QMC_CUDA)
  comm->barrier_and_abort("VMC batched driver is not supported by legacy CUDA builds.");
#endif

  app_summary() << "\n========================================"
                   "\n  Reading VMC driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(input_node_);
  VMCDriverInput vmcdriver_input;
  vmcdriver_input.readXML(input_node_);
  QMCDriverInterface* qmc = nullptr;

  if (vmc_mode_ == 0 || vmc_mode_ == 1) //(0,0,0) (0,0,1)
  {
    qmc = new VMCBatched(project_data, std::move(qmcdriver_input), std::move(vmcdriver_input), std::move(pop), samples,
                         comm);
  }
  else
  {
    throw std::runtime_error("VMCFactoryNew does not support VMC_MODE");
  }

  // TODO: I believe this is redundant and could become a bug.
  qmc->setUpdateMode(vmc_mode_ & 1);
  return qmc;
}
} // namespace qmcplusplus
