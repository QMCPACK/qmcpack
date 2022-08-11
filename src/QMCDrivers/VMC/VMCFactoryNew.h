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

#ifndef QMCPLUSPLUS_VMCFACTORYNEW_H
#define QMCPLUSPLUS_VMCFACTORYNEW_H
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Message/Communicate.h"


namespace qmcplusplus
{
class ParticleSetPool;
class HamiltonianPool;
class MCPopulation;
class ProjectData;

class VMCFactoryNew
{
private:
  const int vmc_mode_;
  // const ?
  xmlNodePtr input_node_;

public:
  VMCFactoryNew(xmlNodePtr cur, const int vmode) : vmc_mode_(vmode), input_node_(cur) {}

  /** create a VMCBatched driver.
   *  \param[in]   project_data   containing so basic options including DriverVersion and max_cpu_seconds
   *  \param[in]   global_emi     optional global estimator manager input passed by value to insure copy,
   *                              a global input should not be consumed by driver.
   */
  std::unique_ptr<QMCDriverInterface> create(const ProjectData& project_data,
                                             const std::optional<EstimatorManagerInput>& global_emi,
                                             WalkerConfigurations& wc,
                                             MCPopulation&& pop,
                                             SampleStack& samples,
                                             Communicate* comm);
};
} // namespace qmcplusplus

#endif
