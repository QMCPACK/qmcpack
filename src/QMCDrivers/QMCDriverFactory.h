//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#define QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#include <bitset>
#include <string>

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCDrivers/DriverTraits.h"
#include "QMCDrivers/MCPopulation.h"
#include "Particle/ParticleSetPool.h"
#include "Estimators/EstimatorManagerInput.h"

class Communicate;

namespace qmcplusplus
{
//forward declaration
class MCWalkerConfiguration;
class QMCDriverInterface;
class WaveFunctionPool;
class HamiltonianPool;
class ProjectData;

class QMCDriverFactory
{
public:
  struct DriverAssemblyState
  {
    std::bitset<QMC_MODE_MAX> what_to_do;
    bool append_run         = false;
    bool enable_profiling   = false;
    std::string traces_tag  = "none";
    QMCRunType new_run_type = QMCRunType::DUMMY;
  };

  /** Application uses this constructor
   *  param[in] project_data    this is stored as a reference and this state controls later behavior.
   *                            For both the driver factory i.e. driver verion. And the drivers it creates
   *                            i.e. the section id and max CPU seconds.
   */
  QMCDriverFactory(const ProjectData& project_data);

  /** default constructor **/
  //QMCDriverFactory() ;

  /** read the current QMC Section
   *  In the application context project data can indicate the input be read in the context of
   *  the batched driver architecture.
   *  param[in] cur            qmc section node
   *  param[in] emi            std::optional<EstimatorManagerInput> if it is there it is the global estimator manager input.
   */
  DriverAssemblyState readSection(xmlNodePtr cur) const;

  /** create a new QMCDriver
   *
   *  Broken out for unit tests
   */
  std::unique_ptr<QMCDriverInterface> createQMCDriver(xmlNodePtr cur,
                                                      DriverAssemblyState& das,
                                                      const std::optional<EstimatorManagerInput>& emi,
                                                      MCWalkerConfiguration& qmc_system,
                                                      ParticleSetPool& particle_pool,
                                                      WaveFunctionPool& wave_function_pool,
                                                      HamiltonianPool& hamiltonian_pool,
                                                      Communicate* comm) const;

private:
  /// project info for accessing global fileroot and series id
  const ProjectData& project_data_;
};
} // namespace qmcplusplus
#endif
