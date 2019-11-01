//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
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
#include "QMCApp/ParticleSetPool.h"

class Communicate;

namespace qmcplusplus
{
//forward declaration
class MCWalkerConfiguration;
class QMCDriverInterface;
class WaveFunctionPool;
class HamiltonianPool;

class QMCDriverFactory
{
public:
  struct DriverAssemblyState
  {
    std::bitset<QMC_MODE_MAX> what_to_do;
    bool append_run         = false;
    std::string traces_tag  = "none";
    QMCRunType new_run_type = QMCRunType::DUMMY;
  };

  /** default constructor **/
  //QMCDriverFactory() ;

  /** read the current QMC Section */
  DriverAssemblyState readSection(int curSeries, xmlNodePtr cur);

  /** set the active qmcDriver
   *  The unique_ptr's take care of killing the old driver if it exists */
  std::unique_ptr<QMCDriverInterface> newQMCDriver(std::unique_ptr<QMCDriverInterface> last_driver,
                                                   int curSeries,
                                                   xmlNodePtr cur,
                                                   DriverAssemblyState& das,
                                                   MCWalkerConfiguration& qmc_system,
                                                   ParticleSetPool& particle_pool,
                                                   WaveFunctionPool& wave_function_pool,
                                                   HamiltonianPool& hamiltonian_pool,
                                                   MCPopulation& population,
                                                   Communicate* comm);

private:
  /** create a new QMCDriver
   *
   *  Broken out for unit tests
   */
  std::unique_ptr<QMCDriverInterface> createQMCDriver(xmlNodePtr cur,
                                                      DriverAssemblyState& das,
                                                      MCWalkerConfiguration& qmc_system,
                                                      ParticleSetPool& particle_pool,
                                                      WaveFunctionPool& wave_function_pool,
                                                      HamiltonianPool& hamiltonian_pool,
                                                      MCPopulation& population,
                                                      Communicate* comm);
};
} // namespace qmcplusplus
#endif
