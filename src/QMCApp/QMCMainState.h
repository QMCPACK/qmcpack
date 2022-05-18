//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// Refactored from QMCDriverFactory.h  -- created by: Jeongnim Kim
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMCMAINSTATE_H
#define QMCPLUSPLUS_QMCMAINSTATE_H
#include <bitset>
#include <string>

#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "Particle/ParticleSetPool.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "Estimators/EstimatorManagerInput.h"
#include "type_traits/template_types.hpp"

class Communicate;

namespace qmcplusplus
{
//forward declaration
class MCWalkerConfiguration;
class QMCDriverInterface;
class WaveFunctionPool;
class HamiltonianPool;

struct QMCMainState : public MPIObjectBase
{
  ///type of qmcdriver
  QMCRunType curRunType;

  ///name of the current QMCriver
  std::string curMethod;

  /** current MCWalkerConfiguration
   */
  MCWalkerConfiguration* qmcSystem;

  /** current QMCDriver
   */
  std::unique_ptr<QMCDriverInterface> last_driver;

  /** ParticleSet Pool
   */
  std::unique_ptr<ParticleSetPool> ptclPool;

  /** TrialWaveFunction Pool
   */
  std::unique_ptr<WaveFunctionPool> psiPool;

  /** QMCHamiltonian Pool
   */
  std::unique_ptr<HamiltonianPool> hamPool;

  /** default constructor **/
  QMCMainState(Communicate* c);

  /** set the active qmcDriver */
  void putCommunicator(xmlNodePtr cur);

  ///Global estimators defined outside of <qmc> nodes
  std::optional<EstimatorManagerInput> estimator_manager_input_;

  /** virtual destructor **/
  virtual ~QMCMainState();
};
} // namespace qmcplusplus
#endif
