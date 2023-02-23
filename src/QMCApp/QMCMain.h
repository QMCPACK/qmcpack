//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCMain.h
 * @brief Declaration of QMCMain class
 */
#ifndef QMCPLUSPLUS_MAINAPPLICATIONS_H
#define QMCPLUSPLUS_MAINAPPLICATIONS_H

#include "Message/MPIObjectBase.h"
#include "QMCApp/QMCAppBase.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "hdf/hdf_error_suppression.h"

namespace qmcplusplus
{
/** @ingroup qmcapp
 * @brief Main application to perform QMC simulations
 *
 * This is a generalized QMC application which can handle multiple ParticleSet,
 * TrialWaveFunction and QMCHamiltonian objects.
 */
class QMCMain : public MPIObjectBase, public QMCAppBase
{
public:
  QMCMain(Communicate* c);

  ~QMCMain() override;

  bool validateXML() override;
  bool execute() override;

  ParticleSetPool& getParticlePool() { return *particle_set_pool_; }

private:
  /// ParticleSet Pool
  std::unique_ptr<ParticleSetPool> particle_set_pool_;

  /// TrialWaveFunction Pool
  std::unique_ptr<WaveFunctionPool> psi_pool_;

  /// QMCHamiltonian Pool
  std::unique_ptr<HamiltonianPool> ham_pool_;

  /// current MCWalkerConfiguration
  MCWalkerConfiguration* qmc_system_;

  ///Global estimators defined outside of <qmc> nodes
  std::optional<EstimatorManagerInput> estimator_manager_input_;

  ///flag to indicate that a qmc is the first QMC
  bool first_qmc_;

  /// the last driver object. Should be in a loop only.
  std::unique_ptr<QMCDriverInterface> last_driver_;

  /// last branch engine used by legacy drivers
  std::unique_ptr<SimpleFixedNodeBranch> last_branch_engine_legacy_driver_;

  /// Suppress HDF5 warning and error messages.
  hdf_error_suppression hide_hdf_errors;

  ///xml mcwalkerset elements for output
  std::vector<xmlNodePtr> walker_set_;

  ///xml mcwalkerset read-in elements
  std::vector<xmlNodePtr> walker_set_in_;

  ///traces xml
  xmlNodePtr traces_xml_;

  ///qmc sections
  std::vector<std::pair<xmlNodePtr, bool>> qmc_action_;

  ///pointer to the last node of the main inputfile
  xmlNodePtr last_input_node_;

  /** execute <qmc/> element
   * @param cur qmc xml node
   * @param reuse if true, reuse the driver built from the last QMC section. This should be used by loop only.
   *
   * Ye: I think this can be merged with executeQMCSection
   */
  bool runQMC(xmlNodePtr cur, bool reuse);

  ///add <mcwalkerset/> elements to continue a run
  bool setMCWalkers(xmlXPathContextPtr cur);

  /** add unique ParticleSet, TrialWaveFunction and QMCHamiltonian elements to Pool objects
   * @param cur xml node
   * @return false, if contains qmc actions
   */
  bool processPWH(xmlNodePtr cur);

  /** execute loop **/
  void executeLoop(xmlNodePtr cur);

  ///execute <debug/> element
  bool executeDebugSection(xmlNodePtr cur);

  /** execute qmc
   * @param cur qmc xml node
   * @param reuse if true, reuse the driver built from the last QMC section. This should be used by loop only.
   * @return true, if a section is successfully executed.
   */
  bool executeQMCSection(xmlNodePtr cur, bool reuse = false);

  ///execute <cmc/> element
  bool executeCMCSection(xmlNodePtr cur);
};
} // namespace qmcplusplus
#endif
