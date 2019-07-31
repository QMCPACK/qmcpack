//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file DMC.h
 * @brief Define OpenMP-able DMC Driver.
 */
#ifndef QMCPLUSPLUS_DMC_H
#define QMCPLUSPLUS_DMC_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief DMC driver using OpenMP paragra
 *
 * This is the main DMC driver with MPI/OpenMP loops over the walkers.
 */
class DMC : public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  DMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool, Communicate* comm);

  bool run();
  bool put(xmlNodePtr cur);
  void setTau(RealType i);
  QMCRunType getRunType() { return QMCRunType::DMC; }

private:
  ///Index to determine what to do when node crossing is detected
  IndexType KillNodeCrossing;
  ///Interval between branching
  IndexType BranchInterval;
  ///hdf5 file name for Branch conditions
  std::string BranchInfo;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  ///input std::string to determine to use reconfiguration
  std::string Reconfiguration;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  ///input std::string to use fast gradient
  std::string UseFastGrad;
  ///input to control maximum age allowed for walkers.
  IndexType mover_MaxAge;

  void resetUpdateEngines();
  /// Copy Constructor (disabled)
  DMC(const DMC&) = delete;
  /// Copy operator (disabled).
  DMC& operator=(const DMC&) = delete;
};
} // namespace qmcplusplus

#endif
