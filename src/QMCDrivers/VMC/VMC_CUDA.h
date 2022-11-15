//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VMC_CUDA_H
#define QMCPLUSPLUS_VMC_CUDA_H
#include "QMCDrivers/QMCDriver.h"
#include "type_traits/CUDATypes.h"
namespace qmcplusplus
{
class QMCUpdateBase;

/** @ingroup QMCDrivers  PbyP
 *@brief Implements the VMC algorithm using particle-by-particle move.
 */
class VMCcuda : public QMCDriver
{
public:
  /// Constructor.
  VMCcuda(const ProjectData& project_data,
          MCWalkerConfiguration& w,
          TrialWaveFunction& psi,
          QMCHamiltonian& h,
          Communicate* comm,
          bool enable_profiling);

  bool run() override;
  bool runWithDrift();

  /// advance walkers without drift
  void advanceWalkers();
  /// advance walkers with drift
  void advanceWalkersWithDrift();

  bool put(xmlNodePtr cur) override;

private:
  using CTS = CUDAGlobalTypes;
  ///previous steps
  int prevSteps;
  ///previous stepsbetweensamples
  int prevStepsBetweenSamples;
  /// tau/mass
  RealType m_tauovermass;
  /// Whether to use drift or not
  std::string UseDrift;
  ///period for walker dump
  int myPeriod4WalkerDump;
  /// Copy Constructor (disabled)
  VMCcuda(const VMCcuda&) = delete;
  /// Copy operator (disabled).
  VMCcuda& operator=(const VMCcuda&) = delete;
  ///hide initialization from the main function
  bool checkBounds(std::vector<PosType>& newpos, std::vector<bool>& valid);
  QMCRunType getRunType() override { return QMCRunType::VMC; }
  void resetRun();
};
} // namespace qmcplusplus

#endif
