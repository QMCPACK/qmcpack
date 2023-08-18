//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "VMCFactory.h"
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#include "RandomNumberControl.h"
#if defined(QMC_BUILD_COMPLETE)
//REMOVE Broken warping
//#if !defined(QMC_COMPLEX)
//#include "QMCDrivers/VMC/VMCMultipleWarp.h"
//#include "QMCDrivers/VMC/VMCPbyPMultiWarp.h"
//#endif
//#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#endif
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
std::unique_ptr<QMCDriverInterface> VMCFactory::create(const ProjectData& project_data,
                                                       MCWalkerConfiguration& w,
                                                       TrialWaveFunction& psi,
                                                       QMCHamiltonian& h,
                                                       Communicate* comm,
                                                       bool enable_profiling)
{
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  std::unique_ptr<QMCDriverInterface> qmc;
  if (VMCMode == 0 || VMCMode == 1) //(0,0,0) (0,0,1)
  {
    qmc = std::make_unique<VMC>(project_data, w, psi, h, RandomNumberControl::Children, comm, enable_profiling);
  }
  else if (VMCMode == 2 || VMCMode == 3)
  {
    qmc = std::make_unique<CSVMC>(project_data, w, psi, h, comm);
  }
  qmc->setUpdateMode(VMCMode & 1);
  return qmc;
}
} // namespace qmcplusplus
