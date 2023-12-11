//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "RMCFactory.h"
#include "QMCDrivers/RMC/RMC.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
std::unique_ptr<QMCDriver> RMCFactory::create(const ProjectData& project_data,
                                              MCWalkerConfiguration& w,
                                              TrialWaveFunction& psi,
                                              QMCHamiltonian& h,
                                              Communicate* comm)
{
  std::unique_ptr<QMCDriver> qmc;

  if (RMCMode == 0 || RMCMode == 1) //(0,0,0) (0,0,1) pbyp and all electron
  {
    qmc = std::make_unique<RMC>(project_data, w, psi, h, comm);
  }
  qmc->setUpdateMode(RMCMode & 1);
  return qmc;
}
} // namespace qmcplusplus
