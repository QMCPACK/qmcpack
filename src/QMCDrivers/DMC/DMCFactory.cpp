//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "DMCFactory.h"
#include "QMCDrivers/DMC/DMC.h"
#include "Concurrency/OpenMP.h"

#ifdef QMC_CUDA
#include "QMCDrivers/DMC/DMC_CUDA.h"
#endif

//#define PETA_DMC_TEST
namespace qmcplusplus
{
std::unique_ptr<QMCDriver> DMCFactory::create(const ProjectData& project_data,
                                              MCWalkerConfiguration& w,
                                              TrialWaveFunction& psi,
                                              QMCHamiltonian& h,
                                              Communicate* comm,
                                              bool enable_profiling)
{
#ifdef QMC_CUDA
  if (GPU)
    return std::make_unique<DMCcuda>(project_data, w, psi, h, comm, enable_profiling);
#endif
  auto qmc = std::make_unique<DMC>(project_data, w, psi, h, comm, enable_profiling);
  qmc->setUpdateMode(PbyPUpdate);
  return qmc;
}
} // namespace qmcplusplus
