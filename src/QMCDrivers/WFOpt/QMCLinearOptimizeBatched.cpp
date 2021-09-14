//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCLinearOptimizeBatched.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>


namespace qmcplusplus
{
QMCLinearOptimizeBatched::QMCLinearOptimizeBatched(const ProjectData& project_data,
                                                   MCWalkerConfiguration& w,
                                                   QMCDriverInput&& qmcdriver_input,
                                                   VMCDriverInput&& vmcdriver_input,
                                                   MCPopulation&& population,
                                                   SampleStack& samples,
                                                   Communicate* comm,
                                                   const std::string& QMC_driver_type)
    : QMCDriverNew(project_data,
                   std::move(qmcdriver_input),
                   std::move(population),
                   "QMCLinearOptimizeBatched::",
                   comm,
                   "QMCLinearOptimizeBatched"),
      param_tol(1e-4),
      vmcdriver_input_(vmcdriver_input),
      samples_(samples)
{
  //     //set the optimization flag
  qmc_driver_mode_.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  m_param.add(param_tol, "alloweddifference");
  //Set parameters for line minimization:
}

} // namespace qmcplusplus
