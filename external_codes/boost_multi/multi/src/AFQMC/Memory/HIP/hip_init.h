///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_HIP_INIT_H
#define AFQMC_HIP_INIT_H
#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"

namespace qmc_hip
{
void HIP_INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed = 911ULL);

}

#endif
