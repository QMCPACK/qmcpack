//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_CUDA_INIT_HPP
#define AFQMC_CUDA_INIT_HPP

#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"

namespace qmc_cuda
{
void CUDA_INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed = 911ULL);

}

#endif
