//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERSET_HPP
#define QMCPLUSPLUS_AFQMC_WALKERSET_HPP

#include "AFQMC/Walkers/SharedWalkerSet.hpp"
#include "AFQMC/Walkers/SerialWalkerSet.hpp"

namespace qmcplusplus
{
namespace afqmc
{
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
using WalkerSet = SerialWalkerSet;
#else
using WalkerSet = SharedWalkerSet;
#endif

} // namespace afqmc

} // namespace qmcplusplus

#endif
