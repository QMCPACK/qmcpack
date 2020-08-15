////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_HSPOTENTIAL_HELPERS_H
#define AFQMC_HSPOTENTIAL_HELPERS_H

#include <vector>
#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"

namespace qmcplusplus
{
namespace afqmc
{
namespace HamHelper
{
std::vector<std::size_t> count_nnz_per_cholvec(double cut, TaskGroup_& TG, SpVType_shm_csr_matrix& V2, int NMO);

std::vector<std::size_t> count_nnz_per_ik(double cut,
                                          TaskGroup_& TG,
                                          SpVType_shm_csr_matrix& V2,
                                          int NMO,
                                          int cv0,
                                          int cvN);

void generateHSPotential(SpVType_shm_csr_matrix& vn,
                         std::vector<int> const& map_,
                         double cut,
                         TaskGroup_& TG,
                         SpVType_shm_csr_matrix& V2,
                         int NMO,
                         int cv0,
                         int cvN);

} // namespace HamHelper

} // namespace afqmc
} // namespace qmcplusplus
#endif
