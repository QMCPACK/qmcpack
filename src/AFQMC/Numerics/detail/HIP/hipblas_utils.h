//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_HIPBLAS_UTILS_H
#define AFQMC_HIPBLAS_UTILS_H

#include "hipblas.h"
#include "rocblas.h"

namespace hipblas {

  // TODO: Temporary hack waiting for upstream version of hipblas
  hipblasStatus_t rocBLASStatusToHIPStatusAFQMC(rocblas_status_ error);

}

#endif
