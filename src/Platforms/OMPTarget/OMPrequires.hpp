//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OMPrequires.hpp
 */

#ifndef QMCPLUSPLUS_OPENMP_REQUIRES_H
#define QMCPLUSPLUS_OPENMP_REQUIRES_H

#ifdef QMC_OFFLOAD_USM
#pragma omp requires unified_shared_memory
#endif

#endif
