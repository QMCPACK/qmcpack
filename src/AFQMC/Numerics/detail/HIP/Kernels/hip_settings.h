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

#ifndef AFQMC_KERNELS_SETTINGS_HPP
#define AFQMC_KERNELS_SETTINGS_HPP

#define BOOST_NO_AUTO_PTR

static const size_t DEFAULT_BLOCK_SIZE  = 32;
static const size_t DOT_BLOCK_SIZE      = 32;
static const size_t REDUCE_BLOCK_SIZE   = 32;
static const size_t MAX_THREADS_PER_DIM = 1024;

#endif
