//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  \brief includes catch.hpp with the correct defines for benchmarking.
 *  If benchmarking macros are going to be used CATCH_CONFIG_ENABLE_BENCHMARKING must be defined
 *  before catch.hpp is included in the file using them.
 *  This could also be done in the cmake.
 */

#ifndef QMCPLUSPLUS_CATCH_BENCHMARK_MAIN_HPP
#define QMCPLUSPLUS_CATCH_BENCHMARK_MAIN_HPP

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#endif
