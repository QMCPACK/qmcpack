//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RUNTIMEOPTIONS_H__
#define QMCPLUSPLUS_RUNTIMEOPTIONS_H__

namespace qmcplusplus
{

// Simple struct that tracks (@TODO: future) runtime options for the project:
// - mixed precision,
// - complex
struct RuntimeOptions
{
  bool is_complex;
  bool is_mixed_precision;

  // Default constructor follows current compile-time options
  // @TODO: remove this at the end of the complex/mixed precision integrations
  RuntimeOptions()
      : is_complex(
#ifdef QMC_COMPLEX
            true
#else
            false
#endif
            ),
        is_mixed_precision(
#ifdef QMC_MIXED_PRECISION
            true
#else
            false
#endif
        )
  {}
};

} // namespace qmcplusplus

#endif
