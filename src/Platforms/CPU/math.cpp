//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


// To make the customized qmcplusplus::isnan/isfinite/isinf always effective, this file must be compiled without -ffast-math.
#include <cmath>

namespace qmcplusplus
{
bool isnan(float a) { return a != a; }
bool isnan(double a) { return a != a; }

bool isfinite(float a) { return std::isfinite(a); }
bool isfinite(double a) { return std::isfinite(a); }

bool isinf(float a) { return std::isinf(a); }
bool isinf(double a) { return std::isinf(a); }
} // namespace qmcplusplus
