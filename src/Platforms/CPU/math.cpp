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


// To make the customized qmcplusplus::numeric_limits/isnan/isfinite/isinf always effective,
// this file must be compiled without -ffast-math.
#include "math.hpp"
#include <cmath>
#include <limits>

namespace qmcplusplus
{
float numeric_limits<float>::infinity() { return std::numeric_limits<float>::infinity(); }
float numeric_limits<float>::quiet_NaN() { return std::numeric_limits<float>::quiet_NaN(); }

double numeric_limits<double>::infinity() { return std::numeric_limits<double>::infinity(); }
double numeric_limits<double>::quiet_NaN() { return std::numeric_limits<double>::quiet_NaN(); }

bool isnan(float a) { return std::isnan(a); }
bool isnan(double a) { return std::isnan(a); }

bool isfinite(float a) { return std::isfinite(a); }
bool isfinite(double a) { return std::isfinite(a); }

bool isinf(float a) { return std::isinf(a); }
bool isinf(double a) { return std::isinf(a); }
} // namespace qmcplusplus
