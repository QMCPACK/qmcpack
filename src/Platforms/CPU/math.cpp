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

namespace qmcplusplus
{
// To make the customized qmcplusplus::isnan always effective, this file must be compiled without -ffast-math.
bool isnan(float a) { return a != a; }
bool isnan(double a) { return a != a; }
} // namespace qmcplusplus
