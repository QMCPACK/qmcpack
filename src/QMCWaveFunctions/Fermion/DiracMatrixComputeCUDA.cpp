//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "DiracMatrixComputeCUDA.hpp"

namespace qmcplusplus
{
template class DiracMatrixComputeCUDA<double, float>;
template class DiracMatrixComputeCUDA<double, double>;
template class DiracMatrixComputeCUDA<std::complex<double>, std::complex<float>>;
template class DiracMatrixComputeCUDA<std::complex<double>, std::complex<double>>;
} // namespace qmcplusplus
