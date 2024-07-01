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

#include "DiracMatrixComputeOMPTarget.hpp"

namespace qmcplusplus
{
template class DiracMatrixComputeOMPTarget<double, float>;
template class DiracMatrixComputeOMPTarget<double, double>;
template class DiracMatrixComputeOMPTarget<std::complex<double>, std::complex<float>>;
template class DiracMatrixComputeOMPTarget<std::complex<double>, std::complex<double>>;
} // namespace qmcplusplus
