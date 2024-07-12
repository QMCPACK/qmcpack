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

#include "DiracMatrixInverterOMPTarget.hpp"

namespace qmcplusplus
{
template class DiracMatrixInverterOMPTarget<double, float>;
template class DiracMatrixInverterOMPTarget<double, double>;
template class DiracMatrixInverterOMPTarget<std::complex<double>, std::complex<float>>;
template class DiracMatrixInverterOMPTarget<std::complex<double>, std::complex<double>>;
} // namespace qmcplusplus
