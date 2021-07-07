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

#include "makeRngSpdMatrix.hpp"

namespace qmcplusplus
{
namespace testing
{
template void makeRngSpdMatrix<double>(Matrix<double>& mat_spd);
template void makeRngSpdMatrix<float>(Matrix<float>& mat_spd);
template void makeRngSpdMatrix<std::complex<double>>(Matrix<std::complex<double>>& mat_spd);
template void makeRngSpdMatrix<std::complex<float>>(Matrix<std::complex<float>>& mat_spd);
} // namespace testing
} // namespace qmcplusplus
