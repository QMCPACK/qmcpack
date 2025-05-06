//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "Common/Queue.hpp"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/InverterAccel.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Containers/tests/makeRngSpdMatrix.hpp"

namespace qmcplusplus
{

template<PlatformKind PL, typename T>
using DeviceAllocator = typename compute::MemManage<PL>::template DeviceAllocator<T>;

template<PlatformKind PL, typename T, typename FP_T>
void benchmark_solver(const size_t N)
{
  typename InverterAccel<PL, FP_T>::Inverter solver;
  compute::Queue<PL> queue;

  Matrix<T> m(N, N);
  Matrix<T> m_invT(N, N);
  Matrix<T> m_invT_CPU(N, N);
  Matrix<T, DeviceAllocator<PL, T>> m_invGPU;
  std::complex<double> log_value;
  m.resize(N, N);
  m_invT.resize(N, N);
  m_invT_CPU.resize(N, N);
  m_invGPU.resize(N, N);

  testing::MakeRngSpdMatrix<T> makeRngSpdMatrix{};
  makeRngSpdMatrix(m);

  BENCHMARK("SolverInverter") { solver.invert_transpose(m, m_invT, m_invGPU, log_value, queue.getNative()); };

  DiracMatrix<T> dmat;
  BENCHMARK("CPU") { dmat.invert_transpose(m, m_invT_CPU, log_value); };

  auto check_matrix_result = checkMatrix(m_invT, m_invT_CPU);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

TEST_CASE("SolverInverter_bench", "[wavefunction][benchmark]")
{
  using Value         = typename QMCTraits::ValueType;
#ifdef QMC_COMPLEX
  using FullPrecValue = std::complex<double>;
#else
  using FullPrecValue = double;
#endif
  benchmark_solver<PlatformKind::CUDA, Value, FullPrecValue>(1024);
}

} // namespace qmcplusplus
