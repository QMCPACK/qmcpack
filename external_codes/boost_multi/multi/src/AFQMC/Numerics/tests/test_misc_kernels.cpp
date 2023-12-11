//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Fionn D. Malone, malone14@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    exit(0);         \
  }

#include <vector>


#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#if defined(ENABLE_CUDA)
#include "AFQMC/Numerics/detail/CUDA/blas_cuda_gpu_ptr.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"
#elif defined(ENABLE_HIP)
#include "AFQMC/Numerics/detail/HIP/blas_hip_gpu_ptr.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"
#endif

#include "multi/array.hpp"
#include "multi/array_ref.hpp"


using boost::multi::array;
using boost::multi::array_ref;
using boost::multi::iextensions;
using std::copy_n;

namespace qmcplusplus
{
using namespace afqmc;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
template<typename T>
using Alloc = device::device_allocator<T>;
#else
template<typename T>
using Alloc = std::allocator<T>;
#endif
template<typename T>
using pointer = typename Alloc<T>::pointer;


template<typename T>
using Tensor1D = array<T, 1, Alloc<T>>;
template<typename T>
using Tensor2D = array<T, 2, Alloc<T>>;
template<typename T>
using Tensor3D = array<T, 3, Alloc<T>>;

TEST_CASE("axpyBatched", "[Numerics][misc_kernels]")
{
  // Only implemented for complex
  Alloc<std::complex<double>> alloc{};
  Tensor2D<std::complex<double>> y({3, 4}, 1.0, alloc);
  Tensor2D<std::complex<double>> x({3, 4}, 1.0, alloc);
  Tensor1D<std::complex<double>> a(iextensions<1u>{3}, 2.0, alloc);
  std::vector<pointer<std::complex<double>>> x_batched, y_batched;
  for (int i = 0; i < std::get<0>(x.sizes()); i++)
  {
    x_batched.emplace_back(x[i].origin());
    y_batched.emplace_back(y[i].origin());
  }
  using ma::axpyBatched;
  axpyBatched(std::get<1>(x.sizes()), to_address(a.origin()), x_batched.data(), 1, y_batched.data(), 1, x_batched.size());
  // 1 + 2 = 3.
  Tensor2D<std::complex<double>> ref({3, 4}, 3.0, alloc);
  verify_approx(y, ref);
}

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
// Not in dispatching routine yet, called directly from AFQMCBasePropagator.
TEST_CASE("construct_X", "[Numerics][misc_kernels]")
{
  Alloc<std::complex<double>> alloc{};
  int ncv                 = 11;
  int nsteps              = 2;
  int nwalk               = 3;
  bool fp                 = false;
  double sqrtdt           = 0.002;
  double vbound           = 40.0;
  std::complex<double> im = std::complex<double>(0.0, 1.0);
  Tensor1D<std::complex<double>> vmf(iextensions<1U>{ncv}, im, alloc);
  Tensor2D<std::complex<double>> vbias({ncv, nwalk}, 1.0, alloc);
  Tensor2D<std::complex<double>> hws({nsteps, nwalk}, -0.2, alloc);
  Tensor2D<std::complex<double>> mf({nsteps, nwalk}, 2.0, alloc);
  Tensor3D<std::complex<double>> x({ncv, nsteps, nwalk}, 0.1, alloc);
  using kernels::construct_X;
  construct_X(ncv, nsteps, nwalk, fp, sqrtdt, vbound, to_address(vmf.origin()), to_address(vbias.origin()),
              to_address(hws.origin()), to_address(mf.origin()), to_address(x.origin()));
  // captured from stdout.
  std::complex<double> ref_val = std::complex<double>(0.102, 0.08);
  Tensor3D<std::complex<double>> ref({ncv, nsteps, nwalk}, ref_val, alloc);
  verify_approx(ref, x);
}

// No cpu equivalent?
TEST_CASE("batchedDot", "[Numerics][misc_kernels]")
{
  Alloc<std::complex<double>> alloc{};
  std::complex<double> im = std::complex<double>(0.0, 1.0);
  int dim                 = 3;
  Tensor1D<std::complex<double>> y(iextensions<1U>{dim}, im, alloc);
  Tensor2D<std::complex<double>> A({dim, dim}, 1.0, alloc);
  Tensor2D<std::complex<double>> B({dim, dim}, -0.2, alloc);
  std::complex<double> alpha(2.0);
  std::complex<double> beta(-1.0);
  using kernels::batchedDot;
  batchedDot(dim, dim, alpha, to_address(A.origin()), dim, to_address(B.origin()), dim, beta, to_address(y.origin()),
             1);
  // from numpy.
  std::complex<double> ref_val(-1.2, -1.0);
  Tensor1D<std::complex<double>> ref(iextensions<1U>{dim}, ref_val, alloc);
  verify_approx(ref, y);
}
#endif

} // namespace qmcplusplus
