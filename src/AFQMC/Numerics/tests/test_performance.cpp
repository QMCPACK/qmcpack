///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////


//#include "catch.hpp"
#include "Configuration.h"

#include <vector>

#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Utilities/myTimer.h"

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

// numpy.arange(size).reshape() / scale
template<typename T>
void fill_matrix(Tensor2D<T> m, T scale)
{
  T count = T(0);
  for (int i = 0; i < m.size(0); i++) {
    for (int j = 0; j < m.size(1); j++) {
      m[i][j] = count;
      count += T(1) / scale;
    }
  }
}

template<typename T>
void time_qr(int nbatch, int m, int n)
{
  Alloc<T> alloc{};
  Tensor3D<T> A({nbatch, m, n}, alloc);

  Tensor3D<T> AT({nbatch, n, m}, alloc);
  Tensor2D<T> T_({nbatch, m}, alloc);
  Tensor2D<T> scl({nbatch, m}, alloc);
  int sz = ma::gqr_optimal_workspace_size(AT[0]);
  Tensor1D<T>  WORK(iextensions<1u>{nbatch * sz}, alloc);
  Alloc<int> ialloc{};
  std::vector<pointer<T>> Aarray;
  Tensor1D<int> IWORK(iextensions<1u>{nbatch * (m + 1)}, ialloc);
  using std::copy_n;
  myTimer timer;
  for (int i = 0; i < nbatch; i++) {
    Tensor2D<T> tmp({m,n}, alloc);
    fill_matrix(tmp, T(100.0));
    copy_n(tmp.origin(), tmp.num_elements(), A[i].origin());
    Aarray.emplace_back(A[i].origin());
  }
  // Actual profile.
  std::string timer_id = std::string("transpose");
  timer.add(timer_id);
  timer.start(timer_id);
  for (int i = 0; i < nbatch; i++)
    ma::transpose(A[i], AT[i]);
  auto ttrans = timer.elapsed(timer_id);
  // careful, expects fortran order
  timer_id = std::string("geqrf");
  timer.add(timer_id);
  timer.start(timer_id);
  geqrfStrided(m, n, AT.origin(), m, m*n, T_.origin(), m, IWORK.origin(), nbatch);
  auto tgeqrf = timer.elapsed(timer_id);
  //timer.stop("transpose_"+timer_id);
  //using ma::determinant_from_geqrf;
  //using ma::scale_columns;
  timer_id = std::string("det_from_geqrf");
  timer.add(timer_id);
  timer.start(timer_id);
  for (int i = 0; i < nbatch; i++)
    determinant_from_geqrf(n, AT[i].origin(), m, scl[i].origin(), T(0.0));
  auto tdet = timer.elapsed(timer_id);
  timer_id = std::string("gqr");
  timer.add(timer_id);
  timer.start(timer_id);
  gqrStrided(m, n, n, AT.origin(), m, m * n, T_.origin(), m, WORK.origin(), sz, IWORK.origin(),
             nbatch);
  auto tgqr = timer.elapsed(timer_id);
  std::cout << nbatch << " " << m << " " << n << " " << ttrans << " " << tgeqrf << " " << tdet << " " << tgqr << std::endl;
  //for (int i = 0; i < nbatch; i++)
  //{
    //ma::transpose(AT[i], *Aarray[i]);
    //scale_columns(NMO, NAEA, (*Aarray[i]).origin(), (*Aarray[i]).stride(0), scl[i].origin());
  //}
}

TEST_CASE("test_performance_qr")
{
  // QR mxn matrix
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
  arch::INIT(node);
  std::vector<int> ms = {312, 468, 702, 1053, 1404, 1872}; 
  std::vector<int> ns = {96, 144, 216, 325, 432, 576};
  time_qr<std::complex<double>>(1, ms[3], ns[3]);
  //for (int nb = 2; nb < 32; nb+=4) {
    for (int i = 0; i < ms.size(); i++) {
      time_qr<std::complex<double>>(5, ms[3], ns[3]);
    }
    //std::cout << std::endl;
  //}
}


TEST_CASE("test_performance_batched_gemm")
{
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
  arch::INIT(node);
  Alloc<float> alloc{};
  int nkmax = 125;
  int nkmin = 8;
  int nbatch_max = nkmax*nkmax;
  int m = 300;
  int nbatch = nkmin*nkmin;
  myTimer timer;
  while (nbatch < nbatch_max) {
    std::string timer_id = std::to_string(nbatch);
    timer.add(timer_id);
    Tensor2D<float> a({m, m}, alloc);
    Tensor2D<float> b({m, m}, alloc);
    Tensor3D<float> c({nbatch, m, m}, 0.0, alloc);
    float scale = float(100.0);
    fill_matrix(a, scale);
    fill_matrix(b, scale);
    std::vector<pointer<float>> A_array;
    std::vector<pointer<float>> B_array;
    std::vector<pointer<float>> C_array;
    float alpha = 1.0;
    float beta  = 0.0;
    for (int i = 0; i < nbatch; i++)
    {
      A_array.emplace_back(a.origin());
      B_array.emplace_back(b.origin());
      C_array.emplace_back(c[i].origin());
    }
    using ma::gemmBatched;
    timer.start(timer_id);
    gemmBatched('N', 'N', m, m, m, alpha, A_array.data(), m, B_array.data(), m, beta, C_array.data(), m, nbatch);
    timer.stop(timer_id);
    std::cout << nbatch << " " << timer.elapsed(timer_id) << std::endl;
    nbatch *= 2;
  }
}
