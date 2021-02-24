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


#include "Configuration.h"

#include <vector>
#include <random>

#include "Utilities/Timer.h"

#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Memory/buffer_managers.h"
#include "AFQMC/Memory/arch.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

using namespace qmcplusplus;
using namespace afqmc;

using std::copy_n;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
template<typename T>
using Alloc = device::device_allocator<T>;
#else
template<typename T>
using Alloc = std::allocator<T>;
#endif
template<typename T>
using pointer = typename Alloc<T>::pointer;


//template<typename T>
//using buffer_alloc_type = device_buffer_type<T>;

//using buffer_ialloc_type = device_buffer_type<int>;

template<typename T>
using Tensor1D = boost::multi::array<T, 1, Alloc<T>>;
template<typename T>
using Tensor2D = boost::multi::array<T, 2, Alloc<T>>;
template<typename T>
using Tensor3D = boost::multi::array<T, 3, Alloc<T>>;

template<typename T>
using Tensor3D_ref = boost::multi::array_ref<T, 3, pointer<T>>;
template<typename T>
using Tensor2D_ref = boost::multi::array_ref<T, 2, pointer<T>>;
template<typename T>
using Tensor1D_ref = boost::multi::array_ref<T, 1, pointer<T>>;

template<typename T>
void fillRandomMatrix(std::vector<T>& vec)
{
  std::mt19937 generator(0);
  std::normal_distribution<T> distribution(0.0, 1.0);
  // avoid uninitialized warning
  T tmp = distribution(generator);
  for (int i = 0; i < vec.size(); i++)
  {
    T val  = distribution(generator);
    vec[i] = val;
  }
}

template<typename T>
void fillRandomMatrix(std::vector<std::complex<T>>& vec)
{
  std::mt19937 generator(0);
  std::normal_distribution<T> distribution(0.0, 1.0);
  T tmp = distribution(generator);
  for (int i = 0; i < vec.size(); i++)
  {
    T re   = distribution(generator);
    T im   = distribution(generator);
    vec[i] = std::complex<T>(re, im);
  }
}

template<class Allocator, class Buff>
void timeBatchedQR(std::ostream& out, Allocator& alloc, Buff& buffer, int nbatch, int m, int n)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor3D_ref<T> A(buffer.origin(), {nbatch, m, n});
  offset += A.num_elements();
  Tensor3D_ref<T> AT(buffer.origin() + offset, {nbatch, n, m});
  offset += AT.num_elements();
  Tensor2D_ref<T> T_(buffer.origin() + offset, {nbatch, m});
  offset += T_.num_elements();
  Tensor2D_ref<T> scl(buffer.origin() + offset, {nbatch, m});
  offset += T_.num_elements();
  int sz = ma::gqr_optimal_workspace_size(AT[0]);
  //std::cout << buffer.num_elements() << " " << 2*nbatch*m*n + 2*nbatch*m + nbatch*sz << " " << offset << std::endl;
  Tensor1D_ref<T> WORK(buffer.origin() + offset, boost::multi::iextensions<1u>{nbatch * sz});
  Alloc<int> ialloc{};
  std::vector<pointer<T>> Aarray;
  Tensor1D<int> IWORK(boost::multi::iextensions<1u>{nbatch * (m + 1)}, ialloc);
  using std::copy_n;
  for (int i = 0; i < nbatch; i++)
  {
    Aarray.emplace_back(A[i].origin());
  }
  // Actual profile.
  Timer timer;
  for (int i = 0; i < nbatch; i++)
    ma::transpose(A[i], AT[i]);
  double ttrans = timer.elapsed();
  timer.restart();
  geqrfStrided(m, n, AT.origin(), m, m * n, T_.origin(), m, IWORK.origin(), nbatch);
  double tgeqrf = timer.elapsed();
  timer.restart();
  for (int i = 0; i < nbatch; i++)
    determinant_from_geqrf(n, AT[i].origin(), m, scl[i].origin(), T(0.0));
  double tdet = timer.elapsed();
  timer.restart();
  gqrStrided(m, n, n, AT.origin(), m, m * n, T_.origin(), m, WORK.origin(), sz, IWORK.origin(), nbatch);
  double tgqr = timer.elapsed();
  out << "  " << std::setw(5) << nbatch << "   " << std::setw(5) << m << " " << std::setw(5) << n << " "
      << std::scientific << ttrans << " " << tgeqrf << " " << tdet << " " << tgqr << "\n";
}

template<class Allocator, class Buff>
void timeQR(std::ostream& out, Allocator& alloc, Buff& buffer, int m)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor2D_ref<T> A(buffer.origin(), {m, m});
  offset += A.num_elements();
  Tensor1D_ref<T> TAU(buffer.origin() + offset, {m});
  offset += TAU.num_elements();
  int sz = ma::gqr_optimal_workspace_size(A);
  Tensor1D_ref<T> WORK(buffer.origin() + offset, boost::multi::iextensions<1u>{sz});
  Timer timer;
  using ma::geqrf;
  geqrf(A, TAU, WORK);
  double tgeqrf = timer.elapsed();
  using ma::gqr;
  timer.restart();
  gqr(A, TAU, WORK);
  double tgqr = timer.elapsed();
  out << "  " << std::setw(5) << m << " " << std::setw(5) << m << " " << std::scientific << tgeqrf << " "
      << " " << tgqr << "\n";
}

template<class Allocator, class Buff>
void timeExchangeKernel(std::ostream& out, Allocator& alloc, Buff& buffer, int nbatch, int nwalk, int nocc, int nchol)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor3D_ref<T> Twabn(buffer.origin(), {2 * nbatch, nwalk * nocc, nocc * nchol});
  offset += Twabn.num_elements();
  Tensor1D_ref<T> scal(buffer.origin() + offset, boost::multi::iextensions<1u>{nbatch});
  offset += scal.num_elements();
  Tensor1D_ref<T> result(buffer.origin() + offset, boost::multi::iextensions<1u>{nwalk});
  using ma::batched_dot_wabn_wban;
  Timer timer;
  batched_dot_wabn_wban(nbatch, nwalk, nocc, nchol, scal.origin(), Twabn.origin(), to_address(result.data()), 1);
  double time = timer.elapsed();
  out << "    " << std::setw(5) << nbatch << " " << std::setw(5) << nwalk << " " << std::setw(5) << nocc << " "
      << std::setw(5) << nchol << "    " << std::scientific << time << "\n";
}

template<class Allocator, class Buff>
void timeBatchedGemm(std::ostream& out, Allocator& alloc, Buff& buffer, int nbatch, int m)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor2D_ref<T> a(buffer.origin(), {m, m});
  offset += a.num_elements();
  Tensor2D_ref<T> b(buffer.origin() + offset, {m, m});
  offset += b.num_elements();
  Tensor3D_ref<T> c(buffer.origin() + offset, {nbatch, m, m});
  //float scale = float(100.0);
  std::vector<pointer<T>> A_array;
  std::vector<pointer<T>> B_array;
  std::vector<pointer<T>> C_array;
  float alpha = 1.0;
  float beta  = 0.0;
  for (int i = 0; i < nbatch; i++)
  {
    A_array.emplace_back(a.origin());
    B_array.emplace_back(b.origin());
    C_array.emplace_back(c[i].origin());
  }
  using ma::gemmBatched;
  Timer timer;
  gemmBatched('N', 'N', m, m, m, alpha, A_array.data(), m, B_array.data(), m, beta, C_array.data(), m, nbatch);
  double tgemm = timer.elapsed();
  out << "  " << std::setw(6) << nbatch << " " << std::setw(5) << m << " " << std::scientific << tgemm << "\n";
}

template<class Allocator, class Buff>
void timeGemm(std::ostream& out, Allocator& alloc, Buff& buffer, int m, int n)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor2D_ref<T> a(buffer.origin(), {m, m});
  offset += a.num_elements();
  Tensor2D_ref<T> b(buffer.origin() + offset, {m, n});
  offset += b.num_elements();
  Tensor2D_ref<T> c(buffer.origin() + offset, {m, n});
  using ma::product;
  Timer timer;
  product(a, b, c);
  double tproduct = timer.elapsed();
  out << "  " << std::setw(6) << m << " " << std::setw(5) << n << " " << std::scientific << tproduct << "\n";
}

template<class Allocator, class Buff>
void timeBatchedMatrixInverse(std::ostream& out, Allocator& alloc, Buff& buffer, int nbatch, int m)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor3D_ref<T> a(buffer.origin(), {nbatch, m, m});
  Tensor3D_ref<T> b(buffer.origin() + a.num_elements(), {nbatch, m, m});
  Alloc<int> ialloc{};
  Tensor1D<int> IWORK(boost::multi::iextensions<1u>{nbatch * (m + 1)}, ialloc);
  std::vector<pointer<T>> A_array, B_array;
  A_array.reserve(nbatch);
  B_array.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
  {
    A_array.emplace_back(a[i].origin());
    B_array.emplace_back(b[i].origin());
  }
  using ma::getrfBatched;
  Timer timer;
  getrfBatched(m, A_array.data(), m, ma::pointer_dispatch(IWORK.origin()),
               ma::pointer_dispatch(IWORK.origin()) + nbatch * m, nbatch);
  double tgetrf = timer.elapsed();
  using ma::getriBatched;
  timer.restart();
  getriBatched(m, A_array.data(), m, ma::pointer_dispatch(IWORK.origin()), B_array.data(), m,
               ma::pointer_dispatch(IWORK.origin()) + nbatch * m, nbatch);
  double tgetri = timer.elapsed();
  out << "  " << std::setw(6) << nbatch << " " << std::setw(5) << m << " " << std::scientific << tgetrf << " " << tgetri
      << "\n";
}

template<class Allocator, class Buff>
void timeMatrixInverse(std::ostream& out, Allocator& alloc, Buff& buffer, int m)
{
  using T    = typename Allocator::value_type;
  int offset = 0;
  Tensor2D_ref<T> a(buffer.origin(), {m, m});
  Tensor1D_ref<T> WORK(buffer.origin() + a.num_elements(), boost::multi::iextensions<1u>{m * m});
  Alloc<int> ialloc{};
  Tensor1D<int> IWORK(boost::multi::iextensions<1u>{m + 1}, ialloc);
  using ma::getrf;
  Timer timer;
  getrf(a, IWORK, WORK);
  double tgetrf = timer.elapsed();
  using ma::getri;
  timer.restart();
  getri(a, IWORK, WORK);
  double tgetri = timer.elapsed();
  out << "  " << std::setw(6) << m << " " << std::scientific << tgetrf << " " << tgetri << "\n";
}

int main(int argc, char* argv[])
{
  boost::mpi3::environment env(argc, argv);
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  arch::INIT(node);
#endif
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
/*
  {
    std::ofstream out;
    out.open("time_batched_zqr.dat");
    std::cout << " - Batched zQR (nbatch, MxN)" << std::endl;
    out << " nbatch       M     N       ttrans       tgeqrf      tgetdet         tgqr\n";
    std::vector<int> batches  = {1, 5, 10, 20};
    std::vector<int> num_rows = {200, 400, 800};
    int max_batch             = batches[batches.size() - 1];
    int max_rows              = num_rows[num_rows.size() - 1];
    int size                  = (2 * max_batch * max_rows * (max_rows / 2.0) + 3 * max_batch * max_rows);
    Alloc<std::complex<double>> alloc{};
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, 1.0, alloc);
    for (auto nb : batches)
    {
      for (auto m : num_rows)
      {
        timeBatchedQR(out, alloc, buffer, nb, m, m / 2);
      }
    }
  }
*/
#endif
  {
    std::ofstream out;
    out.open("time_zqr.dat");
    std::cout << " - zQR (MxM)" << std::endl;
    out << "      M     M      tzgeqrf        tzungqr\n";
    int size = 3 * 1000 * 1000;
    Alloc<std::complex<double>> alloc{};
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, 1.0, alloc);
    std::vector<int> dims = {100, 200, 500, 800, 1000};
    for (auto d : dims)
    {
      timeQR(out, alloc, buffer, d);
    }
  }
  {
    std::ofstream out;
    out.open("time_sgemm.dat");
    std::cout << " - sgemm (MxM)" << std::endl;
    out << "       M     M       tsgemm\n";
    int size = 3 * 8000 * 8000;
    Alloc<float> alloc{};
    Tensor1D<float> buffer(iextensions<1u>{size}, 1.0, alloc);
    std::vector<int> dims = {200, 500, 800, 1000, 2000, 3000, 4000, 8000};
    for (auto d : dims)
    {
      timeGemm(out, alloc, buffer, d, d);
    }
  }
  {
    std::ofstream out;
    out.open("time_batched_sgemm.dat");
    std::cout << " - batched sgemm (nbatch, MxM)" << std::endl;
    out << "  nbatch    M        tsgemm\n";
    Alloc<float> alloc{};
    std::vector<int> num_rows = {100, 200, 300, 400, 500, 600};
    std::vector<int> batches  = {128, 256, 512, 1024};
    int max_batch             = batches[batches.size() - 1];
    int max_rows              = num_rows[num_rows.size() - 1];
    int size                  = 3 * max_batch * max_rows * max_rows;
    Tensor1D<float> buffer(iextensions<1u>{size}, 1.0, alloc);
    for (auto nb : batches)
    {
      for (auto m : num_rows)
      {
        timeBatchedGemm(out, alloc, buffer, nb, m);
      }
    }
  }
  {
    std::ofstream out;
    out.open("time_exchange_kernel.dat");
    std::cout << " - exchange kernel (E[w] = sum_{abn} Twabn Twabn)" << std::endl;
    out << "   nbatch nwalk  nocc nchol tExchangeKernel\n";
    Alloc<std::complex<double>> alloc{};
    int nwalk                = 5;
    int nocc                 = 20;
    int nchol                = 270;
    std::vector<int> batches = {100, 200, 400, 800};
    int nbatch_max           = batches[batches.size() - 1];
    int size                 = 2 * nbatch_max * nwalk * nocc * nocc * nchol + nbatch_max + nwalk;
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, 1.0, alloc);
    for (auto b : batches)
    {
      timeExchangeKernel(out, alloc, buffer, b, nwalk, nocc, nchol);
    }
  }
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  {
    std::ofstream out;
    out.open("time_batched_matrix_inverse.dat");
    std::cout << " - batched matrix inverse (nbatch, MxM)" << std::endl;
    out << "  nbatch     M       tgetrf       tgetri\n";
    Alloc<std::complex<double>> alloc{};
    std::vector<int> batches  = {1, 5, 10, 20};
    std::vector<int> num_rows = {100, 110, 120, 200, 210, 300, 400, 500, 600, 700};
    int max_batch             = batches[batches.size() - 1];
    int max_rows              = num_rows[num_rows.size() - 1];
    int size                  = 2 * max_batch * max_rows * max_rows;
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, alloc);
    {
      std::vector<std::complex<double>> tmp(size);
      fillRandomMatrix(tmp);
      using std::copy_n;
      copy_n(tmp.data(), tmp.size(), buffer.origin());
    }
    for (auto b : batches)
    {
      for (auto m : num_rows)
      {
        timeBatchedMatrixInverse(out, alloc, buffer, b, m);
      }
    }
  }
#endif
  {
    std::ofstream out;
    out.open("time_matrix_inverse.dat");
    std::cout << " - matrix inverse (nbatch, MxM)" << std::endl;
    out << "       M       tgetrf       tgetri\n";
    Alloc<std::complex<double>> alloc{};
    std::vector<int> num_rows = {100, 110, 120, 200, 210, 300, 400, 500, 600, 700, 800, 1000, 2000, 4000};
    int max_rows              = num_rows[num_rows.size() - 1];
    int size                  = 2 * max_rows * max_rows;
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, alloc);
    {
      std::vector<std::complex<double>> tmp(size);
      fillRandomMatrix(tmp);
      using std::copy_n;
      copy_n(tmp.data(), tmp.size(), buffer.origin());
    }
    for (auto m : num_rows)
    {
      timeMatrixInverse(out, alloc, buffer, m);
    }
  }
}
