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

#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Memory/buffer_allocators.h"
#include "AFQMC/Memory/arch.hpp"
#include "AFQMC/Utilities/myTimer.h"

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

template<class Allocator, class Buff>
void timeBatchedQR(Allocator& alloc, Buff& buffer, int nbatch, int m, int n)
{
  using T = typename Allocator::value_type;
  int offset = 0;
  Tensor3D_ref<T> A(buffer.origin(), {nbatch, m, n});
  offset += A.num_elements();
  Tensor3D_ref<T> AT(buffer.origin()+offset, {nbatch, n, m});
  offset += AT.num_elements();
  Tensor2D_ref<T> T_(buffer.origin()+offset, {nbatch, m});
  offset += T_.num_elements();
  Tensor2D_ref<T> scl(buffer.origin()+offset, {nbatch, m});
  offset += T_.num_elements();
  int sz = ma::gqr_optimal_workspace_size(AT[0]);
  //std::cout << buffer.num_elements() << " " << 2*nbatch*m*n + 2*nbatch*m + nbatch*sz << " " << offset << std::endl; 
  Tensor1D_ref<T> WORK(buffer.origin()+offset, boost::multi::iextensions<1u>{nbatch * sz});
  Alloc<int> ialloc{};
  std::vector<pointer<T>> Aarray;
  Tensor1D<int> IWORK(boost::multi::iextensions<1u>{nbatch * (m + 1)}, ialloc);
  using std::copy_n;
  for (int i = 0; i < nbatch; i++) {
    Aarray.emplace_back(A[i].origin());
  }
  // Actual profile.
  myTimer timer;
  std::string timer_id = std::string("transpose");
  timer.add(timer_id);
  timer.start(timer_id);
  for (int i = 0; i < nbatch; i++)
    ma::transpose(A[i], AT[i]);
  auto ttrans = timer.elapsed(timer_id);
  timer_id = std::string("geqrf");
  timer.add(timer_id);
  timer.start(timer_id);
  geqrfStrided(m, n, AT.origin(), m, m*n, T_.origin(), m, IWORK.origin(), nbatch);
  auto tgeqrf = timer.elapsed(timer_id);
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
  std::cout << "  " << std::setw(5) << nbatch << "   " << std::setw(5) << m << " " << std::setw(5)
            << n << " " << std::scientific << ttrans << " " << tgeqrf << " " << tdet << " " <<
            tgqr << std::endl;
}

template<class Allocator, class Buff>
void timeBatchedGemm(Allocator& alloc, Buff& buffer, int nbatch, int m)
{
  using T = typename Allocator::value_type;
  myTimer timer;
  std::string timer_id = std::to_string(nbatch);
  timer.add(timer_id);
  int offset = 0;
  Tensor2D_ref<T> a(buffer.origin(), {m, m});
  offset += a.num_elements();
  Tensor2D_ref<T> b(buffer.origin()+offset, {m, m});
  offset += b.num_elements();
  Tensor3D_ref<T> c(buffer.origin()+offset, {nbatch, m, m});
  //float scale = float(100.0);
  //fill_matrix(a, scale);
  //fill_matrix(b, scale);
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
  timer.start(timer_id);
  gemmBatched('N', 'N', m, m, m, alpha, A_array.data(), m, B_array.data(), m, beta, C_array.data(), m, nbatch);
  timer.stop(timer_id);
  std::cout << "  " << std::setw(6) << nbatch << " " << std::setw(5) << m
            << " " << std::scientific << timer.elapsed(timer_id) << std::endl;
}

template<class Allocator, class Buff>
void timeGemm(Allocator& alloc, Buff& buffer, int m, int n)
{
  using T = typename Allocator::value_type;
  myTimer timer;
  std::string timer_id = std::to_string(m);
  timer.add(timer_id);
  int offset = 0;
  Tensor2D_ref<T> a(buffer.origin(), {m, m});
  offset += a.num_elements();
  Tensor2D_ref<T> b(buffer.origin()+offset, {m, n});
  offset += b.num_elements();
  Tensor2D_ref<T> c(buffer.origin()+offset, {m, n});
  using ma::product;
  timer.start(timer_id);
  product(a, b, c);
  timer.stop(timer_id);
  std::cout << "  " << std::setw(6) << m << " " << std::setw(5) << n
            << " " << std::scientific << timer.elapsed(timer_id) << std::endl;
}

void parse_args(int argc, char* argv[], int &num_walker, int &num_basis, int &num_elec)
{
  char* pend;
  if (argc < 4) {
    std::cout << "usage : afqmc_performance --nwalker nw --nbasis nb --nelec ne" << std::endl;
    exit(1);
  }
  for (int i = 1; i < argc; i++) {
    if (i <= argc)
    {
      if (!std::strcmp(argv[i], "--nwalker"))
      {
        num_walker = strtol(argv[i + 1], &pend, 10);
        i++;
      }
      if (!std::strcmp(argv[i], "--nbasis"))
      {
        num_basis = strtol(argv[i + 1], &pend, 10);
        i++;
      }
      if (!std::strcmp(argv[i], "--nelec"))
      {
        num_elec = strtol(argv[i + 1], &pend, 10);
        i++;
      }
    }
  }
  if (num_elec <= 0) {
    std::cout << "usage : afqmc_performance --nwalk nw --nbasis nb --nelec ne" << std::endl;
    if (num_elec < 0) {
      std::cout << "number of electrons (--nelec) must be greater than zero" << std::endl;
    }
    exit(1);
  }
  if (num_basis <= 0) {
    std::cout << "usage : afqmc_performance --nwalk nw --nbasis nb --nelec ne" << std::endl;
    if (num_basis < 0) {
      std::cout << "number of basis functions (--nbasis) must be greater than zero" << std::endl;
    }
    exit(1);
  }
  if (num_walker <= 0) {
    std::cout << "usage : afqmc_performance --nwalk nw --nbasis nb --nelec ne" << std::endl;
    if (num_walker < 0) {
      std::cout << "number of walkers (--nwalk) must be greater than zero" << std::endl;
    }
    exit(1);
  }
  std::cout << " - nw : " << num_walker << std::endl;
  std::cout << " - nb : " << num_basis << std::endl;
  std::cout << " - ne : " << num_elec << std::endl;
  std::cout << std::endl;
}

int main(int argc, char* argv[])
{
  boost::mpi3::environment env(argc, argv);
  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared(world.rank());
  arch::INIT(node);
  int num_walker = -1;
  int num_basis = -1;
  int num_elec = -1;
  parse_args(argc, argv, num_walker, num_basis, num_elec);
  std::vector<int> num_kpoints = {8, 12, 18, 27, 36, 48};
  int nk_max = num_kpoints[num_kpoints.size()-1];
  int m_max = nk_max * num_basis;
  int n_max = nk_max * num_elec;
  {
    std::cout << " - Batched sgemm (nk^3, mxm)" << std::endl;
    std::cout << std::endl;
    std::cout << "    nk^3     m         time" << std::endl;
    int nb_max = nk_max*nk_max;
    int nb = 10*num_basis;
    int size = nb_max*nb*nb+ 2*nb*nb;
    Alloc<float> alloc{};
    Tensor1D<float> buffer(iextensions<1u>{size}, 1.0, alloc);
    for (auto nk : num_kpoints) {
      int num_batch = nk * nk;
      timeBatchedGemm(alloc, buffer, num_batch, nb);
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  {
    std::cout << " - Batched zQR (nwalk, mxn)" << std::endl;
    std::cout << std::endl;
    std::cout << "  nwalk  nbasis nelec       ttrans       tgeqrf      tgetdet         tgqr" << std::endl;
    //Alloc<std::complex<double>> alloc{};
    // Allocate memory to speed things up.
    int size = (2*num_walker*m_max*n_max + 3*num_walker*m_max);
    Alloc<std::complex<double>> alloc{};
    Tensor1D<std::complex<double>> buffer(iextensions<1u>{size}, 1.0, alloc);
    //auto alloc(device_buffer_generator->template get_allocator<std::complex<double>>());
    for (auto nk : num_kpoints) {
      int m = num_basis * nk;
      int n = num_elec * nk;
      timeBatchedQR(alloc, buffer, num_walker, m, n);
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  {
    std::cout << " - sgemm (MxM)" << std::endl;
    std::cout << std::endl;
    std::cout << "      M    M" << std::endl;
    int size = (3*m_max*m_max);
    Alloc<float> alloc{};
    Tensor1D<float> buffer(iextensions<1u>{size}, 1.0, alloc);
    //auto alloc(device_buffer_generator->template get_allocator<std::complex<double>>());
    for (auto nk : num_kpoints) {
      int m = num_basis * nk;
      timeGemm(alloc, buffer, m, m);
    }
    //std::cout << std::endl;
    //std::cout << std::endl;
  }
}
