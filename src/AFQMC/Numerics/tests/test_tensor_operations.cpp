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
#define APP_ABORT(x) {std::cout << x; exit(0);}

#include <vector>


#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/tensor_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"


using boost::multi::array;
using boost::multi::array_ref;
using std::copy_n;
using boost::multi::iextensions;

namespace qmcplusplus
{

using namespace afqmc;

#ifdef ENABLE_CUDA
template<typename T> using Alloc = device::device_allocator<T>;
#else
template<typename T> using Alloc = std::allocator<T>;
#endif
template<typename T> using pointer = typename Alloc<T>::pointer;


template<typename T> using Tensor1D = array<T,1,Alloc<T>>;
template<typename T> using Tensor2D = array<T,2,Alloc<T>>;
template<typename T> using Tensor3D = array<T,3,Alloc<T>>;
template<typename T> using Tensor4D = array<T,4,Alloc<T>>;
template<typename T> using Tensor5D = array<T,5,Alloc<T>>;

template<typename T> using Tensor2D_ref = array_ref<T,2,pointer<T>>;
template<typename T> using Tensor1D_ref = array_ref<T,1,pointer<T>>;

template <typename T>
void create_data(std::vector<T>& buffer, T scale)
{
  T count = T(0);
  for (int i=0; i<buffer.size(); i++) {
    buffer[i] = count;
    count += T(1) / scale;
  }
}

TEST_CASE("KaKjw_to_KKwaj", "[Numerics][tensor_operations]")
{

  Alloc<double> alloc{};
  int nwalk = 5;
  int occ_k = 5;
  int nmo_k = 10;
  int nk = 3;
  std::vector<double> buffer(nk*occ_k*nk*nmo_k*nwalk);
  create_data(buffer, 1.0);
  Tensor5D<double> KaKjw({nk, occ_k, nk, nmo_k, nwalk}, alloc);
  Tensor5D<double> KKwaj({nk, nk, nwalk, occ_k, nmo_k}, 0.0, alloc);
  Tensor1D<int> nmo_pk = {nmo_k, nmo_k, nmo_k};
  Tensor1D<int> nel_pk = {occ_k, occ_k, occ_k};
  Tensor1D<int> nmo_pk0 = {0, nmo_k, 2*nmo_k};
  Tensor1D<int> nel_pk0 = {0, occ_k, 2*occ_k};
  copy_n(buffer.data(), buffer.size(), KaKjw.origin());
  // KaKjw = numpy.arange(nk*na*nk*nj*nw).reshape((nk,na,nk,nj,nw))
  // KKwaj = numpy.transpose(KaKjw, (0,2,1,3,4))
  using ma::KaKjw_to_KKwaj;
  KaKjw_to_KKwaj(nwalk, nk, nmo_k, nk*nmo_k, occ_k,
                 nmo_pk.origin(), nmo_pk0.origin(),
                 nel_pk.origin(), nel_pk0.origin(),
                 KaKjw.origin(), KKwaj.origin());
  // Reference values from python
  REQUIRE(KKwaj[1][2][3][4][4] == Approx(1473.0));
  REQUIRE(KKwaj[2][1][3][4][4] == Approx(2173.0));
  REQUIRE(KKwaj[2][1][4][3][4] == Approx(2024.0));

}

TEST_CASE("KaKjw_to_QKwaj", "[Numerics][tensor_operations]")
{

  Alloc<double> alloc{};
  int nwalk = 5;
  int occ_k = 5;
  int nmo_k = 10;
  int nk = 3;
  std::vector<double> buffer(nk*occ_k*nk*nmo_k*nwalk);
  create_data(buffer, 1.0);
  Tensor3D<double> KaKjw({nk*occ_k, nk*nmo_k, nwalk}, alloc);
  Tensor5D<double> QKajw({nk, nk, occ_k, nmo_k, nwalk}, 0.0, alloc);
  Tensor1D<int> nmo_pk = {nmo_k, nmo_k, nmo_k};
  Tensor1D<int> nel_pk = {occ_k, occ_k, occ_k};
  Tensor1D<int> nmo_pk0 = {0, nmo_k, 2*nmo_k};
  Tensor1D<int> nel_pk0 = {0, occ_k, 2*occ_k};
  Tensor1D<int> qk_to_k2 = {1, 2, 0};
  copy_n(buffer.data(), buffer.size(), KaKjw.origin());
  using ma::KaKjw_to_QKajw;
  KaKjw_to_QKajw(nwalk, nk, nmo_k, nk*nmo_k, occ_k,
                 nmo_pk.origin(), nmo_pk0.origin(),
                 nel_pk.origin(), nel_pk0.origin(),
                 qk_to_k2.origin(),
                 KaKjw.origin(), QKajw.origin());
  // Just captured from output.
  REQUIRE(QKajw[1][2][3][4][4] == Approx(1974.0));
  REQUIRE(QKajw[2][1][3][4][4] == Approx(1224.0));
  REQUIRE(QKajw[2][1][4][3][4] == Approx(1369.0));

}

TEST_CASE("term_by_term_matrix_vector", "[Numerics][tensor_operations]")
{

  Alloc<ComplexType> alloc{};
  int nrow = 3;
  int ncol = 3;
  Tensor2D<ComplexType> A({nrow,ncol}, alloc);
  Tensor1D<ComplexType> x(iextensions<1u>{ncol}, alloc);
  std::vector<ComplexType> buffer(nrow*ncol);
  create_data(buffer, ComplexType(1.0));
  copy_n(buffer.data(), buffer.size(), A.origin());
  copy_n(buffer.data(), x.size(), x.origin());
  using ma::term_by_term_matrix_vector;
  Tensor2D<ComplexType> ref = {{0,1,2},{4,5,6},{8,9,10}};
  term_by_term_matrix_vector(ma::TOp_PLUS, 0, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  verify_approx(ref,A);
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_PLUS, 1, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,2,4},{3,5,7},{6,8,10}};
  verify_approx(ref,A);
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_MINUS, 0, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,1,2},{2,3,4},{4,5,6}};
  verify_approx(ref,A);
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_MINUS, 1, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,0,0},{3,3,3},{6,6,6}};
  verify_approx(ref,A);
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_MUL, 0, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,0,0},{3,4,5},{12,14,16}};
  verify_approx(ref,A);
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_MUL, 1, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,1,4},{0,4,10},{0,7,16}};
  verify_approx(ref,A);
  // Avoid dividing by zero
  copy_n(buffer.data()+1, x.size(), x.origin());
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_DIV, 0, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0,1,2},{1.5,2.0,2.5},{2,2.333333333333,2.666666666667}};
  verify_approx(ref,A);
  copy_n(buffer.data()+1, x.size(), x.origin());
  copy_n(buffer.data(), buffer.size(), A.origin());
  term_by_term_matrix_vector(ma::TOp_DIV, 1, nrow, ncol,
                             A.origin(), ncol,
                             x.origin(), 1);
  ref = {{0.0,0.5,0.666666666667},
         {3.0,2.0,1.666666666667},
         {6.0,3.5,2.666666666667}};
  verify_approx(ref,A);

}

//TEST_CASE("vKKwij_tovwKiKj", "[Numerics][tensor_operations]")
//{

  //Alloc<double> alloc{};
  //int nwalk = 5;
  //int occ_k = 5;
  //int nmo_k = 10;
  //int nk = 3;
  //std::vector<double> buffer(nk*occ_k*nk*nmo_k*nwalk);
  //create_data(buffer, 1.0);
  //Tensor5D<double> KaKjw({nk, occ_k, nk, nmo_k, nwalk}, alloc);
  //Tensor5D<double> QKajw({nk, nk, occ_k, nmo_k, nwalk}, 0.0, alloc);
  //Tensor1D<int> nmo_pk = {nmo_k, nmo_k, nmo_k};
  //Tensor1D<int> nel_pk = {occ_k, occ_k, occ_k};
  //Tensor1D<int> nmo_pk0 = {0, nmo_k, 2*nmo_k};
  //Tensor1D<int> nel_pk0 = {0, occ_k, 2*occ_k};
  //Tensor1D<int> qk_to_k2 = {1, 2, 0};
  //copy_n(buffer.data(), buffer.size(), KaKjw.origin());
  //using ma::KaKjw_to_QKajw;
  //KaKjw_to_QKajw(nwalk, nk, nmo_k, nk*nmo_k, occ_k,
                 //nmo_pk.origin(), nmo_pk0.origin(),
                 //nel_pk.origin(), nel_pk0.origin(),
                 //qk_to_k2.origin(),
                 //KaKjw.origin(), QKajw.origin());
  //// Just captured from output.
  //REQUIRE(QKajw[1][2][3][4][4] == Approx(1974.0));
  //REQUIRE(QKajw[2][1][3][4][4] == Approx(1224.0));
  //REQUIRE(QKajw[2][1][4][3][4] == Approx(1369.0));

//}

}
