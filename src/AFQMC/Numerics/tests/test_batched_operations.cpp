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
template<typename T>
using Tensor4D = array<T, 4, Alloc<T>>;
template<typename T>
using Tensor5D = array<T, 5, Alloc<T>>;

template<typename T>
using Tensor2D_ref = array_ref<T, 2, pointer<T>>;
template<typename T>
using Tensor1D_ref = array_ref<T, 1, pointer<T>>;

template<typename T>
void create_data(std::vector<T>& buffer, T scale)
{
  T count = T(0);
  for (int i = 0; i < buffer.size(); i++)
  {
    buffer[i] = count;
    count += T(1) / scale;
  }
}

// Test data for assertions was created from captured output.
TEST_CASE("Tab_to_Kl", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  // Tab_to_Kl tests, normally used for Coulomb energy evaluation via
  // E_c[w] = \sum_{na} T[w,aa,n] T[w,aa,n]
  //        = \sum_{n} K[w,n] K[w,n]
  int nwalk = 3;
  int nel   = 3;
  int nchol = 7;
  std::vector<ComplexType> buffer(nwalk * nel * nel * nchol);
  create_data(buffer, ComplexType(1.0));
  Tensor4D<ComplexType> Twban({nwalk, nel, nel, nchol}, alloc);
  copy_n(buffer.data(), buffer.size(), Twban.origin());
  Tensor2D<ComplexType> Kl({nwalk, nchol}, 0.0, alloc);
  using ma::Tab_to_Kl;
  Tab_to_Kl(nwalk, nel, nchol, Twban.origin(), Kl.origin());
  array_ref<ComplexType, 1, pointer<ComplexType>> Kl_(Kl.origin(), iextensions<1u>{nwalk * nchol});
  array<ComplexType, 1, Alloc<ComplexType>> ref = {84.0,  87.0,  90.0,  93.0,  96.0,  99.0,  102.0,
                                                   273.0, 276.0, 279.0, 282.0, 285.0, 288.0, 291.0,
                                                   462.0, 465.0, 468.0, 471.0, 474.0, 477.0, 480.0};
  verify_approx(Kl_, ref);
}

TEST_CASE("batched_Tab_to_Klr", "[batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nwalk              = 3;
  int nel                = 3;
  int nchol              = 11;
  int nbatch             = 4;
  int nchol_max          = 3;
  int ncholQ             = 2;
  int ncholQ0            = 2;
  std::vector<int> kdiag = {0, 1, 2, 3};
  Tensor1D<int> dev_kdiag(iextensions<1u>{nbatch}, alloc);
  copy_n(kdiag.data(), kdiag.size(), dev_kdiag.origin());
  std::vector<ComplexType> buffer(2 * nbatch * nwalk * nel * nel * nchol_max);
  create_data(buffer, ComplexType(1.0));
  Tensor5D<ComplexType> Tab({2 * nbatch, nwalk, nel, nel, nchol_max}, alloc);
  copy_n(buffer.data(), buffer.size(), Tab.origin());
  Tensor2D<ComplexType> Kl({nwalk, 2 * nchol}, 0.0, alloc);
  Tensor2D<ComplexType> Kr({nwalk, 2 * nchol}, 0.0, alloc);
  using ma::batched_Tab_to_Klr;
  batched_Tab_to_Klr(nbatch, nwalk, nel, nchol_max, nchol, ncholQ, ncholQ0, dev_kdiag.origin(), Tab.origin(),
                     Kl.origin(), Kr.origin());
  copy_n(Kr.origin(), Kr.num_elements(), buffer.data());
  //std::cout << std::setprecision(16) << Kl[2][3] << " " << Kl[1][4] << " " << Kr[1][3] << " " << Kr[0][1] << std::endl;
  CHECK(real(buffer[2 * nchol + 3]) == Approx(2262));
}

TEST_CASE("Tanb_to_Kl", "[batched_operations]")
{
  Alloc<ComplexType> alloc{};
  // Tab_to_Kl tests, normally used for Coulomb energy evaluation via
  // E_c[w] = \sum_{na} T[w,aa,n] T[w,aa,n]
  //        = \sum_{n} K[w,n] K[w,n]
  int nwalk = 3;
  int nel   = 3;
  int nchol = 7;
  std::vector<ComplexType> buffer(nwalk * nel * nel * nchol);
  create_data(buffer, ComplexType(1.0));
  Tensor4D<ComplexType> Twanb({nwalk, nel, nchol, nel}, alloc);
  copy_n(buffer.data(), buffer.size(), Twanb.origin());
  Tensor2D<ComplexType> Kl({nwalk, nchol}, 0.0, alloc);
  using ma::Tanb_to_Kl;
  Tanb_to_Kl(nwalk, nel, nchol, nchol, Twanb.origin(), Kl.origin());
  array_ref<ComplexType, 1, pointer<ComplexType>> Kl_(Kl.origin(), iextensions<1u>{nwalk * nchol});
  //std::cout << "{";
  //for (auto i : Kl_)
  //std::cout << "ComplexType(" << real(i) << ")," << std::endl;
  //std::cout << "};" << std::endl;;
  array<ComplexType, 1, Alloc<ComplexType>> ref = {ComplexType(66),  ComplexType(75),  ComplexType(84),
                                                   ComplexType(93),  ComplexType(102), ComplexType(111),
                                                   ComplexType(120), ComplexType(255), ComplexType(264),
                                                   ComplexType(273), ComplexType(282), ComplexType(291),
                                                   ComplexType(300), ComplexType(309), ComplexType(444),
                                                   ComplexType(453), ComplexType(462), ComplexType(471),
                                                   ComplexType(480), ComplexType(489), ComplexType(498)};
  verify_approx(Kl_, ref);
}

TEST_CASE("batched_dot_wabn_wban", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nbatch = 4;
  int nwalk  = 3;
  int nocc   = 7;
  int nchol  = 13;
  // 2*nbatch as routine expects left and right tensors to contract.
  std::vector<ComplexType> buffer(2 * nbatch * nwalk * nocc * nocc * nchol);
  create_data(buffer, ComplexType(1e4));
  Tensor5D<ComplexType> Twabn({2 * nbatch, nwalk, nocc, nocc, nchol}, alloc);
  Tensor1D<ComplexType> scal({nbatch}, 1.0, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  std::vector<pointer<ComplexType>> Aarray;
  array<ComplexType, 1, Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  array<ComplexType, 1, Alloc<ComplexType>> ref = {ComplexType(1693.073254599684), ComplexType(1930.853888599637),
                                                   ComplexType(2189.312510839587)};
  using ma::batched_dot_wabn_wban;
  batched_dot_wabn_wban(nbatch, nwalk, nocc, nchol, scal.origin(), Twabn.origin(), to_address(out.data()), 1);
  //std::cout << std::setprecision(16) << "this: " <<  out[0] << " " << out[1] << " " << out[2] << std::endl;
  verify_approx(ref, out);
}

// Not used.
TEST_CASE("batched_dot_wanb_wbna", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nbatch = 4;
  int nwalk  = 3;
  int nocc   = 7;
  int nchol  = 13;
  // 2*nbatch as routine expects left and right tensors to contract.
  std::vector<ComplexType> buffer(2 * nbatch * nwalk * nocc * nocc * nchol);
  create_data(buffer, ComplexType(1e4));
  Tensor5D<ComplexType> Twabn({2 * nbatch, nwalk, nocc, nchol, nocc}, alloc);
  Tensor1D<ComplexType> scal({nbatch}, 1.0, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  std::vector<pointer<ComplexType>> Aarray;
  array<ComplexType, 1, Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  array<ComplexType, 1, Alloc<ComplexType>> ref = {ComplexType(1692.867783879684), ComplexType(1930.648417879638),
                                                   ComplexType(2189.107040119586)};
  using ma::batched_dot_wanb_wbna;
  batched_dot_wanb_wbna(nbatch, nwalk, nocc, nchol, scal.origin(), Twabn.origin(), to_address(out.data()), 1);
  //std::cout << std::setprecision(16) << "this: " <<  out[0] << " " << out[1] << " " << out[2] << std::endl;
  verify_approx(ref, out);
}


TEST_CASE("dot_wabn", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nwalk = 3;
  int nocc  = 7;
  int nchol = 13;
  std::vector<ComplexType> buffer(nwalk * nocc * nocc * nchol);
  create_data(buffer, ComplexType(100));
  Tensor4D<ComplexType> Twabn({nwalk, nocc, nocc, nchol}, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  array<ComplexType, 1, Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  using ma::dot_wabn;
  dot_wabn(nwalk, nocc, nchol, ComplexType(1.0), Twabn.origin(), to_address(out.origin()), 1);
  array<ComplexType, 1, Alloc<ComplexType>> ref = {ComplexType(7045.35), ComplexType(58699.7), ComplexType(162049.0)};
  verify_approx(out, ref);
}

TEST_CASE("dot_wanb", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nwalk = 3;
  int nocc  = 7;
  int nchol = 13;
  std::vector<ComplexType> buffer(nwalk * nocc * nocc * nchol);
  // T = numpy.arange(nw*nocc*nocc*nchol).reshape(nw,nocc,nchol,nocc) / 100.0
  create_data(buffer, ComplexType(100));
  Tensor4D<ComplexType> Twanb({nwalk, nocc, nchol, nocc}, alloc);
  copy_n(buffer.data(), buffer.size(), Twanb.origin());
  array<ComplexType, 1, Alloc<ComplexType>> out(iextensions<1u>{nwalk}, 0.0, alloc);
  using ma::dot_wanb;
  // out = numpy.einsum('wanb,wbna->w', Twanb, Twanb)
  dot_wanb(nwalk, nocc, nchol, ComplexType(1.0), Twanb.origin(), to_address(out.origin()), 1);
  array<ComplexType, 1, Alloc<ComplexType>> ref = {ComplexType(6531.67), ComplexType(58186.1), ComplexType(161535)};
  verify_approx(out, ref);
}

TEST_CASE("Auwn_Bun_Cuw", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nu = 3;
  int nw = 2;
  int nn = 4;
  Tensor3D<ComplexType> A({nu, nw, nn}, 1.0, alloc);
  Tensor2D<ComplexType> B({nu, nn}, 2.0, alloc);
  Tensor2D<ComplexType> C({nu, nw}, 0.0, alloc);
  ComplexType alpha = 0.5;
  // C = alpha * numpy.einsum('uwn,un->uw', A, B)
  using ma::Auwn_Bun_Cuw;
  Auwn_Bun_Cuw(nu, nw, nn, alpha, A.origin(), B.origin(), C.origin());
  Tensor2D<ComplexType> ref({nu, nw}, 4.0, alloc);
  verify_approx(C, ref);
}

TEST_CASE("Awiu_Biu_Cuw", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nu = 3;
  int nw = 2;
  int nn = 4;
  Tensor3D<ComplexType> A({nw, nn, nu}, 1.0, alloc);
  Tensor2D<ComplexType> B({nn, nu}, 2, alloc);
  B[0][1] = 0.0;
  Tensor2D<ComplexType> C({nu, nw}, 0.0, alloc);
  ComplexType alpha = 0.5;
  // C = alpha * numpy.einsum('wnu,nu->uw', A, B)
  using ma::Awiu_Biu_Cuw;
  Awiu_Biu_Cuw(nu, nw, nn, alpha, A.origin(), B.origin(), std::get<1>(B.sizes()), C.origin(), std::get<1>(C.sizes()));
  Tensor2D<ComplexType> ref({nu, nw}, 4.0, alloc);
  ref[1][0] = 3.0;
  ref[1][1] = 3.0;
  verify_approx(C, ref);
}

TEST_CASE("Aijk_Bkj_Cik", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int ni = 3;
  int nj = 2;
  int nk = 4;
  Tensor3D<ComplexType> A({ni, nj, nk}, 1.0, alloc);
  Tensor2D<ComplexType> B({nk, nj}, 2, alloc);
  B[0][1] = 0.0;
  Tensor2D<ComplexType> C({ni, nk}, 0.0, alloc);
  // C = alpha * numpy.einsum('wnu,nu->uw', A, B)
  using ma::Aijk_Bkj_Cik;
  Aijk_Bkj_Cik(ni, nj, nk, A.origin(), std::get<1>(A.sizes()), A.stride(0), B.origin(), B.stride(0), C.origin(), C.stride(0));
  Tensor2D<ComplexType> ref({ni, nk}, 4.0, alloc);
  ref[0][0] = 2.0;
  ref[1][0] = 2.0;
  ref[2][0] = 2.0;
  verify_approx(C, ref);
}

TEST_CASE("viwj_vwij", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int ni = 3;
  int nj = 2;
  int nw = 4;
  Tensor3D<ComplexType> A({nw, ni, nj}, 0.0, alloc);
  Tensor3D<ComplexType> B({ni, nw, nj}, 0.0, alloc);
  std::vector<ComplexType> buffer(ni * nj * nw);
  create_data(buffer, ComplexType(1.0));
  copy_n(buffer.data(), buffer.size(), B.origin());
  using ma::viwj_vwij;
  //viwj_vwij(nw, ni, 0, nj, B.data(), A.data());
  //std::cout << A[0][1][1] << " " << A[1][2][1] << std::endl;
}

TEST_CASE("element_wise_Aij_Bjk_Ckij", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int ni = 3;
  int nj = 2;
  int nk = 2;
  using ma::element_wise_Aij_Bjk_Ckij;
  {
    Tensor2D<ComplexType> A({ni, nj}, 3.0, alloc);
    Tensor2D<ComplexType> B({nj, nk}, 2.0, alloc);
    Tensor3D<ComplexType> C({nk, ni, nj}, 0.0, alloc);
    element_wise_Aij_Bjk_Ckij('N', ni, nj, nk, A.origin(), A.stride(0), B.origin(), B.stride(0), C.origin(), std::get<1>(C.sizes()),
                              std::get<2>(C.sizes()));
    Tensor3D<ComplexType> ref({nk, ni, nj}, 6.0, alloc);
    verify_approx(C, ref);
  }
  {
    Tensor2D<ComplexType> A({ni, nj}, ComplexType(0.0, -3.0), alloc);
    Tensor2D<ComplexType> B({nj, nk}, ComplexType(1.0, 2.0), alloc);
    Tensor3D<ComplexType> C({nk, ni, nj}, 0.0, alloc);
    element_wise_Aij_Bjk_Ckij('C', ni, nj, nk, A.origin(), A.stride(0), B.origin(), B.stride(0), C.origin(), std::get<1>(C.sizes()),
                              std::get<2>(C.sizes()));
    Tensor3D<ComplexType> ref({nk, ni, nj}, ComplexType(-6.0, 3.0), alloc);
    verify_approx(C, ref);
  }
}

template<typename T1, typename T2>
void test_Aij_Bjk_Ckji()
{
  Alloc<T1> alloc_a{};
  Alloc<T2> alloc_b{};
  int ni = 3;
  int nj = 2;
  int nk = 2;
  using ma::element_wise_Aij_Bjk_Ckji;
  Tensor2D<T1> A({ni, nj}, -3.0, alloc_a);
  Tensor2D<T2> B({nj, nk}, T2(1.0, 2.0), alloc_b);
  Tensor3D<T2> C({nk, nj, ni}, 0.0, alloc_b);
  element_wise_Aij_Bjk_Ckji(ni, nj, nk, A.origin(), A.stride(0), B.origin(), B.stride(0), C.origin(), std::get<2>(C.sizes()),
                            C.stride(0));
  Tensor3D<T2> ref({nk, nj, ni}, T2(-3.0, -6.0), alloc_b);
  verify_approx(C, ref);
}

TEST_CASE("element_wise_Aij_Bjk_Ckji", "[Numerics][batched_operations]")
{
  test_Aij_Bjk_Ckji<double, std::complex<double>>();
  test_Aij_Bjk_Ckji<float, std::complex<float>>();
  test_Aij_Bjk_Ckji<std::complex<double>, std::complex<double>>();
  test_Aij_Bjk_Ckji<std::complex<float>, std::complex<float>>();
}

TEST_CASE("inplace_product", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  Alloc<double> dalloc{};
  int nb = 4;
  int ni = 3;
  int nj = 2;
  Tensor3D<ComplexType> A({nb, ni, nj}, ComplexType(1.0, -3.0), alloc);
  Tensor2D<double> B({ni, nj}, 2.0, dalloc);
  using ma::inplace_product;
  inplace_product(nb, ni, nj, B.origin(), std::get<1>(B.sizes()), A.origin(), std::get<2>(A.sizes()));
  Tensor3D<ComplexType> ref({nb, ni, nj}, ComplexType(2.0, -6.0), alloc);
  verify_approx(A, ref);
}

// Not clear.
//TEST_CASE("batched_ab_ba", "[Numerics][batched_operations]")
//{

//Alloc<ComplexType> alloc{};
//int nbatch = 4;
//int nrows = 7;
//int ncols = 7;
//std::vector<int> packed_dims = {7,7,7,7,7,7,7,7};
//Tensor1D<ComplexType> Buff;
//Buff = std::move(Tensor1D<ComplexType>(iextensions<1u>{2*nrows*ncols+nbatch}, alloc));
//Tensor2D_ref<ComplexType> A(make_device_ptr(Buff.origin()), {nrows,ncols});
//Tensor2D_ref<ComplexType> B(make_device_ptr(Buff.origin()+nrows*ncols), {nrows,ncols});
//Tensor1D_ref<ComplexType> C(make_device_ptr(Buff.origin()+2*nrows*ncols), iextensions<1u>{nbatch});

////Tensor2D<ComplexType> A({nrows, ncols}, alloc);
////Tensor2D<ComplexType> B({ncols, nrows}, alloc);
//std::vector<pointer<ComplexType>> Aarray(nbatch);
//std::vector<pointer<ComplexType>> Barray(nbatch);
////std::vector<pointer<ComplexType>> y(nbatch);
//std::vector<ComplexType> buffer(nrows*ncols);
//create_data(buffer, ComplexType(1.0));
//copy_n(buffer.data(), buffer.size(), A.origin());
//copy_n(buffer.data(), buffer.size(), B.origin());
//for (int i=0; i < nbatch; i++) {
//Aarray.emplace_back(ma::pointer_dispatch(A.origin()));
//Barray.emplace_back(ma::pointer_dispatch(B.origin()));
//}
//using ma::batched_ab_ba;
//batched_ab_ba(packed_dims.data(), A.data(), ncols, B.data(), nrows,
//ComplexType(1.0), C.data(), nbatch);
//}


TEST_CASE("vbias_from_v1", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  int nkpts              = 4;
  int nwalk              = 4;
  int nchol_max          = 87;
  int nchol_tot          = 339;
  Tensor1D<int> qmap     = {1, 2, 3, 4};
  Tensor1D<int> kminus   = {0, 1, 2, 3};
  Tensor1D<int> nchol_pq = {87, 84, 84, 84};
  Tensor1D<int> q2vbias  = {0, 174, 342, 510};
  Tensor3D<ComplexType> v1({2 * nkpts, nchol_max, nwalk}, ComplexType(1.0, -0.007), alloc);
  Tensor2D<ComplexType> vbias({2 * nchol_tot, nwalk}, 0.0, alloc);
  using ma::vbias_from_v1;
  vbias_from_v1(nwalk, nkpts, nchol_max, qmap.origin(), kminus.origin(), nchol_pq.origin(), q2vbias.origin(),
                ComplexType(0.05), v1.origin(), to_address(vbias.origin()));
  array<ComplexType, 2> vbias_host({2 * nchol_tot, nwalk}, 0.0);
  copy_n(vbias.origin(), vbias.num_elements(), vbias_host.origin());
  // Captured from stdout
  CHECK(real(vbias_host[0][0]) == Approx(0.1));
  CHECK(imag(vbias_host[10][0]) == Approx(-0.0007));
  CHECK(imag(vbias_host[100][2]) == Approx(0));
}

} // namespace qmcplusplus
