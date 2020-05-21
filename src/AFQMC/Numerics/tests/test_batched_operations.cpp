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
#include "AFQMC/Numerics/batched_operations.hpp"
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

// Test data for assertions was created from captured output.
TEST_CASE("Tab_to_Kl", "[Numerics][batched_operations]")
{

  Alloc<ComplexType> alloc{};
  // Tab_to_Kl tests, normally used for Coulomb energy evaluation via
  // E_c[w] = \sum_{na} T[w,aa,n] T[w,aa,n]
  //        = \sum_{n} K[w,n] K[w,n]
  int nwalk = 3;
  int nel = 3;
  int nchol = 7;
  std::vector<ComplexType> buffer(nwalk*nel*nel*nchol);
  create_data(buffer, ComplexType(1.0));
  Tensor4D<ComplexType> Twban({nwalk, nel, nel, nchol}, alloc);
  copy_n(buffer.data(), buffer.size(), Twban.origin());
  Tensor2D<ComplexType> Kl({nwalk,nchol}, 0.0, alloc);
  using ma::Tab_to_Kl;
  Tab_to_Kl(nwalk, nel, nchol, Twban.origin(), Kl.origin());
  array_ref<ComplexType,1,pointer<ComplexType>> Kl_(Kl.origin(), iextensions<1u>{nwalk*nchol});
  array<ComplexType,1,Alloc<ComplexType>> ref = {84.0, 87.0, 90.0, 93.0, 96.0, 99.0, 102.0, 273.0, 276.0,
                            279.0, 282.0, 285.0, 288.0, 291.0, 462.0, 465.0, 468.0, 471.0,
                            474.0, 477.0, 480.0};
  verify_approx(Kl_, ref);

}

TEST_CASE("batched_Tab_to_Klr", "[batched_operations]")
{

  Alloc<ComplexType> alloc{};
  int nwalk = 3;
  int nel = 3;
  int nchol = 11;
  int nbatch = 4;
  int nchol_max = 3;
  int ncholQ = 2;
  int ncholQ0 = 2;
  std::vector<int> kdiag = {0,1,2,3};
  Tensor1D<int> dev_kdiag(iextensions<1u>{nbatch}, alloc);
  copy_n(kdiag.data(), kdiag.size(), dev_kdiag.origin());
  std::vector<ComplexType> buffer(2*nbatch*nwalk*nel*nel*nchol_max);
  create_data(buffer, ComplexType(1.0));
  Tensor5D<ComplexType> Tab({2*nbatch, nwalk, nel, nel, nchol_max}, alloc);
  copy_n(buffer.data(), buffer.size(), Tab.origin());
  Tensor2D<ComplexType> Kl({nwalk,2*nchol}, 0.0, alloc);
  Tensor2D<ComplexType> Kr({nwalk,2*nchol}, 0.0, alloc);
  using ma::batched_Tab_to_Klr;
  batched_Tab_to_Klr(nbatch, nwalk, nel, nchol_max, nchol, ncholQ, ncholQ0,
                     dev_kdiag.origin(), Tab.origin(), Kl.origin(), Kr.origin());
  copy_n(Kr.origin(), Kr.num_elements(), buffer.data());
  //std::cout << std::setprecision(16) << Kl[2][3] << " " << Kl[1][4] << " " << Kr[1][3] << " " << Kr[0][1] << std::endl;
  REQUIRE(real(buffer[2*nchol+3]) == Approx(2262));

}

TEST_CASE("Tanb_to_Kl", "[batched_operations]")
{

  Alloc<ComplexType> alloc{};
  // Tab_to_Kl tests, normally used for Coulomb energy evaluation via
  // E_c[w] = \sum_{na} T[w,aa,n] T[w,aa,n]
  //        = \sum_{n} K[w,n] K[w,n]
  int nwalk = 3;
  int nel = 3;
  int nchol = 7;
  std::vector<ComplexType> buffer(nwalk*nel*nel*nchol);
  create_data(buffer, ComplexType(1.0));
  Tensor4D<ComplexType> Twanb({nwalk, nel, nchol, nel}, alloc);
  copy_n(buffer.data(), buffer.size(), Twanb.origin());
  Tensor2D<ComplexType> Kl({nwalk,nchol}, 0.0, alloc);
  using ma::Tanb_to_Kl;
  Tanb_to_Kl(nwalk, nel, nchol, nchol, Twanb.origin(), Kl.origin());
  array_ref<ComplexType,1,pointer<ComplexType>> Kl_(Kl.origin(), iextensions<1u>{nwalk*nchol});
  //std::cout << "{";
  //for (auto i : Kl_)
    //std::cout << "ComplexType(" << real(i) << ")," << std::endl;
  //std::cout << "};" << std::endl;;
  array<ComplexType,1,Alloc<ComplexType>> ref =  {ComplexType(66), ComplexType(75),
    ComplexType(84), ComplexType(93), ComplexType(102), ComplexType(111),
    ComplexType(120), ComplexType(255), ComplexType(264), ComplexType(273),
    ComplexType(282), ComplexType(291), ComplexType(300), ComplexType(309),
    ComplexType(444), ComplexType(453), ComplexType(462), ComplexType(471),
    ComplexType(480), ComplexType(489), ComplexType(498)};
  verify_approx(Kl_, ref);

}

TEST_CASE("batched_dot_wabn_wban", "[Numerics][batched_operations]")
{

  Alloc<ComplexType> alloc{};
  int nbatch = 4;
  int nwalk = 3;
  int nocc = 7;
  int nchol = 13;
  // 2*nbatch as routine expects left and right tensors to contract.
  std::vector<ComplexType> buffer(2*nbatch*nwalk*nocc*nocc*nchol);
  create_data(buffer, ComplexType(1e4));
  Tensor5D<ComplexType> Twabn({2*nbatch, nwalk, nocc, nocc, nchol}, alloc);
  Tensor1D<ComplexType> scal({nbatch}, 1.0, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  std::vector<pointer<ComplexType>> Aarray;
  array<ComplexType,1,Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  array<ComplexType,1,Alloc<ComplexType>> ref = {ComplexType(1693.073254599684),
                                                 ComplexType(1930.853888599637),
                                                 ComplexType(2189.312510839587)};
  using ma::batched_dot_wabn_wban;
  batched_dot_wabn_wban(nbatch, nwalk, nocc, nchol, scal.origin(),
                        Twabn.origin(), to_address(out.data()), 1);
  //std::cout << std::setprecision(16) << "this: " <<  out[0] << " " << out[1] << " " << out[2] << std::endl;
  verify_approx(ref,out);
}

// Not used.
TEST_CASE("batched_dot_wanb_wbna", "[Numerics][batched_operations]")
{

  Alloc<ComplexType> alloc{};
  int nbatch = 4;
  int nwalk = 3;
  int nocc = 7;
  int nchol = 13;
  // 2*nbatch as routine expects left and right tensors to contract.
  std::vector<ComplexType> buffer(2*nbatch*nwalk*nocc*nocc*nchol);
  create_data(buffer, ComplexType(1e4));
  Tensor5D<ComplexType> Twabn({2*nbatch, nwalk, nocc, nchol, nocc}, alloc);
  Tensor1D<ComplexType> scal({nbatch}, 1.0, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  std::vector<pointer<ComplexType>> Aarray;
  array<ComplexType,1,Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  array<ComplexType,1,Alloc<ComplexType>> ref = {ComplexType(1692.867783879684),
                                                 ComplexType(1930.648417879638),
                                                 ComplexType(2189.107040119586)};
  using ma::batched_dot_wanb_wbna;
  batched_dot_wanb_wbna(nbatch, nwalk, nocc, nchol, scal.origin(),
                        Twabn.origin(), to_address(out.data()), 1);
  //std::cout << std::setprecision(16) << "this: " <<  out[0] << " " << out[1] << " " << out[2] << std::endl;
  verify_approx(ref,out);
}


TEST_CASE("dot_wabn", "[Numerics][batched_operations]")
{

  Alloc<ComplexType> alloc{};
  int nwalk = 3;
  int nocc = 7;
  int nchol = 13;
  std::vector<ComplexType> buffer(nwalk*nocc*nocc*nchol);
  create_data(buffer, ComplexType(100));
  Tensor4D<ComplexType> Twabn({nwalk, nocc, nocc, nchol}, alloc);
  copy_n(buffer.data(), buffer.size(), Twabn.origin());
  array<ComplexType,1,Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  using ma::dot_wabn;
  dot_wabn(nwalk, nocc, nchol, ComplexType(1.0), Twabn.origin(), to_address(out.origin()), 1);
  array<ComplexType,1,Alloc<ComplexType>> ref = {ComplexType(7045.35), ComplexType(58699.7), ComplexType(162049.0)};
  verify_approx(out, ref);
}

TEST_CASE("dot_wanb", "[Numerics][batched_operations]")
{

  Alloc<ComplexType> alloc{};
  int nwalk = 3;
  int nocc = 7;
  int nchol = 13;
  std::vector<ComplexType> buffer(nwalk*nocc*nocc*nchol);
  create_data(buffer, ComplexType(100));
  Tensor4D<ComplexType> Twanb({nwalk, nocc, nchol, nocc}, alloc);
  copy_n(buffer.data(), buffer.size(), Twanb.origin());
  array<ComplexType,1,Alloc<ComplexType>> out(iextensions<1u>{nwalk}, alloc);
  using ma::dot_wanb;
  dot_wanb(nwalk, nocc, nchol, ComplexType(1.0), Twanb.origin(), to_address(out.origin()), 1);
  //std::cout << "{";
  //for (auto i : out)
    //std::cout << "ComplexType(" << real(i) << ")," << std::endl;
  //std::cout << "};" << std::endl;;
  array<ComplexType,1,Alloc<ComplexType>> ref = {ComplexType(6531.67), ComplexType(58186.1), ComplexType(161535)};
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
  Awiu_Biu_Cuw(nu, nw, nn, alpha, A.origin(), B.origin(), B.size(1), C.origin(), C.size(1));
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
  Aijk_Bkj_Cik(ni, nj, nk, A.origin(), A.size(1), A.stride(0),
               B.origin(), B.stride(0), C.origin(), C.stride(0));
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
  std::vector<ComplexType> buffer(ni*nj*nw);
  create_data(buffer, ComplexType(1.0));
  copy_n(buffer.data(), buffer.size(), B.origin());
  using ma::viwj_vwij;
  viwj_vwij(nw, ni, 0, nj, B.data(), A.data());
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
    element_wise_Aij_Bjk_Ckij('N', ni, nj, nk,
                              A.origin(), A.stride(0),
                              B.origin(), B.stride(0),
                              C.origin(), C.size(1), C.size(2));
    Tensor3D<ComplexType> ref({nk, ni, nj}, 6.0, alloc);
    verify_approx(C, ref);
  }
  {
    Tensor2D<ComplexType> A({ni, nj}, ComplexType(0.0,-3.0), alloc);
    Tensor2D<ComplexType> B({nj, nk}, ComplexType(1.0,2.0), alloc);
    Tensor3D<ComplexType> C({nk, ni, nj}, 0.0, alloc);
    element_wise_Aij_Bjk_Ckij('C', ni, nj, nk,
                              A.origin(), A.stride(0),
                              B.origin(), B.stride(0),
                              C.origin(), C.size(1), C.size(2));
    Tensor3D<ComplexType> ref({nk, ni, nj}, ComplexType(-6.0,3.0), alloc);
    verify_approx(C, ref);
  }
}

TEST_CASE("element_wise_Aij_Bjk_Ckji", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  Alloc<double> dalloc{};
  int ni = 3;
  int nj = 2;
  int nk = 2;
  using ma::element_wise_Aij_Bjk_Ckij;
  {
    Tensor2D<ComplexType> A({ni, nj}, -3.0, alloc);
    Tensor2D<ComplexType> B({nj, nk}, ComplexType(1.0,2.0), alloc);
    Tensor3D<ComplexType> C({nk, nj, ni}, 0.0, alloc);
    element_wise_Aij_Bjk_Ckji(ni, nj, nk,
                              A.origin(), A.stride(0),
                              B.origin(), B.stride(0),
                              C.origin(), C.size(2), C.stride(0));
    Tensor3D<ComplexType> ref({nk, nj, ni}, ComplexType(-3.0, -6.0), alloc);
    verify_approx(C, ref);
  }
  {
    Tensor2D<double> A({ni, nj}, -3.0, dalloc);
    Tensor2D<ComplexType> B({nj, nk}, ComplexType(1.0,2.0), alloc);
    Tensor3D<ComplexType> C({nk, nj, ni}, 0.0, alloc);
    element_wise_Aij_Bjk_Ckji(ni, nj, nk,
                              A.origin(), A.stride(0),
                              B.origin(), B.stride(0),
                              C.origin(), C.size(2), C.stride(0));
    Tensor3D<ComplexType> ref({nk, nj, ni}, ComplexType(-3.0, -6.0), alloc);
    verify_approx(C, ref);
  }
}

TEST_CASE("inplace_product", "[Numerics][batched_operations]")
{
  Alloc<ComplexType> alloc{};
  Alloc<double> dalloc{};
  int nb = 4;
  int ni = 3;
  int nj = 2;
  Tensor3D<ComplexType> A({nb, ni, nj}, ComplexType(1.0,-3.0), alloc);
  Tensor2D<double> B({ni, nj}, 2.0, dalloc);
  using ma::inplace_product;
  inplace_product(nb, ni, nj, B.origin(), B.size(1), A.origin(), A.size(2));
  Tensor3D<ComplexType> ref({nb, ni, nj}, ComplexType(2.0,-6.0), alloc);
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


// Buggy.
TEST_CASE("vbias_from_v1", "[Numerics][batched_operations]")
{
  //Alloc<ComplexType> alloc{};
  //int nkpts = 3;
  //int nwalk = 4;
  //int nchol_max = 7;
  //int nchol_tot = 7+3+5;
  //std::vector<int> nchol_pq = {7, 3, 5};
  //std::vector<int> Qmap = {1, 2, 3};
  //std::vector<int> Q2vbias = {0, 14, 14+6};
  //std::vector<int> kminus = {0, 1, 2};
  //std::vector<ComplexType> buffer(nkpts*nchol_max*nwalk);
  //Tensor3D<ComplexType> v1({nkpts, nchol_max, nwalk}, alloc);
  //create_data(buffer, ComplexType(100));
  //copy_n(buffer.data(), buffer.size(), v1.origin());
  //Tensor2D<ComplexType> vbias({2*nchol_tot,nwalk}, 0.0, alloc);
  //Tensor1D<int> dev_nchol_pq(iextensions<1u>{nkpts}, alloc);
  //Tensor1D<int> dev_Qmap(iextensions<1u>{nkpts}, alloc);
  //Tensor1D<int> dev_Q2vbias(iextensions<1u>{nkpts}, alloc);
  //Tensor1D<int> dev_kminus(iextensions<1u>{nkpts}, alloc);
  //copy_n(nchol_pq.data(), nchol_pq.size(), dev_nchol_pq.origin());
  //copy_n(Qmap.data(), Qmap.size(), dev_Qmap.origin());
  //copy_n(Q2vbias.data(), Q2vbias.size(), dev_Q2vbias.origin());
  //copy_n(kminus.data(), kminus.size(), dev_kminus.origin());
  //using ma::vbias_from_v1;
  //vbias_from_v1(nwalk, nkpts, nchol_max, dev_Qmap.origin(), dev_kminus.origin(),
                //dev_nchol_pq.origin(), dev_Q2vbias.origin(), ComplexType(1.0),
                //v1.origin(), to_address(vbias.origin()));
  //copy_n(vbias.origin(), buffer.size(), buffer.data());
  //REQUIRE(real(buffer[1]) == Approx(0.01));
  //REQUIRE(imag(buffer[47]) == Approx(-0.12));
}

}
