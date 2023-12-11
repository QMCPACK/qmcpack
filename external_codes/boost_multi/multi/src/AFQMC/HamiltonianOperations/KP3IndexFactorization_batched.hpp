//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_BATCHED_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_KP3INDEXFACTORIZATION_BATCHED_HPP

#include <vector>
#include <type_traits>
#include <random>
#include <algorithm>

#include "Configuration.h"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/tensor_operations.hpp"


namespace qmcplusplus
{
namespace afqmc
{
// testing the use of dynamic data transfer during execution to reduce memory in GPU
// when an approach is found, integrate in original class through additional template parameter

template<class LQKankMatrix>
class KP3IndexFactorization_batched
{
  // allocators
  using Allocator          = device_allocator<ComplexType>;
  using SpAllocator        = device_allocator<SPComplexType>;
  using BAllocator         = device_allocator<bool>;
  using IAllocator         = device_allocator<int>;
  using Allocator_shared   = node_allocator<ComplexType>;
  using SpAllocator_shared = node_allocator<SPComplexType>;
  using IAllocator_shared  = node_allocator<int>;

  using device_alloc_type  = DeviceBufferManager::template allocator_t<SPComplexType>;
  using device_alloc_Itype = DeviceBufferManager::template allocator_t<int>;

  // type defs
  using pointer                 = typename Allocator::pointer;
  using const_pointer           = typename Allocator::const_pointer;
  using sp_pointer              = typename SpAllocator::pointer;
  using const_sp_pointer        = typename SpAllocator::const_pointer;
  using pointer_shared          = typename Allocator_shared::pointer;
  using const_pointer_shared    = typename Allocator_shared::const_pointer;
  using sp_pointer_shared       = typename SpAllocator_shared::pointer;
  using const_sp_pointer_shared = typename SpAllocator_shared::const_pointer;

  using stdIVector = boost::multi::array<int, 1>;

  using IVector    = boost::multi::array<int, 1, IAllocator>;
  using BoolMatrix = boost::multi::array<bool, 2, BAllocator>;
  using CVector    = ComplexVector<Allocator>;
  using IMatrix    = IntegerMatrix<IAllocator>;
  using CMatrix    = ComplexMatrix<Allocator>;
  using C3Tensor   = boost::multi::array<ComplexType, 3, Allocator>;

  using SpVector  = SPComplexVector<SpAllocator>;
  using SpMatrix  = SPComplexMatrix<SpAllocator>;
  using Sp3Tensor = boost::multi::array<SPComplexType, 3, SpAllocator>;

  using CMatrix_cref  = boost::multi::array_ref<ComplexType const, 2, const_pointer>;
  using CVector_ref   = ComplexVector_ref<pointer>;
  using CMatrix_ref   = ComplexMatrix_ref<pointer>;
  using C3Tensor_ref  = Complex3Tensor_ref<pointer>;
  using C4Tensor_ref  = ComplexArray_ref<4, pointer>;
  using C3Tensor_cref = boost::multi::array_ref<ComplexType const, 3, const_pointer>;

  using SpMatrix_cref = boost::multi::array_ref<SPComplexType const, 2, sp_pointer>;
  using SpVector_ref  = SPComplexVector_ref<sp_pointer>;
  using SpMatrix_ref  = SPComplexMatrix_ref<sp_pointer>;
  using Sp3Tensor_ref = SPComplex3Tensor_ref<sp_pointer>;
  using Sp4Tensor_ref = SPComplexArray_ref<4, sp_pointer>;
  using Sp5Tensor_ref = SPComplexArray_ref<5, sp_pointer>;

  using StaticIVector = boost::multi::static_array<int, 1, device_alloc_Itype>;
  using StaticVector  = boost::multi::static_array<SPComplexType, 1, device_alloc_type>;
  using StaticMatrix  = boost::multi::static_array<SPComplexType, 2, device_alloc_type>;
  using Static3Tensor = boost::multi::static_array<SPComplexType, 3, device_alloc_type>;
  using Static4Tensor = boost::multi::static_array<SPComplexType, 4, device_alloc_type>;

  using shmCVector  = ComplexVector<Allocator_shared>;
  using shmCMatrix  = ComplexMatrix<Allocator_shared>;
  using shmIMatrix  = IntegerMatrix<IAllocator_shared>;
  using shmC3Tensor = Complex3Tensor<Allocator_shared>;

  using mpi3C3Tensor = Complex3Tensor<shared_allocator<ComplexType>>;

  using shmSpVector  = SPComplexVector<SpAllocator_shared>;
  using shmSpMatrix  = SPComplexMatrix<SpAllocator_shared>;
  using shmSp3Tensor = SPComplex3Tensor<SpAllocator_shared>;

public:
  static const HamiltonianTypes HamOpType = KPFactorized;
  HamiltonianTypes getHamType() const { return HamOpType; }

  // NOTE: careful with nocc_max, not consistently defined!!!

  // since arrays can be in host, can't assume that types are consistent
  template<class shmCMatrix_, class shmSpMatrix_>
  KP3IndexFactorization_batched(WALKER_TYPES type,
                                afqmc::TaskGroup_& tg_,
                                stdIVector&& nopk_,
                                stdIVector&& ncholpQ_,
                                stdIVector&& kminus_,
                                boost::multi::array<int, 2>&& nelpk_,
                                boost::multi::array<int, 2>&& QKToK2_,
                                mpi3C3Tensor&& hij_,
                                shmCMatrix_&& h1,
                                std::vector<shmSpMatrix_>&& vik,
                                std::vector<shmSpMatrix_>&& vak,
                                std::vector<shmSpMatrix_>&& vakn,
                                std::vector<shmSpMatrix_>&& vbl,
                                std::vector<shmSpMatrix_>&& vbln,
                                stdIVector&& qqm_,
                                mpi3C3Tensor&& vn0_,
                                std::vector<RealType>&& gQ_,
                                int nsampleQ_,
                                ValueType e0_,
                                Allocator const& alloc_,
                                int cv0,
                                int gncv,
                                int bf_size = 4096)
      : TG(tg_),
        allocator_(alloc_),
        sp_allocator_(alloc_),
        device_buffer_manager(),
        walker_type(type),
        global_nCV(gncv),
        global_origin(cv0),
        default_buffer_size_in_MB(bf_size),
        last_nw(-1),
        E0(e0_),
        H1(std::move(hij_)),
        haj(std::move(h1)),
        nopk(std::move(nopk_)),
        ncholpQ(std::move(ncholpQ_)),
        kminus(std::move(kminus_)),
        nelpk(std::move(nelpk_)),
        QKToK2(std::move(QKToK2_)),
        LQKikn(std::move(move_vector<shmSpMatrix>(std::move(vik)))),
        //LQKank(std::move(move_vector<LQKankMatrix>(std::move(vak),TG.Node()))),
        LQKank(std::move(move_vector<LQKankMatrix>(std::move(vak)))),
        //needs_copy(true),
        needs_copy(not std::is_same<decltype(ma::pointer_dispatch(LQKank[0].origin())), sp_pointer>::value),
        LQKakn(std::move(move_vector<shmSpMatrix>(std::move(vakn)))),
        LQKbnl(std::move(move_vector<shmSpMatrix>(std::move(vbl)))),
        LQKbln(std::move(move_vector<shmSpMatrix>(std::move(vbln)))),
        Qmap(std::move(qqm_)),
        Q2vbias(Qmap.size()),
        vn0(std::move(vn0_)),
        nsampleQ(nsampleQ_),
        gQ(std::move(gQ_)),
        Qwn({1, 1}),
        generator(),
        distribution(gQ.begin(), gQ.end()),
        KKTransID({nopk.size(), nopk.size()}, IAllocator{allocator_}),
        dev_nopk(nopk),
        dev_i0pk(typename IVector::extensions_type{nopk.size()}, IAllocator{allocator_}),
        dev_kminus(kminus),
        dev_ncholpQ(ncholpQ),
        dev_Q2vbias(typename IVector::extensions_type{nopk.size()}, IAllocator{allocator_}),
        dev_Qmap(Qmap),
        dev_nelpk(nelpk),
        dev_a0pk(typename IMatrix::extensions_type{std::get<0>(nelpk.sizes()), std::get<1>(nelpk.sizes())}, IAllocator{allocator_}),
        dev_QKToK2(QKToK2),
        EQ(nopk.size() + 2)
  {
    using std::copy_n;
    using std::fill_n;
    nocc_max = *std::max_element(nelpk.origin(), nelpk.origin() + nelpk.num_elements());
    fill_n(EQ.data(), EQ.size(), 0);
    int nkpts = nopk.size();
    // Defines behavior over Q vector:
    //   <0: Ignore (handled by another TG)
    //    0: Calculate, without rho^+ contribution
    //   >0: Calculate, with rho^+ contribution. LQKbln data located at Qmap[Q]-1
    number_of_symmetric_Q = 0;
    number_of_Q_points    = 0;
    local_nCV             = 0;
    std::fill_n(Q2vbias.origin(), nkpts, -1);
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Q > kminus[Q])
      {
        if (Qmap[kminus[Q]] == 0)
        {
          assert(Qmap[Q] == 0);
          Q2vbias[Q] = 2 * local_nCV;
          local_nCV += ncholpQ[Q];
        }
        else
        {
          assert(Qmap[kminus[Q]] < 0);
          assert(Qmap[Q] < 0);
        }
      }
      else if (Qmap[Q] >= 0)
      {
        Q2vbias[Q] = 2 * local_nCV;
        local_nCV += ncholpQ[Q];
        if (Qmap[Q] > 0)
          number_of_symmetric_Q++;
      }
    }
    for (int Q = 0; Q < nkpts; Q++)
    {
      if (Qmap[Q] >= 0)
        number_of_Q_points++;
      if (Qmap[Q] > 0)
      {
        assert(Q == kminus[Q]);
        assert(Qmap[Q] <= number_of_symmetric_Q);
      }
    }
    copy_n(Q2vbias.data(), nkpts, dev_Q2vbias.origin());
    // setup dev integer arrays
    std::vector<int> i0(nkpts);
    // dev_nopk
    i0[0] = 0;
    for (int i = 1; i < nkpts; i++)
      i0[i] = i0[i - 1] + nopk[i - 1];
    copy_n(i0.data(), nkpts, dev_i0pk.origin());
    // dev_nelpk
    for (int n = 0; n < nelpk.size(); n++)
    {
      i0[0] = 0;
      for (int i = 1; i < nkpts; i++)
        i0[i] = i0[i - 1] + nelpk[n][i - 1];
      copy_n(i0.data(), nkpts, dev_a0pk[n].origin());
      if (walker_type == COLLINEAR)
      {
        i0[0] = 0;
        for (int i = 1; i < nkpts; i++)
          i0[i] = i0[i - 1] + nelpk[n][nkpts + i - 1];
        copy_n(i0.data(), nkpts, dev_a0pk[n].origin() + nkpts);
      }
    }
    // setup copy/transpose tags
    // 1: copy from [Ki][Kj] without rho^+ term
    // 2: transpose from [Ki][Kj] without rho^+ term
    // 3: ignore
    // -P: copy from [Ki][Kj] and transpose from [nkpts+P-1][]
    boost::multi::array<int, 2> KKid({nkpts, nkpts});
    std::fill_n(KKid.origin(), KKid.num_elements(), 3); // ignore everything by default
    for (int Q = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Qmap[Q] < 0)
        continue;
      if (Qmap[Q] > 0)
      { // both rho and rho^+
        assert(Q == kminus[Q]);
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK      = QKToK2[Q][K];
          KKid[K][QK] = -Qmap[Q];
        }
      }
      else if (Q <= kminus[Q])
      {
        // since Qmap[Q]==0 here, Q==kminus[Q] means a hermitian L_ik
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK      = QKToK2[Q][K];
          KKid[K][QK] = 1;
        }
      }
      else if (Q > kminus[Q])
      { // use L(-Q)(ki)*
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK      = QKToK2[Q][K];
          KKid[K][QK] = 2;
        }
      }
    }
    copy_n(KKid.origin(), KKid.num_elements(), KKTransID.origin());

    long memank = 0;
    if (needs_copy)
      for (auto& v : LQKank)
        memank = std::max(memank, 2 * v.num_elements());
    else
      for (auto& v : LQKank)
        memank += v.num_elements();

    // report memory usage
    size_t likn(0), lakn(0), lbln(0), misc(0);
    for (auto& v : LQKikn)
      likn += v.num_elements();
    for (auto& v : LQKakn)
      lakn += v.num_elements();
    for (auto& v : LQKbln)
      lbln += v.num_elements();
    for (auto& v : LQKbnl)
      lbln += v.num_elements();
    app_log() << "****************************************************************** \n";
    if (needs_copy)
      app_log() << "  Using out of core storage of LQKakn \n";
    else
      app_log() << "  Using device storage of LQKakn \n";
    app_log() << "  Static memory usage by KP3IndexFactorization_batched (node 0 in MB) \n"
              << "    L[Q][K][ikn]: " << likn * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n"
              << "    L[Q][K][akn]: " << (lakn + memank) * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n"
              << "    L[Q][K][bln]: " << lbln * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n";
    memory_report();
  }

  ~KP3IndexFactorization_batched() {}

  KP3IndexFactorization_batched(const KP3IndexFactorization_batched& other) = delete;
  KP3IndexFactorization_batched& operator=(const KP3IndexFactorization_batched& other) = delete;
  KP3IndexFactorization_batched(KP3IndexFactorization_batched&& other)                 = default;
  KP3IndexFactorization_batched& operator=(KP3IndexFactorization_batched&& other) = default;

  // must have the same signature as shared classes, so keeping it with std::allocator
  // NOTE: THIS SHOULD USE mpi3::shm!!!
  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(TaskGroup_& TG_,
                                                                 boost::multi::array<ComplexType, 1> const& vMF)
  {
    int nkpts = nopk.size();
    int NMO   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int npol  = (walker_type == NONCOLLINEAR) ? 2 : 1;

    CVector vMF_(vMF);
    CVector P0D(iextensions<1u>{NMO * NMO});
    fill_n(P0D.origin(), P0D.num_elements(), ComplexType(0));
    vHS(vMF_, P0D);
    if (TG_.TG().size() > 1)
      TG_.TG().all_reduce_in_place_n(to_address(P0D.origin()), P0D.num_elements(), std::plus<>());

    boost::multi::array<ComplexType, 2> P0({NMO, NMO});
    copy_n(P0D.origin(), NMO * NMO, P0.origin());

    boost::multi::array<ComplexType, 2> P1({npol * NMO, npol * NMO});
    std::fill_n(P1.origin(), P1.num_elements(), ComplexType(0.0));

    // add spin-dependent H1
    for (int K = 0, nk0 = 0; K < nkpts; ++K)
    {
      for (int i = 0, I = nk0; i < nopk[K]; i++, I++)
      {
        for (int p = 0; p < npol; ++p)
          P1[p * NMO + I][p * NMO + I] += H1[K][p * nopk[K] + i][p * nopk[K] + i];
        for (int j = i + 1, J = I + 1; j < nopk[K]; j++, J++)
        {
          for (int p = 0; p < npol; ++p)
          {
            P1[p * NMO + I][p * NMO + J] += H1[K][p * nopk[K] + i][p * nopk[K] + j];
            P1[p * NMO + J][p * NMO + I] += H1[K][p * nopk[K] + j][p * nopk[K] + i];
          }
        }
        if (walker_type == NONCOLLINEAR)
        {
          // offdiagonal piece
          for (int j = 0, J = nk0; j < nopk[K]; j++, J++)
          {
            P1[I][NMO + J] += H1[K][i][nopk[K] + j];
            P1[NMO + J][I] += H1[K][nopk[K] + j][i];
          }
        }
      }
      nk0 += nopk[K];
    }

    // add P0 (diagonal in spin)
    for (int p = 0; p < npol; ++p)
      for (int I = 0; I < NMO; I++)
        for (int J = 0; J < NMO; J++)
          P1[p * NMO + I][p * NMO + J] += P0[I][J];

    // add vn0 (diagonal in spin)
    for (int K = 0, nk0 = 0; K < nkpts; ++K)
    {
      for (int i = 0, I = nk0; i < nopk[K]; i++, I++)
      {
        for (int p = 0; p < npol; ++p)
          P1[p * NMO + I][p * NMO + I] += vn0[K][i][i];
        for (int j = i + 1, J = I + 1; j < nopk[K]; j++, J++)
        {
          for (int p = 0; p < npol; ++p)
          {
            P1[p * NMO + I][p * NMO + J] += vn0[K][i][j];
            P1[p * NMO + J][p * NMO + I] += vn0[K][j][i];
          }
        }
      }
      nk0 += nopk[K];
    }

    using ma::conj;
    // symmetrize
    for (int I = 0; I < npol * NMO; I++)
    {
      for (int J = I + 1; J < npol * NMO; J++)
      {
        // This is really cutoff dependent!!!
#if defined(MIXED_PRECISION)
        if (std::abs(P1[I][J] - ma::conj(P1[J][I])) * 2.0 > 1e-5)
        {
#else
        if (std::abs(P1[I][J] - ma::conj(P1[J][I])) * 2.0 > 1e-6)
        {
#endif
          app_error() << " WARNING in getOneBodyPropagatorMatrix. H1 is not hermitian. \n";
          app_error() << I << " " << J << " " << P1[I][J] << " " << P1[J][I] << std::endl;
          //<< H1[K][i][j] << " "
          //<< H1[K][j][i] << " " << vn0[K][i][j] << " " << vn0[K][j][i] << std::endl;
        }
        P1[I][J] = 0.5 * (P1[I][J] + ma::conj(P1[J][I]));
        P1[J][I] = ma::conj(P1[I][J]);
      }
    }
    return P1;
  }

  template<class Mat, class MatB>
  void energy(Mat&& E, MatB const& G, int k = 0, bool addH1 = true, bool addEJ = true, bool addEXX = true)
  {
    MatB* Kr(nullptr);
    MatB* Kl(nullptr);
    energy(E, G, k, Kl, Kr, addH1, addEJ, addEXX);
  }

  // KEleft and KEright must be in shared memory for this to work correctly
  template<
      class Mat,
      class MatB,
      class MatC,
      class MatD
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<Mat>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatC>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatD>()))
      >
  void energy(Mat&& E,
              MatB const& Gc,
              int nd,
              MatC* KEleft,
              MatD* KEright,
              bool addH1  = true,
              bool addEJ  = true,
              bool addEXX = true)
  {
    if (nsampleQ > 0)
      energy_sampleQ(E, Gc, nd, KEleft, KEright, addH1, addEJ, addEXX);
    else
      energy_exact(E, Gc, nd, KEleft, KEright, addH1, addEJ, addEXX);
  }

  // KEleft and KEright must be in shared memory for this to work correctly
  template<class Mat, class MatB, class MatC, class MatD>
  void energy_exact(Mat&& E,
                    MatB const& Gc,
                    int nd,
                    MatC* KEleft,
                    MatD* KEright,
                    bool addH1  = true,
                    bool addEJ  = true,
                    bool addEXX = true)
  {
    using std::copy_n;
    using std::fill_n;
    int nkpts = nopk.size();
    assert(std::get<1>(E.sizes()) >= 3);
    assert(nd >= 0 && nd < nelpk.size());

    int nwalk     = std::get<1>(Gc.sizes());
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int npol      = (walker_type == NONCOLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nocca_tot = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + nkpts, 0);
    int nocca_max = *std::max_element(nelpk[nd].begin(), nelpk[nd].begin() + nkpts);
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    int noccb_tot = 0;
    if (walker_type == COLLINEAR)
      noccb_tot = std::accumulate(nelpk[nd].begin() + nkpts, nelpk[nd].begin() + 2 * nkpts, 0);
    int getKr = KEright != nullptr;
    int getKl = KEleft != nullptr;
    if (std::get<0>(E.sizes()) != nwalk || std::get<1>(E.sizes()) < 3)
      APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix "
                "dimensions \n");

    // take from BufferManager.
    //      long default_buffer_size_in_MB(4L*1024L);
    long batch_size(0);
    if (addEXX)
    {
      long Bytes = long(default_buffer_size_in_MB) * 1024L * 1024L;
      Bytes /= size_t(nwalk * nocc_max * nocc_max * nchol_max * sizeof(SPComplexType));
      long bz0 = std::max(2L, Bytes);
      // batch_size includes the factor of 2 from Q/Qm pair
      batch_size = std::min(bz0, long(2 * number_of_Q_points * nkpts));
      // make sure batch_size is even
      batch_size = batch_size - (batch_size % 2L);
      assert(batch_size % 2L == 0);
    }

    long Knr = 0, Knc = 0;
    if (addEJ)
    {
      Knr = nwalk;
      Knc = local_nCV;
      if (getKr)
      {
        assert(std::get<0>(KEright->sizes()) == nwalk && std::get<1>(KEright->sizes()) == local_nCV);
        assert(KEright->stride(0) == std::get<1>(KEright->sizes()));
      }
      if (getKl)
      {
        assert(std::get<0>(KEleft->sizes()) == nwalk && std::get<1>(KEleft->sizes()) == local_nCV);
        assert(KEleft->stride(0) == std::get<1>(KEleft->sizes()));
      }
    }
    else if (getKr or getKl)
    {
      APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
    }
    StaticMatrix Kl({Knr, Knc}, device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    StaticMatrix Kr({Knr, Knc}, device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    fill_n(Kr.origin(), Knr * Knc, SPComplexType(0.0));
    fill_n(Kl.origin(), Knr * Knc, SPComplexType(0.0));

    for (int n = 0; n < nwalk; n++)
      fill_n(E[n].origin(), 3, ComplexType(0.));

    assert(Gc.num_elements() == nwalk * (nocca_tot + noccb_tot) * npol * nmo_tot);
    C3Tensor_cref G3Da(make_device_ptr(Gc.origin()), {nocca_tot * npol, nmo_tot, nwalk});
    C3Tensor_cref G3Db(make_device_ptr(Gc.origin()) + G3Da.num_elements() * (nspin - 1), {noccb_tot, nmo_tot, nwalk});

    // later on, rewrite routine to loop over spins, to avoid storage of both spin
    // components simultaneously
    Static4Tensor GKK({nspin, nkpts, nkpts, nwalk * npol * nmo_max * nocc_max},
                      device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    GKaKjw_to_GKKwaj(G3Da, GKK[0], nelpk[nd].sliced(0, nkpts), dev_nelpk[nd], dev_a0pk[nd]);
    if (walker_type == COLLINEAR)
      GKaKjw_to_GKKwaj(G3Db, GKK[1], nelpk[nd].sliced(nkpts, 2 * nkpts), dev_nelpk[nd].sliced(nkpts, 2 * nkpts),
                       dev_a0pk[nd].sliced(nkpts, 2 * nkpts));
    // one-body contribution
    // haj[ndet*nkpts][nocc*nmo]
    // not parallelized for now, since it would require customization of Wfn
    if (addH1)
    {
      for (int n = 0; n < nwalk; n++)
        fill_n(E[n].origin(), 1, ComplexType(E0));
      // must use Gc since GKK is is SP
      int na = 0, nk = 0, nb = 0;
      for (int K = 0; K < nkpts; ++K)
      {
#if defined(MIXED_PRECISION)
        int ni(nopk[K]);
        CMatrix_ref haj_K(make_device_ptr(haj[nd * nkpts + K].origin()), {nocc_max, npol * nmo_max});
        for (int a = 0; a < nelpk[nd][K]; ++a)
          for (int pol = 0; pol < npol; ++pol)
            ma::product(ComplexType(1.), ma::T(G3Da[(na + a) * npol + pol].sliced(nk, nk + ni)),
                        haj_K[a].sliced(pol * ni, pol * ni + ni), ComplexType(1.), E({0, nwalk}, 0));
        na += nelpk[nd][K];
        if (walker_type == COLLINEAR)
        {
          boost::multi::array_ref<ComplexType, 2, pointer> haj_Kb(haj_K.origin() + haj_K.num_elements(),
                                                                  {nocc_max, nmo_max});
          for (int b = 0; b < nelpk[nd][nkpts + K]; ++b)
            ma::product(ComplexType(1.), ma::T(G3Db[nb + b].sliced(nk, nk + ni)), haj_Kb[b].sliced(0, ni),
                        ComplexType(1.), E({0, nwalk}, 0));
          nb += nelpk[nd][nkpts + K];
        }
        nk += ni;
#else
        nk = nopk[K];
        {
          na = nelpk[nd][K];
          CVector_ref haj_K(make_device_ptr(haj[nd * nkpts + K].origin()), {nocc_max * npol * nmo_max});
          SpMatrix_ref Gaj(GKK[0][K][K].origin(), {nwalk, nocc_max * npol * nmo_max});
          ma::product(ComplexType(1.), Gaj, haj_K, ComplexType(1.), E({0, nwalk}, 0));
        }
        if (walker_type == COLLINEAR)
        {
          na = nelpk[nd][nkpts + K];
          CVector_ref haj_K(make_device_ptr(haj[nd * nkpts + K].origin()) + nocc_max * nmo_max, {nocc_max * nmo_max});
          SpMatrix_ref Gaj(GKK[1][K][K].origin(), {nwalk, nocc_max * nmo_max});
          ma::product(ComplexType(1.), Gaj, haj_K, ComplexType(1.), E({0, nwalk}, 0));
        }
#endif
      }
    }

    // move calculation of H1 here
    // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
    //       Not sure how to do it for COLLINEAR.
    if (addEXX)
    {
      int batch_cnt(0);
      using ma::gemmBatched;
      std::vector<sp_pointer> Aarray;
      std::vector<sp_pointer> Barray;
      std::vector<sp_pointer> Carray;
      Aarray.reserve(batch_size);
      Barray.reserve(batch_size);
      Carray.reserve(batch_size);
      std::vector<SPComplexType> scl_factors;
      scl_factors.reserve(batch_size);
      std::vector<int> kdiag;
      kdiag.reserve(batch_size);

      StaticIVector IMats(iextensions<1u>{batch_size},
                          device_buffer_manager.get_generator().template get_allocator<int>());
      fill_n(IMats.origin(), IMats.num_elements(), 0);
      StaticVector dev_scl_factors(iextensions<1u>{batch_size},
                                   device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      Static3Tensor T1({batch_size, nwalk * nocc_max, nocc_max * nchol_max},
                       device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      SPRealType scl = (walker_type == CLOSED ? 2.0 : 1.0);

      // I WANT C++17!!!!!!
      long mem_ank(0);
      if (needs_copy)
        mem_ank = nkpts * nocc_max * nchol_max * npol * nmo_max;
      StaticVector LBuff(iextensions<1u>{2 * mem_ank},
                         device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      sp_pointer LQptr(nullptr), LQmptr(nullptr);
      if (needs_copy)
      {
        // data will be copied here
        LQptr  = LBuff.origin();
        LQmptr = LBuff.origin() + mem_ank;
      }

      for (int spin = 0; spin < nspin; ++spin)
      {
        for (int Q = 0; Q < nkpts; ++Q)
        {
          if (Qmap[Q] < 0)
            continue;
          bool haveKE = false;
          int Qm      = kminus[Q];

          // simple implementation for now
          Aarray.clear();
          Barray.clear();
          Carray.clear();
          scl_factors.clear();
          kdiag.clear();
          batch_cnt = 0;

          // choose source of data depending on whether data needs to be copied or not
          if (!needs_copy)
          {
            // set to local array origin
            LQptr  = make_device_ptr(LQKank[nd * nspin * nkpts + spin * nkpts + Q].origin());
            LQmptr = make_device_ptr(LQKank[nd * nspin * nkpts + spin * nkpts + Qm].origin());
          }

          SpMatrix_ref LQ(LQptr, LQKank[nd * nspin * nkpts + spin * nkpts + Q].extensions());
          SpMatrix_ref LQm(LQmptr, LQKank[nd * nspin * nkpts + spin * nkpts + Qm].extensions());

          if (needs_copy)
          {
            copy_n(to_address(LQKank[nd * nspin * nkpts + spin * nkpts + Q].origin()), LQ.num_elements(), LQ.origin());
            if (Q != Qm)
              copy_n(to_address(LQKank[nd * nspin * nkpts + spin * nkpts + Qm].origin()), LQm.num_elements(),
                     LQm.origin());
          }

          for (int Ka = 0; Ka < nkpts; ++Ka)
          {
            int K0 = ((Qmap[Q] > 0) ? 0 : Ka);
            for (int Kb = K0; Kb < nkpts; ++Kb)
            {
              int Kl_ = QKToK2[Qm][Kb];
              int Kk  = QKToK2[Q][Ka];

              if (addEJ && Ka == Kb)
                kdiag.push_back(batch_cnt);

              if (Qmap[Q] > 0)
                Aarray.push_back(sp_pointer(
                    LQKbnl[nd * nspin * number_of_symmetric_Q + spin * number_of_symmetric_Q + Qmap[Q] - 1][Kb]
                        .origin()));
              else
                Aarray.push_back(sp_pointer(LQm[Kb].origin()));

              Barray.push_back(GKK[spin][Ka][Kl_].origin());
              Carray.push_back(T1[batch_cnt++].origin());
              Aarray.push_back(sp_pointer(LQ[Ka].origin()));
              Barray.push_back(GKK[spin][Kb][Kk].origin());
              Carray.push_back(T1[batch_cnt++].origin());

              if (Qmap[Q] > 0 || Ka == Kb)
                scl_factors.push_back(SPComplexType(-scl * 0.5));
              else
                scl_factors.push_back(SPComplexType(-scl));

              if (batch_cnt >= batch_size)
              {
                gemmBatched('T', 'N', nocc_max * nchol_max, nwalk * nocc_max, npol * nmo_max, SPComplexType(1.0),
                            Aarray.data(), npol * nmo_max, Barray.data(), npol * nmo_max, SPComplexType(0.0),
                            Carray.data(), nocc_max * nchol_max, Aarray.size());

                copy_n(scl_factors.data(), scl_factors.size(), dev_scl_factors.origin());
                using ma::batched_dot_wabn_wban;
                batched_dot_wabn_wban(scl_factors.size(), nwalk, nocc_max, nchol_max, dev_scl_factors.origin(),
                                      T1.origin(), to_address(E[0].origin()) + 1, E.stride(0));

                if (addEJ)
                {
                  int nc0 = Q2vbias[Q] / 2; //std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
                  copy_n(kdiag.data(), kdiag.size(), IMats.origin());
                  using ma::batched_Tab_to_Klr;
                  batched_Tab_to_Klr(kdiag.size(), nwalk, nocc_max, nchol_max, local_nCV, ncholpQ[Q], nc0,
                                     IMats.origin(), T1.origin(), Kl.origin(), Kr.origin());
                }

                // reset
                Aarray.clear();
                Barray.clear();
                Carray.clear();
                scl_factors.clear();
                kdiag.clear();
                batch_cnt = 0;
              }
            }
          }

          if (batch_cnt > 0)
          {
            gemmBatched('T', 'N', nocc_max * nchol_max, nwalk * nocc_max, npol * nmo_max, SPComplexType(1.0),
                        Aarray.data(), npol * nmo_max, Barray.data(), npol * nmo_max, SPComplexType(0.0), Carray.data(),
                        nocc_max * nchol_max, Aarray.size());

            copy_n(scl_factors.data(), scl_factors.size(), dev_scl_factors.origin());
            using ma::batched_dot_wabn_wban;
            batched_dot_wabn_wban(scl_factors.size(), nwalk, nocc_max, nchol_max, dev_scl_factors.origin(), T1.origin(),
                                  to_address(E[0].origin()) + 1, E.stride(0));

            if (addEJ)
            {
              int nc0 = Q2vbias[Q] / 2; //std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
              copy_n(kdiag.data(), kdiag.size(), IMats.origin());
              using ma::batched_Tab_to_Klr;
              batched_Tab_to_Klr(kdiag.size(), nwalk, nocc_max, nchol_max, local_nCV, ncholpQ[Q], nc0, IMats.origin(),
                                 T1.origin(), Kl.origin(), Kr.origin());
            }
          }
        } // Q
      }   // COLLINEAR
    }

    if (addEJ)
    {
      if (not addEXX)
      {
        // calculate Kr
        APP_ABORT(" Error: Finish addEJ and not addEXX");
      }
      RealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
      using ma::adotpby;
      for (int n = 0; n < nwalk; ++n)
      {
        adotpby(SPComplexType(0.5 * scl * scl), Kl[n], Kr[n], ComplexType(0.0), E[n].origin() + 2);
      }
      if (getKr)
        copy_n_cast(Kr.origin(), Kr.num_elements(), make_device_ptr(KEright->origin()));
      if (getKl)
        copy_n_cast(Kl.origin(), Kl.num_elements(), make_device_ptr(KEleft->origin()));
    }
  }

  // KEleft and KEright must be in shared memory for this to work correctly
  template<class Mat, class MatB, class MatC, class MatD>
  void energy_sampleQ(Mat&& E,
                      MatB const& Gc,
                      int nd,
                      MatC* KEleft,
                      MatD* KEright,
                      bool addH1  = true,
                      bool addEJ  = true,
                      bool addEXX = true)
  {
    APP_ABORT("Error: energy_sampleQ not yet implemented in batched routine.\n");
    /*
      using std::fill_n;
      int nkpts = nopk.size(); 
      assert(E.size(1)>=3);
      assert(nd >= 0 && nd < nelpk.size());  

      int nwalk = Gc.size(1);
      int nspin = (walker_type==COLLINEAR?2:1);
      int nmo_tot = std::accumulate(nopk.begin(),nopk.end(),0);
      int nmo_max = *std::max_element(nopk.begin(),nopk.end());
      int nocca_tot = std::accumulate(nelpk[nd].begin(),nelpk[nd].begin()+nkpts,0);
      int nocca_max = *std::max_element(nelpk[nd].begin(),nelpk[nd].begin()+nkpts);
      int nchol_max = *std::max_element(ncholpQ.begin(),ncholpQ.end());
      int noccb_tot = 0;
      if(walker_type==COLLINEAR) noccb_tot = std::accumulate(nelpk[nd].begin()+nkpts,
                                      nelpk[nd].begin()+2*nkpts,0);
      int getKr = KEright!=nullptr;
      int getKl = KEleft!=nullptr;
      if(E.size(0) != nwalk || E.size(1) < 3)
        APP_ABORT(" Error in AFQMC/HamiltonianOperations/sparse_matrix_energy::calculate_energy(). Incorrect matrix dimensions \n");

      size_t mem_needs(nwalk*nkpts*nkpts*nspin*nocca_max*nmo_max);
      size_t cnt(0);  
      if(addEJ) { 
#if defined(MIXED_PRECISION)
        mem_needs += 2*nwalk*local_nCV;
#else
        if(not getKr) mem_needs += nwalk*local_nCV;
        if(not getKl) mem_needs += nwalk*local_nCV;
#endif
      }
      set_buffer(mem_needs);

      // messy
      sp_pointer Krptr(nullptr), Klptr(nullptr);
      long Knr=0, Knc=0;
      if(addEJ) {
        Knr=nwalk;
        Knc=local_nCV;
        cnt=0;
#if defined(MIXED_PRECISION)
        if(getKr) {
          assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV);
          assert(KEright->stride(0) == KEright->size(1));
        }
#else
        if(getKr) {
          assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV);
          assert(KEright->stride() == KEright->size(1));
          Krptr = make_device_ptr(KEright->origin());
        } else 
#endif
        {
          Krptr = BTMats.origin();
          cnt += nwalk*local_nCV;
        }
#if defined(MIXED_PRECISION)
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
          assert(KEleft->stride(0) == KEleft->size(1));
        }
#else
        if(getKl) {
          assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
          assert(KEleft->stride(0) == KEleft->size(1));
          Klptr = make_device_ptr(KEleft->origin());
        } else 
#endif
        {
          Klptr = BTMats.origin()+cnt;
          cnt += nwalk*local_nCV;
        }
        fill_n(Krptr,Knr*Knc,SPComplexType(0.0));
        fill_n(Klptr,Knr*Knc,SPComplexType(0.0));
      } else if(getKr or getKl) {
        APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
      }
      SpMatrix_ref Kl(Klptr,{Knr,Knc});
      SpMatrix_ref Kr(Krptr,{Knr,Knc});

      for(int n=0; n<nwalk; n++) 
        fill_n(E[n].origin(),3,ComplexType(0.));

      assert(Gc.num_elements() == nwalk*(nocca_tot+noccb_tot)*nmo_tot);
      C3Tensor_cref G3Da(make_device_ptr(Gc.origin()),{nocca_tot,nmo_tot,nwalk} );
      C3Tensor_cref G3Db(make_device_ptr(Gc.origin())+G3Da.num_elements()*(nspin-1),
                            {noccb_tot,nmo_tot,nwalk} );

      Sp4Tensor_ref GKK(BTMats.origin()+cnt,
                        {nspin,nkpts,nkpts,nwalk*nmo_max*nocca_max});
      cnt+=GKK.num_elements();
      GKaKjw_to_GKKwaj(G3Da,GKK[0],nelpk[nd].sliced(0,nkpts),dev_nelpk[nd],dev_a0pk[nd]);
      if(walker_type==COLLINEAR)  
        GKaKjw_to_GKKwaj(G3Db,GKK[1],nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_nelpk[nd].sliced(nkpts,2*nkpts),
                                     dev_a0pk[nd].sliced(nkpts,2*nkpts));

      // one-body contribution
      // haj[ndet*nkpts][nocc*nmo]
      // not parallelized for now, since it would require customization of Wfn 
      if(addH1) {
        // must use Gc since GKK is is SP
        int na=0, nk=0, nb=0;
        for(int n=0; n<nwalk; n++)
          E[n][0] = E0;  
        for(int K=0; K<nkpts; ++K) {
#if defined(MIXED_PRECISION) 
          CMatrix_ref haj_K(make_device_ptr(haj[nd*nkpts+K].origin()),{nocc_max,nmo_max});
          for(int a=0; a<nelpk[nd][K]; ++a)
            ma::product(ComplexType(1.),ma::T(G3Da[na+a].sliced(nk,nk+nopk[K])),
                                        haj_K[a].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
          na+=nelpk[nd][K];
          if(walker_type==COLLINEAR) {
            boost::multi::array_ref<ComplexType,2,pointer> haj_Kb(haj_K.origin()+haj_K.num_elements(),
                                                      {nocc_max,nmo_max});
            for(int b=0; b<nelpk[nd][nkpts+K]; ++b)
              ma::product(ComplexType(1.),ma::T(G3Db[nb+b].sliced(nk,nk+nopk[K])),
                                        haj_Kb[b].sliced(0,nopk[K]),
                        ComplexType(1.),E({0,nwalk},0));
            nb+=nelpk[nd][nkpts+K];
          }
          nk+=nopk[K];
#else
          nk = nopk[K];
          {
            na = nelpk[nd][K];
            CVector_ref haj_K(make_device_ptr(haj[nd*nkpts+K].origin()),{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[0][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }
          if(walker_type==COLLINEAR) {
            na = nelpk[nd][nkpts+K];
            CVector_ref haj_K(make_device_ptr(haj[nd*nkpts+K].origin())+nocc_max*nmo_max,{nocc_max*nmo_max});
            SpMatrix_ref Gaj(GKK[1][K][K].origin(),{nwalk,nocc_max*nmo_max});
            ma::product(ComplexType(1.),Gaj,haj_K,ComplexType(1.),E({0,nwalk},0));
          }
#endif
        }
      }

      // move calculation of H1 here	
      // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
      //       Not sure how to do it for COLLINEAR.
      if(addEXX) {  

        if(Qwn.size(0) != nwalk || Qwn.size(1) != nsampleQ)
          Qwn = std::move(boost::multi::array<int,2>({nwalk,nsampleQ})); 
        {
          for(int n=0; n<nwalk; ++n) 
            for(int nQ=0; nQ<nsampleQ; ++nQ) {
              Qwn[n][nQ] = distribution(generator);
              RealType drand = distribution(generator);
              RealType s(0.0);
              bool found=false;
              for(int Q=0; Q<nkpts; Q++) {
                s += gQ[Q];
                if( drand < s ) {
                  Qwn[n][nQ] = Q;
                  found=true;
                  break;
                }
              } 
              if(not found) 
                APP_ABORT(" Error: sampleQ Qwn. \n");  
            }
        }
        size_t local_memory_needs = 2*nocca_max*nocca_max*nchol_max; 
        if(TMats.num_elements() < local_memory_needs) { 
          TMats = std::move(SpVector(iextensions<1u>{local_memory_needs})); 
          using std::fill_n;
          fill_n(TMats.origin(),TMats.num_elements(),SPComplexType(0.0));
        }
        size_t local_cnt=0; 
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        size_t nqk=1;  
        for(int n=0; n<nwalk; ++n) {
          for(int nQ=0; nQ<nsampleQ; ++nQ) {
            int Q = Qwn[n][nQ];
            for(int Ka=0; Ka<nkpts; ++Ka) {
              for(int Kb=0; Kb<nkpts; ++Kb) {
                { 
                  int nchol = ncholpQ[Q];
                  int Qm = kminus[Q];
                  int Kl = QKToK2[Qm][Kb];
                  int Kk = QKToK2[Q][Ka];
                  int nl = nopk[Kl];
                  int nb = nelpk[nd][Kb];
                  int na = nelpk[nd][Ka];
                  int nk = nopk[Kk];

                  SpMatrix_ref Gal(GKK[0][Ka][Kl].origin()+n*na*nl,{na,nl});
                  SpMatrix_ref Gbk(GKK[0][Kb][Kk].origin()+n*nb*nk,{nb,nk});
                  SpMatrix_ref Lank(sp_pointer(LQKank[nd*nspin*nkpts+Q][Ka].origin()),
                                                 {na*nchol,nk});
                  auto bnl_ptr(sp_pointer(LQKank[nd*nspin*nkpts+Qm][Kb].origin()));
                  if( Q == Qm ) bnl_ptr = sp_pointer(LQKbnl[nd*nspin*number_of_symmetric_Q+Qmap[Q]-1][Kb].origin());
                  SpMatrix_ref Lbnl(bnl_ptr,{nb*nchol,nl});

                  SpMatrix_ref Tban(TMats.origin()+local_cnt,{nb,na*nchol});
                  Sp3Tensor_ref T3Dban(TMats.origin()+local_cnt,{nb,na,nchol});
                  SpMatrix_ref Tabn(Tban.origin()+Tban.num_elements(),{na,nb*nchol});
                  Sp3Tensor_ref T3Dabn(Tban.origin()+Tban.num_elements(),{na,nb,nchol});

                  ma::product(Gal,ma::T(Lbnl),Tabn);
                  ma::product(Gbk,ma::T(Lank),Tban);

                  SPComplexType E_(0.0);
                  for(int a=0; a<na; ++a)
                    for(int b=0; b<nb; ++b)
                      E_ += ma::dot(T3Dabn[a][b],T3Dban[b][a]);
                  E[n][1] -= scl*0.5*static_cast<ComplexType>(E_)/gQ[Q]/double(nsampleQ);

                } // if

                if(walker_type==COLLINEAR) {

                  { 
                    int nchol = ncholpQ[Q];
                    int Qm = kminus[Q];
                    int Kl = QKToK2[Qm][Kb];
                    int Kk = QKToK2[Q][Ka];
                    int nl = nopk[Kl];
                    int nb = nelpk[nd][nkpts+Kb];
                    int na = nelpk[nd][nkpts+Ka];
                    int nk = nopk[Kk];

                    SpMatrix_ref Gal(GKK[1][Ka][Kl].origin()+n*na*nl,{na,nl});
                    SpMatrix_ref Gbk(GKK[1][Kb][Kk].origin()+n*nb*nk,{nb,nk});
                    SpMatrix_ref Lank(sp_pointer(LQKank[(nd*nspin+1)*nkpts+Q][Ka].origin()),
                                                 {na*nchol,nk});
                    auto bnl_ptr(sp_pointer(LQKank[nd*nspin*nkpts+Qm][Kb].origin()));
                    if( Q == Qm ) bnl_ptr = sp_pointer(LQKbnl[(nd*nspin+1)*number_of_symmetric_Q+
                                                                Qmap[Q]-1][Kb].origin());
                    SpMatrix_ref Lbnl(bnl_ptr,{nb*nchol,nl});

                    SpMatrix_ref Tban(TMats.origin()+local_cnt,{nb,na*nchol});
                    Sp3Tensor_ref T3Dban(TMats.origin()+local_cnt,{nb,na,nchol});
                    SpMatrix_ref Tabn(Tban.origin()+Tban.num_elements(),{na,nb*nchol});
                    Sp3Tensor_ref T3Dabn(Tban.origin()+Tban.num_elements(),{na,nb,nchol});
  
                    ma::product(Gal,ma::T(Lbnl),Tabn);
                    ma::product(Gbk,ma::T(Lank),Tban);
  
                    SPComplexType E_(0.0);
                    for(int a=0; a<na; ++a)
                      for(int b=0; b<nb; ++b)
                        E_ += ma::dot(T3Dabn[a][b],T3Dban[b][a]);
                    E[n][1] -= scl*0.5*static_cast<ComplexType>(E_)/gQ[Q]/double(nsampleQ);

                  } // if
                } // COLLINEAR 
              } // Kb 
            } // Ka
          } // nQ
        } // n 
      }  

      if(addEJ) {
        size_t local_memory_needs = 2*nchol_max*nwalk; 
        if(TMats.num_elements() < local_memory_needs) { 
          TMats = std::move(SpVector(iextensions<1u>{local_memory_needs}));
          using std::fill_n;
          fill_n(TMats.origin(),TMats.num_elements(),SPComplexType(0.0));
        }
        cnt=0; 
        SpMatrix_ref Kr_local(TMats.origin(),{nwalk,nchol_max}); 
        cnt+=Kr_local.num_elements();
        SpMatrix_ref Kl_local(TMats.origin()+cnt,{nwalk,nchol_max}); 
        cnt+=Kl_local.num_elements();
        fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
        fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
        size_t nqk=1;  
        for(int Q=0; Q<nkpts; ++Q) {
          bool haveKE=false;
          for(int Ka=0; Ka<nkpts; ++Ka) {
            { 
              haveKE=true;
              int nchol = ncholpQ[Q];
              int Qm = kminus[Q];
              int Kl = QKToK2[Qm][Ka];
              int Kk = QKToK2[Q][Ka];
              int nl = nopk[Kl];
              int na = nelpk[nd][Ka];
              int nk = nopk[Kk];

              Sp3Tensor_ref Gwal(GKK[0][Ka][Kl].origin(),{nwalk,na,nl});
              Sp3Tensor_ref Gwbk(GKK[0][Ka][Kk].origin(),{nwalk,na,nk});
              Sp3Tensor_ref Lank(sp_pointer(LQKank[nd*nspin*nkpts+Q][Ka].origin()),
                                                 {na,nchol,nk});
              auto bnl_ptr(sp_pointer(LQKank[nd*nspin*nkpts+Qm][Ka].origin()));
              if( Q == Qm ) bnl_ptr = sp_pointer(LQKbnl[nd*nspin*number_of_symmetric_Q+Qmap[Q]-1][Ka].origin());
              Sp3Tensor_ref Lbnl(bnl_ptr,{na,nchol,nl});

              // Twan = sum_l G[w][a][l] L[a][n][l]
              for(int n=0; n<nwalk; ++n) 
                for(int a=0; a<na; ++a)  
                  ma::product(SPComplexType(1.0),Lbnl[a],Gwal[n][a],
                              SPComplexType(1.0),Kl_local[n]);
              for(int n=0; n<nwalk; ++n) 
                for(int a=0; a<na; ++a)  
                  ma::product(SPComplexType(1.0),Lank[a],Gwbk[n][a],
                              SPComplexType(1.0),Kr_local[n]);
            } // if

            if(walker_type==COLLINEAR) {

              { 
                haveKE=true;
                int nchol = ncholpQ[Q];
                int Qm = kminus[Q];
                int Kl = QKToK2[Qm][Ka];
                int Kk = QKToK2[Q][Ka];
                int nl = nopk[Kl];
                int na = nelpk[nd][nkpts+Ka];
                int nk = nopk[Kk];

                Sp3Tensor_ref Gwal(GKK[1][Ka][Kl].origin(),{nwalk,na,nl});
                Sp3Tensor_ref Gwbk(GKK[1][Ka][Kk].origin(),{nwalk,na,nk});
                Sp3Tensor_ref Lank(sp_pointer(LQKank[(nd*nspin+1)*nkpts+Q][Ka].origin()),
                                                 {na,nchol,nk});
                auto bnl_ptr(sp_pointer(LQKank[(nd*nspin+1)*nkpts+Qm][Ka].origin()));
                if( Q == Qm ) bnl_ptr = sp_pointer(LQKbnl[(nd*nspin+1)*number_of_symmetric_Q+Qmap[Q]-1][Ka].origin());
                Sp3Tensor_ref Lbnl(bnl_ptr,{na,nchol,nl});

                // Twan = sum_l G[w][a][l] L[a][n][l]
                for(int n=0; n<nwalk; ++n)
                  for(int a=0; a<na; ++a)  
                    ma::product(SPComplexType(1.0),Lbnl[a],Gwal[n][a],
                                SPComplexType(1.0),Kl_local[n]);
                for(int n=0; n<nwalk; ++n)
                  for(int a=0; a<na; ++a)  
                    ma::product(SPComplexType(1.0),Lank[a],Gwbk[n][a],
                                SPComplexType(1.0),Kr_local[n]);

              } // if
            } // COLLINEAR
          } // Ka
          if(haveKE) {
            int nc0 = Q2vbias[Q]/2; //std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);  
            using ma::axpy;
            for(int n=0; n<nwalk; n++) {
              axpy(SPComplexType(1.0),Kr_local[n].sliced(0,ncholpQ[Q]),
                                        Kr[n].sliced(nc0,nc0+ncholpQ[Q])); 
              axpy(SPComplexType(1.0),Kl_local[n].sliced(0,ncholpQ[Q]),
                                        Kl[n].sliced(nc0,nc0+ncholpQ[Q])); 
            }
          } // to release the lock
          if(haveKE) { 
            fill_n(Kr_local.origin(),Kr_local.num_elements(),SPComplexType(0.0));
            fill_n(Kl_local.origin(),Kl_local.num_elements(),SPComplexType(0.0));
          }  
        } // Q
        nqk=0;  
        RealType scl = (walker_type==CLOSED?2.0:1.0);
        for(int n=0; n<nwalk; ++n) {
          for(int Q=0; Q<nkpts; ++Q) {      // momentum conservation index   
            {
              int nc0 = Q2vbias[Q]/2; //std::accumulate(ncholpQ.begin(),ncholpQ.begin()+Q,0);
              E[n][2] += 0.5*scl*scl*static_cast<ComplexType>(ma::dot(Kl[n]({nc0,nc0+ncholpQ[Q]}),
                                            Kr[n]({nc0,nc0+ncholpQ[Q]})));  
            }
          }
        }
      }
*/
  }

  template<class... Args>
  void fast_energy(Args&&... args)
  {
    APP_ABORT(" Error: fast_energy not implemented in KP3IndexFactorization_batched. \n");
  }

  template<
      class MatA,
      class MatB,
      typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
      typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
      typename = void>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    using BType = typename std::decay<MatB>::type::element;
    using AType = typename std::decay<MatA>::type::element;
    boost::multi::array_ref<AType, 2, decltype(X.origin())> X_(X.origin(), {X.size(), 1});
    boost::multi::array_ref<BType, 2, decltype(v.origin())> v_(v.origin(), {1, v.size()});
    return vHS(X_, v_, a, c);
  }

  template<
      class MatA,
      class MatB,
      typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
      typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>()))
      >
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    int nkpts = nopk.size();
    int nwalk = std::get<1>(X.sizes());
    assert(v.size() == nwalk);
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    assert(X.num_elements() == nwalk * 2 * local_nCV);
    assert(v.num_elements() == nwalk * nmo_tot * nmo_tot);
    SPComplexType one(1.0, 0.0);
    SPComplexType im(0.0, 1.0);
    SPComplexType halfa(0.5 * a, 0.0);
    SPComplexType minusimhalfa(0.0, -0.5 * a);
    SPComplexType imhalfa(0.0, 0.5 * a);

    Static3Tensor vKK({nkpts + number_of_symmetric_Q, nkpts, nwalk * nmo_max * nmo_max},
                      device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    fill_n(vKK.origin(), vKK.num_elements(), SPComplexType(0.0));
    Static4Tensor XQnw({nkpts, 2, nchol_max, nwalk},
                       device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    fill_n(XQnw.origin(), XQnw.num_elements(), SPComplexType(0.0));

    // "rotate" X
    //  XIJ = 0.5*a*(Xn+ -i*Xn-), XJI = 0.5*a*(Xn+ +i*Xn-)
#if defined(MIXED_PRECISION)
    StaticMatrix Xdev(X.extensions(), device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
    copy_n_cast(make_device_ptr(X.origin()), X.num_elements(), Xdev.origin());
#else
    SpMatrix_ref Xdev(make_device_ptr(X.origin()), X.extensions());
#endif
    for (int Q = 0; Q < nkpts; ++Q)
    {
      if (Qmap[Q] < 0)
        continue;
      int nq = Q2vbias[Q];
      auto&& Xp(Xdev.sliced(nq, nq + ncholpQ[Q]));
      auto&& Xm(Xdev.sliced(nq + ncholpQ[Q], nq + 2 * ncholpQ[Q]));
      ma::add(halfa, Xp, minusimhalfa, Xm, XQnw[Q][0].sliced(0, ncholpQ[Q]));
      ma::add(halfa, Xp, imhalfa, Xm, XQnw[Q][1].sliced(0, ncholpQ[Q]));
      nq += 2 * ncholpQ[Q];
    }
    //  then combine Q/(-Q) pieces
    //  X(Q)np = (X(Q)np + X(-Q)nm)
    for (int Q = 0; Q < nkpts; ++Q)
    {
      if (Qmap[Q] == 0)
      {
        int Qm = kminus[Q];
        ma::axpy(SPComplexType(1.0), XQnw[Qm][1], XQnw[Q][0]);
      }
    }
    {
      // assuming contiguous
      ma::scal(c, v);
    }

    int nmo_max2 = nmo_max * nmo_max;
    using ma::gemmBatched;
    std::vector<sp_pointer> Aarray;
    std::vector<sp_pointer> Barray;
    std::vector<sp_pointer> Carray;
    Aarray.reserve(nkpts * nkpts);
    Barray.reserve(nkpts * nkpts);
    Carray.reserve(nkpts * nkpts);
    for (int Q = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Qmap[Q] < 0)
        continue;
      // v[nw][i(in K)][k(in Q(K))] += sum_n LQK[i][k][n] X[Q][0][n][nw]
      if (Q <= kminus[Q])
      {
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK = QKToK2[Q][K];
          Aarray.push_back(sp_pointer(LQKikn[Q][K].origin()));
          Barray.push_back(XQnw[Q][0].origin());
          Carray.push_back(vKK[K][QK].origin());
        }
      }
    }
    // C: v = T(X) * T(Lik) --> F: T(Lik) * T(X) = v
    gemmBatched('T', 'T', nmo_max2, nwalk, nchol_max, SPComplexType(1.0), Aarray.data(), nchol_max, Barray.data(),
                nwalk, SPComplexType(0.0), Carray.data(), nmo_max2, Aarray.size());


    Aarray.clear();
    Barray.clear();
    Carray.clear();
    for (int Q = 0; Q < nkpts; ++Q)
    { // momentum conservation index
      if (Qmap[Q] < 0)
        continue;
      // v[nw][i(in K)][k(in Q(K))] += sum_n LQK[i][k][n] X[Q][0][n][nw]
      if (Q > kminus[Q])
      { // use L(-Q)(ki)*
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK = QKToK2[Q][K];
          Aarray.push_back(sp_pointer(LQKikn[kminus[Q]][QK].origin()));
          Barray.push_back(XQnw[Q][0].origin());
          Carray.push_back(vKK[K][QK].origin());
        }
      }
      else if (Qmap[Q] > 0)
      { // rho(Q)^+ term
        for (int K = 0; K < nkpts; ++K)
        { // K is the index of the kpoint pair of (i,k)
          int QK = QKToK2[Q][K];
          Aarray.push_back(sp_pointer(LQKikn[Q][K].origin()));
          Barray.push_back(XQnw[Q][1].origin());
          Carray.push_back(vKK[nkpts + Qmap[Q] - 1][QK].origin());
        }
      }
    }
    // C: v = T(X) * T(Lik) --> F: T(Lik) * T(X) = v
    gemmBatched('C', 'T', nmo_max2, nwalk, nchol_max, SPComplexType(1.0), Aarray.data(), nchol_max, Barray.data(),
                nwalk, SPComplexType(0.0), Carray.data(), nmo_max2, Aarray.size());


    using vType = typename std::decay<MatB>::type::element;
    boost::multi::array_ref<vType, 3, decltype(make_device_ptr(v.origin()))> v3D(make_device_ptr(v.origin()),
                                                                                 {nwalk, nmo_tot, nmo_tot});
    vKKwij_to_vwKiKj(vKK, v3D);
    // do I need to "rotate" back, can be done if necessary
  }

  template<
      class MatA,
      class MatB,
      typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
      typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatA&>())),
      //             typename = decltype(boost::multi::static_array_cast<ComplexType, pointer>(std::declval<MatB>())),
      typename = void>
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    using BType = typename std::decay<MatB>::type::element;
    using AType = typename std::decay<MatA>::type::element;
    boost::multi::array_ref<BType, 2, decltype(v.origin())> v_(v.origin(), {v.size(), 1});
    boost::multi::array_ref<AType const, 2, decltype(G.origin())> G_(G.origin(), {G.size(), 1});
    return vbias(G_, v_, a, c, k);
  }

  /*
    template<class MatA, class MatB,
             typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality==2)>,
             typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality==1)>,
//             typename = typename std::enable_if_t<(std::is_convertible<typename std::decay<MatA>::type::pointer,pointer>::value)>,  
//             typename = typename std::enable_if_t<(not std::is_convertible<typename std::decay<MatB>::type::pointer,pointer>::value)>,  
              typename = void,  
              typename = void,  
              typename = void,  
              typename = void,  
              typename = void  
            >
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int nd=0) {
    }
*/

  template<
      class MatA,
      class MatB,
      typename = std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
      typename = std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>
      //             typename = std::enable_if_t<(std::is_convertible<typename std::decay<MatA>::type::element_ptr,pointer>::value)>,
      //             typename = std::enable_if_t<(std::is_convertible<typename std::decay<MatB>::type::element_ptr,pointer>::value)>
      >
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int nd = 0)
  {
    using ma::gemmBatched;

    int nkpts = nopk.size();
    assert(nd >= 0 && nd < nelpk.size());
    int nwalk = std::get<1>(G.sizes());
    assert(std::get<0>(v.sizes()) == 2 * local_nCV);
    assert(std::get<1>(v.sizes()) == nwalk);
    int nspin     = (walker_type == COLLINEAR ? 2 : 1);
    int npol      = (walker_type == NONCOLLINEAR ? 2 : 1);
    int nmo_tot   = std::accumulate(nopk.begin(), nopk.end(), 0);
    int nmo_max   = *std::max_element(nopk.begin(), nopk.end());
    int nocca_tot = std::accumulate(nelpk[nd].begin(), nelpk[nd].begin() + nkpts, 0);
    int nocca_max = *std::max_element(nelpk[nd].begin(), nelpk[nd].begin() + nkpts);
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());
    int noccb_tot = 0;
    if (walker_type == COLLINEAR)
    {
      noccb_tot = std::accumulate(nelpk[nd].begin() + nkpts, nelpk[nd].begin() + 2 * nkpts, 0);
    }
    RealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
    SPComplexType one(1.0, 0.0);
    SPComplexType halfa(0.5 * a * scl, 0.0);
    SPComplexType minusimhalfa(0.0, -0.5 * a * scl);
    SPComplexType imhalfa(0.0, 0.5 * a * scl);

    assert(G.num_elements() == nwalk * (nocca_tot + noccb_tot) * npol * nmo_tot);
    // MAM: use reshape when available, then no need to deal with types
    using GType = typename std::decay<MatA>::type::element;
    boost::multi::array_ref<GType const, 3, decltype(make_device_ptr(G.origin()))> G3Da(make_device_ptr(G.origin()),
                                                                                        {nocca_tot * npol, nmo_tot,
                                                                                         nwalk});
    boost::multi::array_ref<GType const, 3, decltype(make_device_ptr(G.origin()))> G3Db(make_device_ptr(G.origin()) +
                                                                                            G3Da.num_elements() *
                                                                                                (nspin - 1),
                                                                                        {noccb_tot, nmo_tot, nwalk});

    // assuming contiguous
    ma::scal(c, v);

    for (int spin = 0; spin < nspin; spin++)
    {
      size_t cnt(0);
      Static3Tensor v1({nkpts + number_of_symmetric_Q, nchol_max, nwalk},
                       device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      Static3Tensor GQ({nkpts, nkpts * nocc_max * npol * nmo_max, nwalk},
                       device_buffer_manager.get_generator().template get_allocator<SPComplexType>());
      fill_n(v1.origin(), v1.num_elements(), SPComplexType(0.0));
      fill_n(GQ.origin(), GQ.num_elements(), SPComplexType(0.0));

      if (spin == 0)
        GKaKjw_to_GQKajw(G3Da, GQ, nelpk[nd], dev_nelpk[nd], dev_a0pk[nd]);
      else
        GKaKjw_to_GQKajw(G3Db, GQ, nelpk[nd].sliced(nkpts, 2 * nkpts), dev_nelpk[nd].sliced(nkpts, 2 * nkpts),
                         dev_a0pk[nd].sliced(nkpts, 2 * nkpts));

      // can use productStridedBatched if LQKakn is changed to a 3Tensor array
      int Kak = nkpts * nocc_max * npol * nmo_max;
      std::vector<sp_pointer> Aarray;
      std::vector<sp_pointer> Barray;
      std::vector<sp_pointer> Carray;
      Aarray.reserve(nkpts + number_of_symmetric_Q);
      Barray.reserve(nkpts + number_of_symmetric_Q);
      Carray.reserve(nkpts + number_of_symmetric_Q);
      for (int Q = 0; Q < nkpts; ++Q)
      { // momentum conservation index
        if (Qmap[Q] < 0)
          continue;
        // v_[Q][n][w] = sum_Kak LQ[Kak][n]*G[Q][Kak][w]
        //             F: -->   G[Kak][w] * LQ[Kak][n]
        Aarray.push_back(GQ[Q].origin());
        Barray.push_back(sp_pointer(LQKakn[nd * nspin * nkpts + spin * nkpts + Q].origin()));
        Carray.push_back(v1[Q].origin());
        if (Qmap[Q] > 0)
        {
          Aarray.push_back(GQ[Q].origin());
          Barray.push_back(sp_pointer(
              LQKbln[nd * nspin * number_of_symmetric_Q + spin * number_of_symmetric_Q + Qmap[Q] - 1].origin()));
          Carray.push_back(v1[nkpts + Qmap[Q] - 1].origin());
        }
      }
      gemmBatched('N', 'T', nwalk, nchol_max, Kak, SPComplexType(1.0), Aarray.data(), nwalk, Barray.data(), nchol_max,
                  SPComplexType(0.0), Carray.data(), nwalk, Aarray.size());
      // optimize later, right now it adds contributions from Q's not assigned
      vbias_from_v1(halfa, v1, v);
    }
  }

  template<class Mat, class MatB>
  void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
  {
    APP_ABORT(" Error: generalizedFockMatrix not implemented for this hamiltonian.\n");
  }

  bool distribution_over_cholesky_vectors() const { return true; }
  int number_of_ke_vectors() const { return local_nCV; }
  int local_number_of_cholesky_vectors() const { return 2 * local_nCV; }
  int global_number_of_cholesky_vectors() const { return global_nCV; }
  int global_origin_cholesky_vector() const { return global_origin; }

  // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
  bool transposed_G_for_vbias() const { return false; }
  bool transposed_G_for_E() const { return false; }
  // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
  bool transposed_vHS() const { return true; }

  bool fast_ph_energy() const { return false; }

  boost::multi::array<ComplexType, 2> getHSPotentials() { return boost::multi::array<ComplexType, 2>{}; }

private:
  int nocc_max;

  afqmc::TaskGroup_& TG;

  Allocator allocator_;
  SpAllocator sp_allocator_;
  DeviceBufferManager device_buffer_manager;

  WALKER_TYPES walker_type;

  int global_nCV;
  int local_nCV;
  int global_origin;

  int default_buffer_size_in_MB;
  int last_nw;

  ValueType E0;

  // bare one body hamiltonian
  mpi3C3Tensor H1;

  // (potentially half rotated) one body hamiltonian
  shmCMatrix haj;
  //std::vector<shmCVector> haj;

  // number of orbitals per k-point
  boost::multi::array<int, 1> nopk;

  // number of cholesky vectors per Q-point
  boost::multi::array<int, 1> ncholpQ;

  // position of (-K) in kp-list for every K
  boost::multi::array<int, 1> kminus;

  // number of electrons per k-point
  // nelpk[ndet][nspin*nkpts]
  //shmIMatrix nelpk;
  boost::multi::array<int, 2> nelpk;

  // maps (Q,K) --> k2
  //shmIMatrix QKToK2;
  boost::multi::array<int, 2> QKToK2;

  //Cholesky Tensor Lik[Q][nk][i][k][n]
  std::vector<shmSpMatrix> LQKikn;

  // half-transformed Cholesky tensor
  std::vector<LQKankMatrix> LQKank;
  const bool needs_copy;

  // half-transformed Cholesky tensor
  std::vector<shmSpMatrix> LQKakn;

  // half-transformed Cholesky tensor
  std::vector<shmSpMatrix> LQKbnl;

  // half-transformed Cholesky tensor
  std::vector<shmSpMatrix> LQKbln;

  // number of Q vectors that satisfy Q==-Q
  int number_of_symmetric_Q;

  // number of Q points assigned to this task
  int number_of_Q_points;

  // Defines behavior over Q vector:
  //   <0: Ignore (handled by another TG)
  //    0: Calculate, without rho^+ contribution
  //   >0: Calculate, with rho^+ contribution. LQKbln data located at Qmap[Q]-1
  stdIVector Qmap;

  // maps Q (only for those with Qmap >=0) to the corresponding sector in vbias
  stdIVector Q2vbias;

  // one-body piece of Hamiltonian factorization
  mpi3C3Tensor vn0;

  int nsampleQ;
  std::vector<RealType> gQ;
  boost::multi::array<int, 2> Qwn;
  std::default_random_engine generator;
  std::discrete_distribution<int> distribution;

  IMatrix KKTransID;
  IVector dev_nopk;
  IVector dev_i0pk;
  IVector dev_kminus;
  IVector dev_ncholpQ;
  IVector dev_Q2vbias;
  IVector dev_Qmap;
  IMatrix dev_nelpk;
  IMatrix dev_a0pk;
  IMatrix dev_QKToK2;

  //    std::vector<std::unique_ptr<shared_mutex>> mutex;

  //    boost::multi::array<ComplexType,3> Qave;
  //    int cntQave=0;
  std::vector<ComplexType> EQ;
  //    std::default_random_engine generator;
  //    std::uniform_real_distribution<RealType> distribution(RealType(0.0),Realtype(1.0));

  template<class MatA, class MatB, class IVec, class IVec2>
  void GKaKjw_to_GKKwaj(MatA const& GKaKj, MatB&& GKKaj, IVec&& nocc, IVec2&& dev_no, IVec2&& dev_a0)
  {
    int npol    = (walker_type == NONCOLLINEAR) ? 2 : 1;
    int nmo_max = *std::max_element(nopk.begin(), nopk.end());
    //      int nocc_max = *std::max_element(nocc.begin(),nocc.end());
    int nmo_tot = std::get<1>(GKaKj.sizes());
    int nwalk   = std::get<2>(GKaKj.sizes());
    int nkpts   = nopk.size();
    assert(GKKaj.num_elements() >= nkpts * nkpts * nwalk * nocc_max * npol * nmo_max);

    using ma::KaKjw_to_KKwaj;
    KaKjw_to_KKwaj(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, dev_nopk.origin(), dev_i0pk.origin(),
                   dev_no.origin(), dev_a0.origin(), GKaKj.origin(), GKKaj.origin());
  }

  template<class MatA, class MatB, class IVec, class IVec2>
  void GKaKjw_to_GQKajw(MatA const& GKaKj, MatB&& GQKaj, IVec&& nocc, IVec2&& dev_no, IVec2&& dev_a0)
  {
    int npol    = (walker_type == NONCOLLINEAR) ? 2 : 1;
    int nmo_max = *std::max_element(nopk.begin(), nopk.end());
    //      int nocc_max = *std::max_element(nocc.begin(),nocc.end());
    int nmo_tot = std::get<1>(GKaKj.sizes());
    int nwalk   = std::get<2>(GKaKj.sizes());
    int nkpts   = nopk.size();
    assert(GQKaj.num_elements() >= nkpts * nkpts * nwalk * nocc_max * npol * nmo_max);

    using ma::KaKjw_to_QKajw;
    KaKjw_to_QKajw(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, dev_nopk.origin(), dev_i0pk.origin(),
                   dev_no.origin(), dev_a0.origin(), dev_QKToK2.origin(), GKaKj.origin(), GQKaj.origin());
  }


  /*
     *   vKiKj({nwalk,nmo_tot,nmo_tot});
     *   vKK({nkpts,nkpts,nwalk*nmo_max*nmo_max} );
     */
  template<class MatA, class MatB>
  void vKKwij_to_vwKiKj(MatA const& vKK, MatB&& vKiKj)
  {
    int nmo_max = *std::max_element(nopk.begin(), nopk.end());
    int nwalk   = std::get<0>(vKiKj.sizes());
    int nmo_tot = std::get<1>(vKiKj.sizes());
    int nkpts   = nopk.size();

    using ma::vKKwij_to_vwKiKj;
    vKKwij_to_vwKiKj(nwalk, nkpts, nmo_max, nmo_tot, KKTransID.origin(), dev_nopk.origin(), dev_i0pk.origin(),
                     vKK.origin(), vKiKj.origin());
  }

  template<class MatA, class MatB>
  void vbias_from_v1(ComplexType a, MatA const& v1, MatB&& vbias)
  {
    using BType   = typename std::decay<MatB>::type::element;
    int nwalk     = std::get<1>(vbias.sizes());
    int nkpts     = nopk.size();
    int nchol_max = *std::max_element(ncholpQ.begin(), ncholpQ.end());

    using ma::vbias_from_v1;
    // using make_device_ptr(vbias.origin()) to catch errors here
    vbias_from_v1(nwalk, nkpts, nchol_max, dev_Qmap.origin(), dev_kminus.origin(), dev_ncholpQ.origin(),
                  dev_Q2vbias.origin(), static_cast<BType>(a), v1.origin(),
                  to_address(make_device_ptr(vbias.origin())));
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
