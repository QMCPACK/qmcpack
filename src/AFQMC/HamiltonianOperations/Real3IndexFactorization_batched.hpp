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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_BATCHED_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_BATCHED_HPP

#include <vector>
#include <type_traits>
#include <random>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/tensor_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// Custom implementation for real build
class Real3IndexFactorization_batched
{
  // allocators
  using Allocator           = device_allocator<ComplexType>;
  using SpAllocator         = device_allocator<SPComplexType>;
  using Allocator_shared    = node_allocator<ComplexType>;
  using SpAllocator_shared  = node_allocator<SPComplexType>;
  using SpRAllocator_shared = node_allocator<SPRealType>;

  // type defs
  using pointer            = typename Allocator::pointer;
  using sp_pointer         = typename SpAllocator::pointer;
  using pointer_shared     = typename Allocator_shared::pointer;
  using sp_pointer_shared  = typename SpAllocator_shared::pointer;
  using sp_rpointer_shared = typename SpRAllocator_shared::pointer;

  using CVector        = ComplexVector<Allocator>;
  using SpVector       = SPComplexVector<SpAllocator>;
  using SpCMatrix      = SPComplexMatrix<SpAllocator>;
  using CVector_ref    = ComplexVector_ref<pointer>;
  using SpCVector_ref  = SPComplexVector_ref<sp_pointer>;
  using CMatrix_ref    = ComplexMatrix_ref<pointer>;
  using SpCMatrix_ref  = SPComplexMatrix_ref<sp_pointer>;
  using SpC4Tensor_ref = boost::multi::array_ref<ComplexType, 4, sp_pointer>;

  using shmCMatrix    = ComplexMatrix<Allocator_shared>;
  using shmSpC3Tensor = SPComplex3Tensor<SpAllocator_shared>;
  using shmSpCMatrix  = SPComplexMatrix<SpAllocator_shared>;
  using shmSpRMatrix  = boost::multi::array<SPRealType, 2, SpRAllocator_shared>;

  using mpi3RMatrix = boost::multi::array<RealType, 2, shared_allocator<RealType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;

public:
  static const HamiltonianTypes HamOpType = RealDenseFactorized;
  HamiltonianTypes getHamType() const { return HamOpType; }

  template<class shmCMatrix_, class shmSpRMatrix_, class shmSpCMatrix_, class shmSpC3Tensor_>
  Real3IndexFactorization_batched(WALKER_TYPES type,
                                  mpi3RMatrix&& hij_,
                                  shmCMatrix_&& haj_,
                                  shmSpRMatrix_&& vik,
                                  shmSpCMatrix_&& vak,
                                  std::vector<shmSpC3Tensor_>&& vank,
                                  mpi3CMatrix&& vn0_,
                                  ValueType e0_,
                                  Allocator const& alloc_,
                                  int cv0,
                                  int gncv)
      : allocator_(alloc_),
        sp_allocator_(alloc_),
        walker_type(type),
        global_origin(cv0),
        global_nCV(gncv),
        local_nCV(0),
        E0(e0_),
        hij(std::move(hij_)),
        haj(std::move(haj_)),
        Likn(std::move(vik)),
        Lank(std::move(move_vector<shmSpC3Tensor>(std::move(vank)))),
        Lakn(std::move(vak)),
        vn0(std::move(vn0_)),
        TBuff(iextensions<1u>{1}, sp_allocator_)
  {
    local_nCV = Likn.size(1);
    size_t lank(0);
    for (auto& v : Lank)
      lank += v.num_elements();
    app_log() << "****************************************************************** \n"
              << "  Static memory usage by Real3IndexFactorization_batched (node 0 in MB) \n"
              << "  Likn: " << Likn.num_elements() * sizeof(SPRealType) / 1024.0 / 1024.0 << " \n"
              << "  Lakn: " << Lakn.num_elements() * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n"
              << "  Lank: " << lank * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n";
    memory_report();
  }

  ~Real3IndexFactorization_batched() {}

  Real3IndexFactorization_batched(const Real3IndexFactorization_batched& other) = delete;
  Real3IndexFactorization_batched& operator=(const Real3IndexFactorization_batched& other) = delete;
  Real3IndexFactorization_batched(Real3IndexFactorization_batched&& other)                 = default;
  Real3IndexFactorization_batched& operator=(Real3IndexFactorization_batched&& other) = delete;

  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(TaskGroup_& TG,
                                                                 boost::multi::array<ComplexType, 1> const& vMF)
  {
    int NMO = hij.size(0);
    // in non-collinear case with SO, keep SO matrix here and add it
    // for now, stay collinear

    CVector vMF_(vMF);
    CVector P1D(iextensions<1u>{NMO * NMO});
    fill_n(P1D.origin(), P1D.num_elements(), ComplexType(0));
    vHS(vMF_, P1D);
    if (TG.TG().size() > 1)
      TG.TG().all_reduce_in_place_n(to_address(P1D.origin()), P1D.num_elements(), std::plus<>());

    boost::multi::array<ComplexType, 2> P1({NMO, NMO});
    copy_n(P1D.origin(), NMO * NMO, P1.origin());

    using ma::conj;

    for (int i = 0; i < NMO; i++)
    {
      P1[i][i] += hij[i][i] + vn0[i][i];
      for (int j = i + 1; j < NMO; j++)
      {
        P1[i][j] += hij[i][j] + vn0[i][j];
        P1[j][i] += hij[j][i] + vn0[j][i];
        // This is really cutoff dependent!!!
        if (std::abs(P1[i][j] - ma::conj(P1[j][i])) > 1e-6)
        {
          app_error() << " WARNING in getOneBodyPropagatorMatrix. P1 is not hermitian. \n";
          app_error() << i << " " << j << " " << P1[i][j] << " " << P1[j][i] << " " << hij[i][j] << " " << hij[j][i]
                      << " " << vn0[i][j] << " " << vn0[j][i] << std::endl;
          //APP_ABORT("Error in getOneBodyPropagatorMatrix. P1 is not hermitian. \n");
        }
        P1[i][j] = 0.5 * (P1[i][j] + ma::conj(P1[j][i]));
        P1[j][i] = ma::conj(P1[i][j]);
      }
    }

    return P1;
  }

  template<class Mat, class MatB>
  void energy(Mat&& E, MatB const& G, int k, bool addH1 = true, bool addEJ = true, bool addEXX = true)
  {
    MatB* Kr(nullptr);
    MatB* Kl(nullptr);
    energy(E, G, k, Kl, Kr, addH1, addEJ, addEXX);
  }

  // KEleft and KEright must be in shared memory for this to work correctly
  template<class Mat, class MatB, class MatC, class MatD>
  void energy(Mat&& E,
              MatB const& Gc,
              int nd,
              MatC* KEleft,
              MatD* KEright,
              bool addH1  = true,
              bool addEJ  = true,
              bool addEXX = true)
  {
    assert(E.size(1) >= 3);
    assert(nd >= 0);
    assert(nd < haj.size());
    if (walker_type == COLLINEAR)
      assert(2 * nd + 1 < Lank.size());
    else
      assert(nd < Lank.size());

    int nwalk = Gc.size(0);
    int nspin = (walker_type == COLLINEAR ? 2 : 1);
    int NMO   = hij.size(0);
    int nel[2];
    nel[0] = Lank[nspin * nd].size(0);
    nel[1] = ((nspin == 2) ? Lank[nspin * nd + 1].size(0) : 0);
    assert(Lank[nspin * nd].size(1) == local_nCV);
    assert(Lank[nspin * nd].size(2) == NMO);
    if (nspin == 2)
    {
      assert(Lank[nspin * nd + 1].size(1) == local_nCV);
      assert(Lank[nspin * nd + 1].size(2) == NMO);
    }
    assert(Gc.num_elements() == nwalk * (nel[0] + nel[1]) * NMO);

    int getKr = KEright != nullptr;
    int getKl = KEleft != nullptr;
    if (E.size(0) != nwalk || E.size(1) < 3)
      APP_ABORT(" Error in AFQMC/HamiltonianOperations/Real3IndexFactorization_batched::energy(...). Incorrect matrix "
                "dimensions \n");

    // T[nwalk][nup][nup][local_nCV] + D[nwalk][nwalk][local_nCV]
    size_t mem_needs(0);
    size_t cnt(0);
    if (addEJ)
    {
#if defined(MIXED_PRECISION)
      mem_needs += nwalk * local_nCV;
#else
      if (not getKl)
        mem_needs += nwalk * local_nCV;
#endif
    }
    if (addEXX)
    {
      mem_needs += nwalk * nel[0] * nel[0] * local_nCV;
#if defined(MIXED_PRECISION)
      mem_needs += nwalk * nel[0] * NMO;
#else
      if (nspin == 2)
        mem_needs += nwalk * nel[0] * NMO;
#endif
    }
    set_buffer(mem_needs);

    // messy
    sp_pointer Klptr(nullptr);
    long Knr = 0, Knc = 0;
    if (addEJ)
    {
      Knr = nwalk;
      Knc = local_nCV;
      if (getKr)
      {
        assert(KEright->size(0) == nwalk && KEright->size(1) == local_nCV);
        assert(KEright->stride(0) == KEright->size(1));
      }
#if defined(MIXED_PRECISION)
      if (getKl)
      {
        assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
        assert(KEleft->stride(0) == KEleft->size(1));
      }
#else
      if (getKl)
      {
        assert(KEleft->size(0) == nwalk && KEleft->size(1) == local_nCV);
        assert(KEleft->stride(0) == KEleft->size(1));
        Klptr = make_device_ptr(KEleft->origin());
      }
      else
#endif
      {
        Klptr = TBuff.origin() + cnt;
        cnt += Knr * Knc;
      }
      fill_n(Klptr, Knr * Knc, SPComplexType(0.0));
    }
    else if (getKr or getKl)
    {
      APP_ABORT(" Error: Kr and/or Kl can only be calculated with addEJ=true.\n");
    }
    SpCMatrix_ref Kl(Klptr, {long(Knr), long(Knc)});

    for (int n = 0; n < nwalk; n++)
      std::fill_n(E[n].origin(), 3, ComplexType(0.));


    // one-body contribution
    // haj[ndet][nocc*nmo]
    // not parallelized for now, since it would require customization of Wfn
    if (addH1)
    {
      CVector_ref haj_ref(make_device_ptr(haj[nd].origin()), iextensions<1u>{haj[nd].num_elements()});
      ma::product(ComplexType(1.), Gc, haj_ref, ComplexType(1.), E(E.extension(0), 0));
      for (int i = 0; i < nwalk; i++)
        E[i][0] += E0;
    }

    // move calculation of H1 here
    // NOTE: For CLOSED/NONCOLLINEAR, can do all walkers simultaneously to improve perf. of GEMM
    //       Not sure how to do it for COLLINEAR.
    if (addEXX)
    {
      SPRealType scl = (walker_type == CLOSED ? 2.0 : 1.0);

      for (int ispin = 0, is0 = 0; ispin < nspin; ispin++)
      {
        size_t cnt_(cnt);
        sp_pointer ptr(nullptr);
#if defined(MIXED_PRECISION)
        ptr = TBuff.origin() + cnt_;
        cnt_ += nwalk * nel[ispin] * NMO;
        for (int n = 0; n < nwalk; ++n)
        {
          copy_n_cast(make_device_ptr(Gc[n].origin()) + is0, nel[ispin] * NMO, ptr + n * nel[ispin] * NMO);
        }
#else
        if (nspin == 1)
        {
          ptr = make_device_ptr(Gc.origin());
        }
        else
        {
          ptr = TBuff.origin() + cnt_;
          cnt_ += nwalk * nel[ispin] * NMO;
          for (int n = 0; n < nwalk; ++n)
          {
            using std::copy_n;
            copy_n(make_device_ptr(Gc[n].origin()) + is0, nel[ispin] * NMO, ptr + n * nel[ispin] * NMO);
          }
        }
#endif

        SpCMatrix_ref GF(ptr, {nwalk * nel[ispin], NMO});
        SpCMatrix_ref Lan(make_device_ptr(Lank[nd * nspin + ispin].origin()), {nel[ispin] * local_nCV, NMO});
        SpCMatrix_ref Twban(TBuff.origin() + cnt_, {nwalk * nel[ispin], nel[ispin] * local_nCV});
        SpC4Tensor_ref T4Dwban(Twban.origin(), {nwalk, nel[ispin], nel[ispin], local_nCV});

        ma::product(GF, ma::T(Lan), Twban);

        using ma::dot_wabn;
        dot_wabn(nwalk, nel[ispin], local_nCV, SPComplexType(-0.5 * scl), Twban.origin(), to_address(E[0].origin()) + 1,
                 E.stride(0));
        //          for(int n=0, an=0; n<nwalk; ++n) {
        //            ComplexType E_(0.0);
        //            for(int a=0; a<nel[ispin]; ++a, an++) {
        //              for(int b=0; b<nel[ispin]; ++b)
        //                E_ += static_cast<ComplexType>(ma::dot(T4Dwban[n][a][b],T4Dwban[n][b][a]));
        //            }
        //            E[n][1] -= 0.5*scl*E_;
        //          }

        if (addEJ)
        {
          //            for(int n=0; n<nwalk; ++n) {
          //              for(int a=0; a<nel[ispin]; ++a)
          //                ma::axpy(SPComplexType(1.0),T4Dwban[n][a][a],Kl[n]);
          //            }
          using ma::Tab_to_Kl;
          Tab_to_Kl(nwalk, nel[ispin], local_nCV, Twban.origin(), Kl.origin());
        }
        is0 += nel[ispin] * NMO;

      } // if
    }

    if (addEJ)
    {
      if (not addEXX)
      {
        // calculate Kr
        APP_ABORT(" Error: Finish addEJ and not addEXX");
      }
      SPRealType scl = (walker_type == CLOSED ? 2.0 : 1.0);
      for (int n = 0; n < nwalk; ++n)
        E[n][2] += 0.5 * static_cast<ComplexType>(scl * scl * ma::dot(Kl[n], Kl[n]));
#if defined(MIXED_PRECISION)
      if (getKl)
        copy_n_cast(Klptr, KEleft->num_elements(), make_device_ptr(KEleft->origin()));
#endif
      if (getKr)
        copy_n_cast(Klptr, KEright->num_elements(), make_device_ptr(KEright->origin()));
    }
  }

  template<class... Args>
  void fast_energy(Args&&... args)
  {
    APP_ABORT(" Error: fast_energy not implemented in Real3IndexFactorization_batched. \n");
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    assert(Likn.size(1) == X.size(0));
    assert(Likn.size(0) == v.size(0));
#if defined(MIXED_PRECISION)
    size_t mem_needs = X.num_elements() + v.num_elements();
    set_buffer(mem_needs);
    SpCVector_ref vsp(TBuff.origin(), v.extensions());
    SpCVector_ref Xsp(vsp.origin() + vsp.num_elements(), X.extensions());
    copy_n_cast(make_device_ptr(X.origin()), X.num_elements(), Xsp.origin());
    ma::product(SPValueType(a), Likn, Xsp, SPValueType(c), vsp);
    copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
    ma::product(SPValueType(a), Likn, X, SPValueType(c), v);
#endif
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    assert(Likn.size(1) == X.size(0));
    assert(Likn.size(0) == v.size(0));
    assert(X.size(1) == v.size(1));
#if defined(MIXED_PRECISION)
    size_t mem_needs = X.num_elements() + v.num_elements();
    set_buffer(mem_needs);
    SpCMatrix_ref vsp(TBuff.origin(), v.extensions());
    SpCMatrix_ref Xsp(vsp.origin() + vsp.num_elements(), X.extensions());
    copy_n_cast(make_device_ptr(X.origin()), X.num_elements(), Xsp.origin());
    ma::product(SPValueType(a), Likn, Xsp, SPValueType(c), vsp);
    copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
    ma::product(SPValueType(a), Likn, X, SPValueType(c), v);
#endif
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    if (haj.size(0) == 1)
    {
      assert(Lakn.size(0) == G.size(0));
      assert(Lakn.size(1) == v.size(0));

#if defined(MIXED_PRECISION)
      size_t mem_needs = G.num_elements() + v.num_elements();
      set_buffer(mem_needs);
      SpCVector_ref vsp(TBuff.origin(), v.extensions());
      SpCVector_ref Gsp(vsp.origin() + vsp.num_elements(), G.extensions());
      copy_n_cast(make_device_ptr(G.origin()), G.num_elements(), Gsp.origin());
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPComplexType(a), ma::T(Lakn), Gsp, SPComplexType(c), vsp);
      copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPComplexType(a), ma::T(Lakn), G, SPComplexType(c), v);
#endif
    }
    else
    {
      // multideterminant is not half-rotated, so use Likn
      assert(Likn.size(0) == G.size(0));
      assert(Likn.size(1) == v.size(0));

#if defined(MIXED_PRECISION)
      size_t mem_needs = G.num_elements() + v.num_elements();
      set_buffer(mem_needs);
      SpCVector_ref vsp(TBuff.origin(), v.extensions());
      SpCVector_ref Gsp(vsp.origin() + vsp.num_elements(), G.extensions());
      copy_n_cast(make_device_ptr(G.origin()), G.num_elements(), Gsp.origin());
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPValueType(a), ma::T(Likn), Gsp, SPValueType(c), vsp);
      copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPValueType(a), ma::T(Likn), G, SPValueType(c), v);
#endif
    }
  }

  // v(n,w) = sum_ak L(ak,n) G(w,ak)
  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    if (haj.size(0) == 1)
    {
      assert(Lakn.size(0) == G.size(0));
      assert(Lakn.size(1) == v.size(0));
      assert(G.size(1) == v.size(1));

#if defined(MIXED_PRECISION)
      size_t mem_needs = G.num_elements() + v.num_elements();
      set_buffer(mem_needs);
      SpCMatrix_ref vsp(TBuff.origin(), v.extensions());
      SpCMatrix_ref Gsp(vsp.origin() + vsp.num_elements(), G.extensions());
      copy_n_cast(make_device_ptr(G.origin()), G.num_elements(), Gsp.origin());
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPComplexType(a), ma::T(Lakn), Gsp, SPComplexType(c), vsp);
      copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPComplexType(a), ma::T(Lakn), G, SPComplexType(c), v);
#endif
    }
    else
    {
      // multideterminant is not half-rotated, so use Likn
      assert(Likn.size(0) == G.size(0));
      assert(Likn.size(1) == v.size(0));
      assert(G.size(1) == v.size(1));

#if defined(MIXED_PRECISION)
      size_t mem_needs = G.num_elements() + v.num_elements();
      set_buffer(mem_needs);
      SpCMatrix_ref vsp(TBuff.origin(), v.extensions());
      SpCMatrix_ref Gsp(vsp.origin() + vsp.num_elements(), G.extensions());
      copy_n_cast(make_device_ptr(G.origin()), G.num_elements(), Gsp.origin());
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPValueType(a), ma::T(Likn), Gsp, SPValueType(c), vsp);
      copy_n_cast(vsp.origin(), vsp.num_elements(), make_device_ptr(v.origin()));
#else
      if (walker_type == CLOSED)
        a *= 2.0;
      ma::product(SPValueType(a), ma::T(Likn), G, SPValueType(c), v);
#endif
    }
  }

  template<class Mat, class MatB>
  void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
  {
    APP_ABORT(" Error: generalizedFockMatrix not implemented for this hamiltonian.\n");
  }

  bool distribution_over_cholesky_vectors() const { return true; }
  int number_of_ke_vectors() const { return local_nCV; }
  int local_number_of_cholesky_vectors() const { return local_nCV; }
  int global_number_of_cholesky_vectors() const { return global_nCV; }
  int global_origin_cholesky_vector() const { return global_origin; }

  // transpose=true means G[nwalk][ik], false means G[ik][nwalk]
  bool transposed_G_for_vbias() const { return false; }
  bool transposed_G_for_E() const { return true; }
  // transpose=true means vHS[nwalk][ik], false means vHS[ik][nwalk]
  bool transposed_vHS() const { return false; }

  bool fast_ph_energy() const { return false; }

  boost::multi::array<ComplexType, 2> getHSPotentials() { return boost::multi::array<ComplexType, 2>{}; }

private:
  Allocator allocator_;
  SpAllocator sp_allocator_;

  WALKER_TYPES walker_type;

  int global_origin;
  int global_nCV;
  int local_nCV;

  ValueType E0;

  // bare one body hamiltonian
  mpi3RMatrix hij;

  // (potentially half rotated) one body hamiltonian
  shmCMatrix haj;

  //Cholesky Tensor Lik[i][k][n]
  shmSpRMatrix Likn;

  // permuted half-transformed Cholesky tensor
  // Lank[ 2*idet + ispin ]
  std::vector<shmSpC3Tensor> Lank;

  // half-transformed Cholesky tensor
  // only used in single determinant case, haj.size(0)==1.
  shmSpCMatrix Lakn;

  // one-body piece of Hamiltonian factorization
  mpi3CMatrix vn0;

  // shared buffer space
  // using matrix since there are issues with vectors
  SpVector TBuff;

  myTimer Timer;

  void set_buffer(size_t N)
  {
    if (TBuff.num_elements() < N)
    {
      app_log() << " Resizing buffer space in Real3IndexFactorization_batched to "
                << N * sizeof(SPComplexType) / 1024.0 / 1024.0 << " MBs. \n";
      {
        TBuff = std::move(SpVector(iextensions<1u>{N}));
      }
      memory_report();
      using std::fill_n;
      fill_n(TBuff.origin(), N, SPComplexType(0.0));
    }
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
