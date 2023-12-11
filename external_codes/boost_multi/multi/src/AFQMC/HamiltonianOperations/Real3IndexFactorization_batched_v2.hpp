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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_BATCHED_V2_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_REAL3INDEXFACTORIZATION_BATCHED_V2_HPP

#include <vector>
#include <type_traits>
#include <random>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Numerics/tensor_operations.hpp"
#include "AFQMC/Utilities/myTimer.h"

namespace qmcplusplus
{
namespace afqmc
{
// Custom implementation for real build
class Real3IndexFactorization_batched_v2
{
  // allocators
  using Allocator           = device_allocator<ComplexType>;
  using SpAllocator         = device_allocator<SPComplexType>;
  using SpRAllocator        = device_allocator<SPRealType>;
  using Allocator_shared    = node_allocator<ComplexType>;
  using SpAllocator_shared  = node_allocator<SPComplexType>;
  using SpRAllocator_shared = node_allocator<SPRealType>;

  // type defs
  using pointer            = typename Allocator::pointer;
  using sp_pointer         = typename SpAllocator::pointer;
  using const_sp_pointer   = typename SpAllocator::const_pointer;
  using sp_rpointer        = typename SpRAllocator::pointer;
  using pointer_shared     = typename Allocator_shared::pointer;
  using sp_pointer_shared  = typename SpAllocator_shared::pointer;
  using sp_rpointer_shared = typename SpRAllocator_shared::pointer;

  using device_alloc_type  = DeviceBufferManager::template allocator_t<SPComplexType>;
  using device_alloc_Rtype = DeviceBufferManager::template allocator_t<SPRealType>;

  using CVector        = ComplexVector<Allocator>;
  using SpVector       = SPComplexVector<SpAllocator>;
  using SpCMatrix      = SPComplexMatrix<SpAllocator>;
  using CVector_ref    = ComplexVector_ref<pointer>;
  using SpCVector_ref  = SPComplexVector_ref<sp_pointer>;
  using CMatrix_ref    = ComplexMatrix_ref<pointer>;
  using SpCMatrix_ref  = SPComplexMatrix_ref<sp_pointer>;
  using SpRMatrix_ref  = SPComplexMatrix_ref<sp_rpointer>;
  using SpCTensor_ref  = boost::multi::array_ref<SPComplexType, 3, sp_pointer>;
  using SpC4Tensor_ref = boost::multi::array_ref<SPComplexType, 4, sp_pointer>;
  using C4Tensor_ref   = boost::multi::array_ref<ComplexType, 4, pointer>;

  using StaticVector   = boost::multi::static_array<SPComplexType, 1, device_alloc_type>;
  using StaticMatrix   = boost::multi::static_array<SPComplexType, 2, device_alloc_type>;
  using Static3Tensor  = boost::multi::static_array<SPComplexType, 3, device_alloc_type>;
  using Static4Tensor  = boost::multi::static_array<SPComplexType, 4, device_alloc_type>;
  using StaticRVector  = boost::multi::static_array<SPRealType, 1, device_alloc_Rtype>;
  using StaticRMatrix  = boost::multi::static_array<SPRealType, 2, device_alloc_Rtype>;
  using Static3RTensor = boost::multi::static_array<SPRealType, 3, device_alloc_Rtype>;
  using Static4RTensor = boost::multi::static_array<SPRealType, 4, device_alloc_Rtype>;

  using shmCMatrix    = ComplexMatrix<Allocator_shared>;
  using shmSpC3Tensor = SPComplex3Tensor<SpAllocator_shared>;
  using shmSpCMatrix  = SPComplexMatrix<SpAllocator_shared>;
  using shmSpRMatrix  = boost::multi::array<SPRealType, 2, SpRAllocator_shared>;

  using mpi3RMatrix = boost::multi::array<RealType, 2, shared_allocator<RealType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;

public:
  static const HamiltonianTypes HamOpType = RealDenseFactorized;
  HamiltonianTypes getHamType() const { return HamOpType; }

  template<class shmCMatrix_, class shmSpRMatrix_, class shmSpC3Tensor_>
  Real3IndexFactorization_batched_v2(WALKER_TYPES type,
                                     TaskGroup_& TG,
                                     mpi3RMatrix&& hij_,
                                     shmCMatrix_&& haj_,
                                     shmSpRMatrix_&& vik,
                                     std::vector<shmSpC3Tensor_>&& vnak,
                                     mpi3CMatrix&& vn0_,
                                     ValueType e0_,
                                     Allocator const& alloc_,
                                     int cv0,
                                     int gncv,
                                     long maxMem = 2000)
      : allocator_(alloc_),
        sp_allocator_(alloc_),
        buffer_manager(),
        walker_type(type),
        max_memory_MB(maxMem),
        global_origin(cv0),
        global_nCV(gncv),
        local_nCV(0),
        E0(e0_),
        hij(std::move(hij_)),
        hij_dev(hij.extensions(), make_node_allocator<ComplexType>(TG)),
        haj(std::move(haj_)),
        Likn(std::move(vik)),
        Lnak(std::move(move_vector<shmSpC3Tensor>(std::move(vnak)))),
        vn0(std::move(vn0_))
  {
    local_nCV = std::get<1>(Likn.sizes());
    size_t lnak(0);
    for (auto& v : Lnak)
      lnak += v.num_elements();
    for (int i = 0; i < std::get<0>(hij.sizes()); i++)
    {
      for (int j = 0; j < std::get<1>(hij.sizes()); j++)
      {
        hij_dev[i][j] = ComplexType(hij[i][j]);
      }
    }
    app_log() << "****************************************************************** \n"
              << "  Static memory usage by Real3IndexFactorization_batched_v2 (node 0 in MB) \n"
              << "  Likn: " << Likn.num_elements() * sizeof(SPRealType) / 1024.0 / 1024.0 << " \n"
              << "  Lnak: " << lnak * sizeof(SPComplexType) / 1024.0 / 1024.0 << " \n"
              << "  Buffer memory limited to (not yet allocated) :" << max_memory_MB << " MB. \n";
    memory_report();
  }

  ~Real3IndexFactorization_batched_v2() {}

  Real3IndexFactorization_batched_v2(const Real3IndexFactorization_batched_v2& other) = delete;
  Real3IndexFactorization_batched_v2& operator=(const Real3IndexFactorization_batched_v2& other) = delete;
  Real3IndexFactorization_batched_v2(Real3IndexFactorization_batched_v2&& other)                 = default;
  Real3IndexFactorization_batched_v2& operator=(Real3IndexFactorization_batched_v2&& other) = delete;

  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(TaskGroup_& TG,
                                                                 boost::multi::array<ComplexType, 1> const& vMF)
  {
    int NMO = hij.size();
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
    assert(std::get<1>(E.sizes()) >= 3);
    assert(nd >= 0);
    assert(nd < haj.size());
    if (walker_type == COLLINEAR)
      assert(2 * nd + 1 < Lnak.size());
    else
      assert(nd < Lnak.size());

    int nwalk = Gc.size();
    int nspin = (walker_type == COLLINEAR ? 2 : 1);
    int NMO   = hij.size();
    int nel[2];
    nel[0] = std::get<1>(Lnak[nspin * nd].sizes());
    nel[1] = ((nspin == 2) ? std::get<1>(Lnak[nspin * nd + 1].sizes()) : 0);
    assert(std::get<0>(Lnak[nspin * nd].sizes()) == local_nCV);
    assert(std::get<2>(Lnak[nspin * nd].sizes()) == NMO);
    if (nspin == 2)
    {
      assert(std::get<0>(Lnak[nspin * nd + 1].sizes()) == local_nCV);
      assert(std::get<2>(Lnak[nspin * nd + 1].sizes()) == NMO);
    }
    assert(Gc.num_elements() == nwalk * (nel[0] + nel[1]) * NMO);

    int getKr = KEright != nullptr;
    int getKl = KEleft != nullptr;
    if (std::get<0>(E.sizes()) != nwalk || std::get<1>(E.sizes()) < 3)
      APP_ABORT(" Error in AFQMC/HamiltonianOperations/Real3IndexFactorization_batched_v2::energy(...). Incorrect "
                "matrix dimensions \n");

    int max_nCV = 0;
    if (addEXX)
    {
      long LBytes = max_memory_MB * 1024L * 1024L;
      int Bytes   = int(LBytes / long(nwalk * nel[0] * nel[0] * sizeof(SPComplexType)));
      max_nCV     = std::min(std::max(1, Bytes), local_nCV);
      assert(max_nCV > 1 && max_nCV <= local_nCV);
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
    StaticMatrix Kl({Knr, Knc}, buffer_manager.get_generator().template get_allocator<SPComplexType>());
    fill_n(Kl.origin(), Kl.num_elements(), SPComplexType(0.0));

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

      long mem_needs(0);
#if defined(MIXED_PRECISION)
      mem_needs = nwalk * nel[0] * NMO;
#else
      if (nspin > 1)
        mem_needs = nwalk * std::max(nel[0], nel[1]) * NMO;
#endif
      StaticVector T1(iextensions<1u>{mem_needs},
                      buffer_manager.get_generator().template get_allocator<SPComplexType>());

      for (int ispin = 0, is0 = 0; ispin < nspin; ispin++)
      {
        sp_pointer ptr(nullptr);
#if defined(MIXED_PRECISION)
        ptr = T1.origin();
        for (int n = 0; n < nwalk; ++n)
          copy_n_cast(make_device_ptr(Gc[n].origin()) + is0, nel[ispin] * NMO, ptr + n * nel[ispin] * NMO);
#else
        if (nspin == 1)
        {
          ptr = make_device_ptr(Gc.origin());
        }
        else
        {
          ptr = T1.origin();
          for (int n = 0; n < nwalk; ++n)
          {
            using std::copy_n;
            copy_n(make_device_ptr(Gc[n].origin()) + is0, nel[ispin] * NMO, ptr + n * nel[ispin] * NMO);
          }
        }
#endif

        SpCMatrix_ref GF(ptr, {nwalk * nel[ispin], NMO});

        int nCV = 0;
        while (nCV < local_nCV)
        {
          int nvecs = std::min(local_nCV - nCV, max_nCV);
          SpCMatrix_ref Lna(make_device_ptr(Lnak[nd * nspin + ispin][nCV].origin()), {nvecs * nel[ispin], NMO});
          StaticMatrix Twbna({nwalk * nel[ispin], nvecs * nel[ispin]},
                             buffer_manager.get_generator().template get_allocator<SPComplexType>());
          SpC4Tensor_ref T4Dwbna(Twbna.origin(), {nwalk, nel[ispin], nvecs, nel[ispin]});

          ma::product(GF, ma::T(Lna), Twbna);

          using ma::dot_wanb;
          dot_wanb(nwalk, nel[ispin], nvecs, SPComplexType(-0.5 * scl), Twbna.origin(), to_address(E[0].origin()) + 1,
                   E.stride(0));

          if (addEJ)
          {
            using ma::Tanb_to_Kl;
            Tanb_to_Kl(nwalk, nel[ispin], nvecs, local_nCV, Twbna.origin(), Kl.origin() + nCV);
          }

          nCV += max_nCV;
        }
        is0 += nel[ispin] * NMO;

      } // ispin
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
      if (getKl)
        copy_n_cast(Kl.origin(), KEleft->num_elements(), make_device_ptr(KEleft->origin()));
      if (getKr)
        copy_n_cast(Kl.origin(), KEright->num_elements(), make_device_ptr(KEright->origin()));
    }
  }

  template<class... Args>
  void fast_energy(Args&&... args)
  {
    APP_ABORT(" Error: fast_energy not implemented in Real3IndexFactorization_batched_v2. \n");
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    using BType = typename std::decay<MatB>::type::element;
    using AType = typename std::decay<MatA>::type::element;
    boost::multi::array_ref<BType, 2, decltype(v.origin())> v_(v.origin(), {std::get<0>(v.sizes()), 1});
    boost::multi::array_ref<AType, 2, decltype(X.origin())> X_(X.origin(), {std::get<0>(X.sizes()), 1});
    return vHS(X_, v_, a, c);
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vHS(MatA& X, MatB&& v, double a = 1., double c = 0.)
  {
    using XType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay<MatB>::type::element;
    assert(std::get<1>(Likn.sizes()) == std::get<0>(X.sizes()));
    assert(std::get<0>(Likn.sizes()) == std::get<0>(v.sizes()));
    assert(std::get<1>(X.sizes()) == std::get<1>(v.sizes()));
    // setup buffer space if changing precision in X or v
    size_t vmem(0), Xmem(0);
    if (not std::is_same<XType, SPComplexType>::value)
      Xmem = X.num_elements();
    if (not std::is_same<vType, SPComplexType>::value)
      vmem = v.num_elements();
    StaticVector SPBuff(iextensions<1u>(Xmem + vmem),
                        buffer_manager.get_generator().template get_allocator<SPComplexType>());
    sp_pointer vptr(nullptr);
    const_sp_pointer Xptr(nullptr);
    // setup origin of Gsp and copy_n_cast if necessary
    using qmcplusplus::afqmc::pointer_cast;
    if (std::is_same<XType, SPComplexType>::value)
    {
      Xptr = pointer_cast<SPComplexType const>(make_device_ptr(X.origin()));
    }
    else
    {
      copy_n_cast(make_device_ptr(X.origin()), X.num_elements(), make_device_ptr(SPBuff.origin()));
      Xptr = make_device_ptr(SPBuff.origin());
    }
    // setup origin of vsp and copy_n_cast if necessary
    if (std::is_same<vType, SPComplexType>::value)
    {
      vptr = pointer_cast<SPComplexType>(make_device_ptr(v.origin()));
    }
    else
    {
      vptr = make_device_ptr(SPBuff.origin()) + Xmem;
      if (std::abs(c) > 1e-12)
        copy_n_cast(make_device_ptr(v.origin()), v.num_elements(), vptr);
    }
    // work
    boost::multi::array_cref<SPComplexType const, 2, const_sp_pointer> Xsp(Xptr, X.extensions());
    boost::multi::array_ref<SPComplexType, 2, sp_pointer> vsp(vptr, v.extensions());
    ma::product(SPValueType(a), Likn, Xsp, SPValueType(c), vsp);
    if (not std::is_same<vType, SPComplexType>::value)
    {
      copy_n_cast(make_device_ptr(vsp.origin()), v.num_elements(), make_device_ptr(v.origin()));
    }
  }

  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 1)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 1)>,
           typename = void>
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    using BType = typename std::decay<MatB>::type::element;
    using AType = typename std::decay<MatA>::type::element;
    boost::multi::array_ref<BType, 2, decltype(v.origin())> v_(v.origin(), {v.size(), 1});
    boost::multi::array_ref<AType const, 2, decltype(G.origin())> G_(G.origin(), {G.size(), 1});
    return vbias(G_, v_, a, c, k);
  }

  // v(n,w) = sum_ak L(ak,n) G(w,ak)
  template<class MatA,
           class MatB,
           typename = typename std::enable_if_t<(std::decay<MatA>::type::dimensionality == 2)>,
           typename = typename std::enable_if_t<(std::decay<MatB>::type::dimensionality == 2)>>
  void vbias(const MatA& G, MatB&& v, double a = 1., double c = 0., int k = 0)
  {
    using GType = typename std::decay_t<typename MatA::element>;
    using vType = typename std::decay_t<MatB>::element;
    if (walker_type == CLOSED)
      a *= 2.0;
    // setup buffer space if changing precision in G or v
    size_t vmem(0), Gmem(0);
    if (not std::is_same<GType, SPComplexType>::value)
      Gmem = G.num_elements();
    if (not std::is_same<vType, SPComplexType>::value)
      vmem = v.num_elements();
    StaticVector SPBuff(iextensions<1u>(Gmem + vmem),
                        buffer_manager.get_generator().template get_allocator<SPComplexType>());
    sp_pointer vptr(nullptr);
    const_sp_pointer Gptr(nullptr);
    // setup origin of Gsp and copy_n_cast if necessary
    using qmcplusplus::afqmc::pointer_cast;
    if (std::is_same<GType, SPComplexType>::value)
    {
      Gptr = pointer_cast<SPComplexType const>(make_device_ptr(G.origin()));
    }
    else
    {
      copy_n_cast(make_device_ptr(G.origin()), G.num_elements(), make_device_ptr(SPBuff.origin()));
      Gptr = make_device_ptr(SPBuff.origin());
    }
    // setup origin of vsp and copy_n_cast if necessary
    if (std::is_same<vType, SPComplexType>::value)
    {
      vptr = pointer_cast<SPComplexType>(make_device_ptr(v.origin()));
    }
    else
    {
      vptr = make_device_ptr(SPBuff.origin()) + Gmem;
      if (std::abs(c) > 1e-12)
        copy_n_cast(make_device_ptr(v.origin()), v.num_elements(), vptr);
    }
    // setup array references
    boost::multi::array_cref<SPComplexType const, 2, const_sp_pointer> Gsp(Gptr, G.extensions());
    boost::multi::array_ref<SPComplexType, 2, sp_pointer> vsp(vptr, v.extensions());

    if (haj.size() == 1)
    {
      int nwalk = std::get<1>(v.sizes());
      if (walker_type == COLLINEAR)
      {
        assert(std::get<1>(G.sizes()) == std::get<1>(v.sizes()));
        int NMO, nel[2];
        NMO    = std::get<2>(Lnak[0].sizes());
        nel[0] = std::get<1>(Lnak[0].sizes());
        nel[1] = std::get<1>(Lnak[1].sizes());
        double c_[2];
        c_[0] = c;
        c_[1] = c;
        if (std::abs(c) < 1e-8)
          c_[1] = 1.0;
        assert((nel[0]+nel[1])*NMO == std::get<0>(G.sizes()));
        for (int ispin = 0, is0 = 0; ispin < 2; ispin++)
        {
          assert(std::get<0>(Lnak[ispin].sizes()) == std::get<0>(v.sizes()));
          SpCMatrix_ref Ln(make_device_ptr(Lnak[ispin].origin()), {local_nCV, nel[ispin] * NMO});
          ma::product(SPComplexType(a), Ln, Gsp.sliced(is0, is0 + nel[ispin] * NMO), SPComplexType(c_[ispin]), vsp);
          is0 += nel[ispin] * NMO;
        }
      }
      else
      {
        assert(std::get<1>(G.sizes()) == std::get<1>(v.sizes()));
        assert(std::get<1>(Lnak[0].sizes()) * std::get<2>(Lnak[0].sizes()) == std::get<0>(G.sizes()));
        assert(std::get<0>(Lnak[0].sizes()) == std::get<0>(v.sizes()));
        SpCMatrix_ref Ln(make_device_ptr(Lnak[0].origin()), {local_nCV, std::get<1>(Lnak[0].sizes()) * std::get<2>(Lnak[0].sizes())});
        ma::product(SPComplexType(a), Ln, Gsp, SPComplexType(c), vsp);
      }
    }
    else
    {
      // multideterminant is not half-rotated, so use Likn
      assert(std::get<0>(Likn.sizes()) == std::get<0>(G.sizes()));
      assert(std::get<1>(Likn.sizes()) == std::get<0>(v.sizes()));
      assert(std::get<1>(G.sizes()) == std::get<1>(v.sizes()));

      ma::product(SPValueType(a), ma::T(Likn), Gsp, SPValueType(c), vsp);
    }
    if (not std::is_same<vType, SPComplexType>::value)
    {
      copy_n_cast(make_device_ptr(vsp.origin()), v.num_elements(), make_device_ptr(v.origin()));
    }
  }

  template<class Mat, class MatB>
  void generalizedFockMatrix(Mat&& G, MatB&& Fp, MatB&& Fm)
  {
    int nwalk = G.size();
    int nspin = (walker_type == COLLINEAR ? 2 : 1);
    int NMO   = hij.size();
    int nel[2];
    assert(Fp.size() == nwalk);
    assert(Fm.size() == nwalk);
    assert(G[0].num_elements() == nspin * NMO * NMO);
    assert(Fp[0].num_elements() == nspin * NMO * NMO);
    assert(Fm[0].num_elements() == nspin * NMO * NMO);

    // Rwn[nwalk][nCV]: 1+nspin copies
    // Twpqn[nwalk][NMO][NMO][nCV]: 1+nspin copies
    // extra copies

    // can you find out how much memory is available on the buffer?
    long LBytes = max_memory_MB * 1024L * 1024L / long(sizeof(SPComplexType));
#if defined(MIXED_PRECISION)
    LBytes -= long((3 * nspin + 1) * nwalk * NMO * NMO); // G, Fp, Fm and Gt
#else
    LBytes -= long((1 + nspin) * nwalk * NMO * NMO); //  G and Gt
#endif
    LBytes *= long(sizeof(SPComplexType));
    int Bytes = int(LBytes / long(2 * (NMO * NMO + 1) * local_nCV * sizeof(SPComplexType)));
    int nwmax = std::min(std::max(1, Bytes), nwalk);
    assert(nwmax >= 1 && nwmax <= nwmax);

#if defined(MIXED_PRECISION)
    StaticMatrix Fp_({nwalk, nspin * NMO * NMO},
                     buffer_manager.get_generator().template get_allocator<SPComplexType>());
    StaticMatrix Fm_({nwalk, nspin * NMO * NMO},
                     buffer_manager.get_generator().template get_allocator<SPComplexType>());
#else
    SpCMatrix_ref Fp_(make_device_ptr(Fp.origin()), {nwalk, nspin * NMO * NMO});
    SpCMatrix_ref Fm_(make_device_ptr(Fm.origin()), {nwalk, nspin * NMO * NMO});
#endif
    fill_n(Fp_.origin(), Fp_.num_elements(), SPComplexType(0.0));
    fill_n(Fm_.origin(), Fm_.num_elements(), SPComplexType(0.0));

    SPComplexType scl = (walker_type == CLOSED ? 2.0 : 1.0);
    std::vector<sp_pointer> Aarray;
    std::vector<sp_pointer> Barray;
    std::vector<sp_pointer> Carray;
    Aarray.reserve(nwalk);
    Barray.reserve(nwalk);
    Carray.reserve(nwalk);

    long gsz(0);
#if defined(MIXED_PRECISION)
    gsz = nspin * nwmax * NMO * NMO;
#else
    if (nspin > 1)
      gsz = nspin * nwmax * NMO * NMO;
#endif
    StaticVector GBuff(iextensions<1u>{gsz}, buffer_manager.get_generator().template get_allocator<SPComplexType>());

    int nw0(0);
    while (nw0 < nwalk)
    {
      int nw = std::min(nwalk - nw0, nwmax);

      sp_pointer ptr(nullptr);
      // transpose/cast G
#if defined(MIXED_PRECISION)
      ptr = GBuff.origin();
      for (int ispin = 0, is0 = 0, ip = 0; ispin < nspin; ispin++, is0 += NMO * NMO)
        for (int n = 0; n < nw; ++n, ip += NMO * NMO)
          copy_n_cast(make_device_ptr(G[nw0 + n].origin()) + is0, NMO * NMO, ptr + ip);
#else
      if (nspin == 1)
      {
        ptr = make_device_ptr(G[nw0].origin());
      }
      else
      {
        ptr = GBuff.origin();
        using std::copy_n;
        for (int ispin = 0, is0 = 0, ip = 0; ispin < nspin; ispin++, is0 += NMO * NMO)
          for (int n = 0; n < nw; ++n, ip += NMO * NMO)
            copy_n(make_device_ptr(G[nw0 + n].origin()) + is0, NMO * NMO, ptr + ip);
      }
#endif
      SpCTensor_ref GF(ptr, {nspin, nw * NMO, NMO}); // now contains G in the correct structure [spin][w][i][j]
      StaticMatrix Gt({NMO * NMO, nw}, buffer_manager.get_generator().template get_allocator<SPComplexType>());
      fill_n(Gt.origin(), Gt.num_elements(), SPComplexType(0.0));

      StaticMatrix Rnw({local_nCV, nw}, buffer_manager.get_generator().template get_allocator<SPComplexType>());
      // calculate Rwn
      for (int ispin = 0; ispin < nspin; ispin++)
      {
        SpCMatrix_ref G_(GF[ispin].origin(), {nw, NMO * NMO});
        ma::add(SPComplexType(1.0), Gt, SPComplexType(1.0), ma::T(G_), Gt);
      }
      // R[n,w] = \sum_ik L[n,ik] G[ik,w]
      ma::product(SPValueType(1.0), ma::T(Likn), Gt, SPValueType(0.0), Rnw);
      StaticMatrix Rwn({nw, local_nCV}, buffer_manager.get_generator().template get_allocator<SPComplexType>());
      ma::transpose(Rnw, Rwn);

      // add coulomb contribution of <pr||qs>Grs term to Fp, reuse Gt for temporary storage
      // Fp[p,t] = \sum_{jl} L[p,t,n] L[j,l,n] P[j,l]
      // Fp[pt,w] = \sum_n L[pt,n] R[n,w]
      ma::product(SPValueType(1.0), Likn, Rnw, SPValueType(0.0), Gt);
      for (int ispin = 0; ispin < nspin; ispin++)
      {
        ma::add(SPComplexType(1.0), Fp_({nw0, nw0 + nw}, {ispin * NMO * NMO, (ispin + 1) * NMO * NMO}),
                SPComplexType(scl), ma::T(Gt), Fp_({nw0, nw0 + nw}, {ispin * NMO * NMO, (ispin + 1) * NMO * NMO}));
      }

      // L[i,kn]
      SpRMatrix_ref Ln(make_device_ptr(Likn.origin()), {NMO, NMO * local_nCV});
      // T[w,p,t,n] = \sum_{l} L[l,t,n] P[w,l,p]
      StaticMatrix Twptn({nw * NMO, NMO * local_nCV},
                         buffer_manager.get_generator().template get_allocator<SPComplexType>());
      // transpose for faster contraction
      StaticMatrix Taux({nw * NMO, NMO * local_nCV},
                        buffer_manager.get_generator().template get_allocator<SPComplexType>());
      SpCTensor_ref Taux3D(Taux.origin(), {nw, NMO, NMO * local_nCV});
      SpCTensor_ref Twptn3D(Twptn.origin(), {nw, NMO, NMO * local_nCV});
      SpCTensor_ref Twptn3D_(Twptn.origin(), {nw, NMO * NMO, local_nCV});
      SpCMatrix_ref Ttnwp(Taux.origin(), {NMO * local_nCV, nw * NMO});
      SpCMatrix_ref Gt_(Gt.origin(), {NMO, nw * NMO});

      for (int ispin = 0, is0 = 0; ispin < nspin; ispin++, is0 += NMO * NMO)
      {
        SpCMatrix_ref G_(GF[ispin].origin(), {nw * NMO, NMO});
        ma::transpose(G_, Gt_);

        // J = \sum_{iklr} L[i,k,n] L[q,l,n] P[s,p,l] P[r,i,k]
        // R[n] = \sum_{ik} L[i,k,n] P[r,i,k]
        // Here T[tn,wp] = \sum_{l} L[tn,l] P[l,wp]
        // T[ln,p] = T[npl] = L[nkl] P[p,l]
        ma::product(SPValueType(1.0), ma::T(Ln), Gt_, SPValueType(0.0), Ttnwp);
        // T[wp,tn]
        ma::transpose(Ttnwp, Twptn);

        // transpose Twptn -> Twtpn=Taux
        using ma::transpose_wabn_to_wban;
        // T[wt,pn]
        transpose_wabn_to_wban(nw, NMO, NMO, local_nCV, Twptn.origin(), Taux.origin());

        // add exchange component to Fm_
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int w = 0; w < nw; w++)
        {
          Aarray.push_back(Taux3D[w].origin());
          Barray.push_back(Twptn3D[w].origin());
          Carray.push_back(Fm_[w].origin() + is0);
        }
        using ma::gemmBatched;
        // careful with expected Fortran ordering here!!!
        // K[p,q] = \sum_{ln} T[n,l,p] T[n,q,l]
        //          \sum_{ln} T[nl,p] T[nl,q]
        gemmBatched('T', 'N', NMO, NMO, NMO * local_nCV, SPComplexType(1.0), Aarray.data(), NMO * local_nCV,
                    Barray.data(), NMO * local_nCV, SPComplexType(1.0), Carray.data(), NMO, nw);

        // add coulomb component to Fm_
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int w = 0; w < nw; w++)
        {
          Aarray.push_back(Twptn3D_[w].origin());
          Barray.push_back(Rwn[w].origin());
          Carray.push_back(Fm_[w].origin() + is0);
        }
        using ma::gemmBatched;
        // careful with expected Fortran ordering here!!!
        // J[w][pq] = \sum_{n} T[w][pq,n] R[w][n]
        gemmBatched('T', 'N', NMO * NMO, 1, local_nCV, SPComplexType(-1.0) * scl, Aarray.data(), local_nCV,
                    Barray.data(), local_nCV, SPComplexType(1.0), Carray.data(), NMO * NMO, nw);

        // Fp
        // Need Gt_[i][wj]
        transpose_wabn_to_wban(1, nw, NMO, NMO, G_.origin(), Gt_.origin());
        ma::product(SPValueType(1.0), ma::T(Ln), Gt_, SPValueType(0.0), Ttnwp);
        ma::transpose(Ttnwp, Twptn);
        transpose_wabn_to_wban(nw, NMO, NMO, local_nCV, Twptn.origin(), Taux.origin());
        // add coulomb component to Fp_, same as Fm_ above
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int w = 0; w < nw; w++)
        {
          Aarray.push_back(Twptn3D_[w].origin());
          Barray.push_back(Rwn[w].origin());
          Carray.push_back(Fp_[w].origin() + is0);
        }
        using ma::gemmBatched;
        // careful with expected Fortran ordering here!!!
        // Coulomb component
        gemmBatched('T', 'N', NMO * NMO, 1, local_nCV, SPComplexType(-1.0) * scl, Aarray.data(), local_nCV,
                    Barray.data(), local_nCV, SPComplexType(1.0), Carray.data(), NMO * NMO, nw);

        // add exchange component of Fp_
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int w = 0; w < nw; w++)
        {
          Aarray.push_back(Taux3D[w].origin());
          Barray.push_back(Twptn3D[w].origin());
          Carray.push_back(Fp_[w].origin() + is0);

          // add exchange contribution of <pr||qs>Grs term by adding Lptn to Twptn
          // dispatch directly from here to be able to add to the real part only
          // K1B[p,q] = -\sum_{jl} L[jt,n] L[pl,n] P[jl]
          using ma::axpy;
          axpy(Likn.num_elements(), SPRealType(-1.0), ma::pointer_dispatch(Likn.origin()), 1,
               pointer_cast<SPRealType>(ma::pointer_dispatch(Twptn3D_[w].origin())), 2);
        }
        using ma::gemmBatched;
        // careful with expected Fortran ordering here!!!
        gemmBatched('T', 'N', NMO, NMO, NMO * local_nCV, SPComplexType(1.0), Aarray.data(), NMO * local_nCV,
                    Barray.data(), NMO * local_nCV, SPComplexType(1.0), Carray.data(), NMO, nw);

      } // ispin

      nw0 += nw;
    }

#if defined(MIXED_PRECISION)
    copy_n_cast(Fp_.origin(), Fp_.num_elements(), make_device_ptr(Fp.origin()));
    copy_n_cast(Fm_.origin(), Fm_.num_elements(), make_device_ptr(Fm.origin()));
#endif

    //fill_n(Fp.origin(),Fp.num_elements(),SPComplexType(0.0));
    //fill_n(Fm.origin(),Fm.num_elements(),SPComplexType(0.0));
    // add one body terms now
    {
      std::vector<pointer> Aarr;
      std::vector<pointer> Barr;
      std::vector<pointer> Carr;
      Aarr.reserve(nspin * nwalk);
      Barr.reserve(nspin * nwalk);
      Carr.reserve(nspin * nwalk);
      // Fm -= G[w][p][r] * h[q][r]
      Aarr.clear();
      Barr.clear();
      Carr.clear();
      for (int ispin = 0, is0 = 0; ispin < nspin; ispin++, is0 += NMO * NMO)
      {
        for (int w = 0; w < nwalk; w++)
        {
          Aarr.push_back(make_device_ptr(hij_dev.origin()));
          Barr.push_back(make_device_ptr(G[w].origin()) + is0);
          Carr.push_back(make_device_ptr(Fm[w].origin()) + is0);
        }
      }
      using ma::gemmBatched;
      // careful with expected Fortran ordering here!!!
      gemmBatched('T', 'N', NMO, NMO, NMO, ComplexType(-1.0), Aarr.data(), NMO, Barr.data(), NMO, ComplexType(1.0),
                  Carr.data(), NMO, Aarr.size());


      // Fp -= G[w][r][p] * h[q][r]
      Aarr.clear();
      Barr.clear();
      Carr.clear();
      C4Tensor_ref Fp4D(make_device_ptr(Fp.origin()), {nwalk, nspin, NMO, NMO});
      for (int ispin = 0, is0 = 0; ispin < nspin; ispin++, is0 += NMO * NMO)
      {
        for (int w = 0; w < nwalk; w++)
        {
          Aarr.push_back(make_device_ptr(hij_dev.origin()));
          Barr.push_back(make_device_ptr(G[w].origin()) + is0);
          Carr.push_back(make_device_ptr(Fp[w].origin()) + is0);

          // add diagonal contribution to Fp
          ma::add(ComplexType(1.0), Fp4D[w][ispin], ComplexType(1.0), ma::T(hij_dev), Fp4D[w][ispin]);
        }
      }
      using ma::gemmBatched;
      // careful with expected Fortran ordering here!!!
      gemmBatched('T', 'T', NMO, NMO, NMO, ComplexType(-1.0), Aarr.data(), NMO, Barr.data(), NMO, ComplexType(1.0),
                  Carr.data(), NMO, Aarr.size());
    }
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
  DeviceBufferManager buffer_manager;

  WALKER_TYPES walker_type;

  long max_memory_MB;
  int global_origin;
  int global_nCV;
  int local_nCV;

  ValueType E0;

  // bare one body hamiltonian
  mpi3RMatrix hij;

  // one body hamiltonian
  shmCMatrix hij_dev;

  // (potentially half rotated) one body hamiltonian
  shmCMatrix haj;

  //Cholesky Tensor Lik[i][k][n]
  shmSpRMatrix Likn;

  // permuted half-transformed Cholesky tensor
  // Lnak[ 2*idet + ispin ]
  std::vector<shmSpC3Tensor> Lnak;

  // one-body piece of Hamiltonian factorization
  mpi3CMatrix vn0;

  myTimer Timer;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
