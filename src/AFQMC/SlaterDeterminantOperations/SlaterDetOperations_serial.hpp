//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SERIAL_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SERIAL_HPP

#include <fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"
#include "AFQMC/SlaterDeterminantOperations/apply_expM.hpp"

#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations_base.hpp"

#include "mpi3/shared_communicator.hpp"
#include "type_traits/complex_help.hpp"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
// Implementation that doesn't use shared memory
// This version is designed for GPU and/or threading with OpenMP/libraries
template<class Type, class BufferManager>
class SlaterDetOperations_serial : public SlaterDetOperations_base<Type, BufferManager>
{
public:
  using communicator = boost::mpi3::shared_communicator;
  using Base         = SlaterDetOperations_base<Type, BufferManager>;

  using T             = typename Base::T;
  using buffer_type_T = typename Base::buffer_type_T;
  using IVector       = typename Base::IVector;
  using TVector       = typename Base::TVector;
  using TMatrix       = typename Base::TMatrix;
  using TTensor       = boost::multi::static_array<T, 3, buffer_type_T>;

  using Base::MixedDensityMatrix;
  using Base::MixedDensityMatrix_noHerm;
  using Base::MixedDensityMatrixForWoodbury;
  using Base::MixedDensityMatrixFromConfiguration;
  using Base::Orthogonalize;
  using Base::Overlap;
  using Base::Overlap_noHerm;
  using Base::OverlapForWoodbury;
  using Base::Propagate;

  SlaterDetOperations_serial(BufferManager b) : SlaterDetOperations_base<Type, BufferManager>(b) {}

  SlaterDetOperations_serial(int NMO, int NAEA, BufferManager b)
      : SlaterDetOperations_base<Type, BufferManager>(NMO, NAEA, b)
  {}

  ~SlaterDetOperations_serial() {}

  SlaterDetOperations_serial(const SlaterDetOperations_serial& other) = delete;
  SlaterDetOperations_serial(SlaterDetOperations_serial&& other)      = default;
  SlaterDetOperations_serial& operator=(const SlaterDetOperations_serial& other) = delete;
  SlaterDetOperations_serial& operator=(SlaterDetOperations_serial&& other) = default;

  // C must live in shared memory for this routine to work as expected
  template<class MatA, class MatB, class MatC>
  T MixedDensityMatrix(const MatA& hermA,
                       const MatB& B,
                       MatC&& C,
                       T LogOverlapFactor,
                       communicator& comm,
                       bool compact = false,
                       bool herm    = true)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
    return Base::MixedDensityMatrix(hermA, B, C, LogOverlapFactor, compact);
  }

  template<class integer, class MatA, class MatB, class MatC, class MatQ>
  T MixedDensityMatrixForWoodbury(const MatA& hermA,
                                  const MatB& B,
                                  MatC&& C,
                                  T LogOverlapFactor,
                                  integer* ref,
                                  MatQ&& QQ0,
                                  communicator& comm,
                                  bool compact = false)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
    return Base::MixedDensityMatrixForWoodbury(hermA, B, std::forward<MatC>(C), LogOverlapFactor, ref,
                                               std::forward<MatQ>(QQ0), compact);
  }

  template<class MatA, class MatB>
  T Overlap(const MatA& hermA, const MatB& B, T LogOverlapFactor, communicator& comm, bool herm = true)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
    return Base::Overlap(hermA, B, LogOverlapFactor);
  }

  template<typename integer, class MatA, class MatB, class MatC>
  T OverlapForWoodbury(const MatA& hermA,
                       const MatB& B,
                       T LogOverlapFactor,
                       integer* ref,
                       MatC&& QQ0,
                       communicator& comm)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
    return Base::OverlapForWoodbury(hermA, B, LogOverlapFactor, ref, std::forward<MatC>(QQ0));
  }

  template<class Mat, class MatP1, class MatV>
  void Propagate(Mat&& A,
                 const MatP1& P1,
                 const MatV& V,
                 communicator& comm,
                 int order         = 6,
                 char TA           = 'N',
                 bool noncollinear = false)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
    Base::Propagate(std::forward<Mat>(A), P1, V, order, TA, noncollinear);
  }

  template<class MatA, class MatP1, class MatV>
  void BatchedPropagate(std::vector<MatA>& Ai,
                        const MatP1& P1,
                        const MatV& V,
                        int order         = 6,
                        char TA           = 'N',
                        bool noncollinear = false)
  {
    static_assert(pointedType<MatA>::dimensionality == 2, " dimenionality == 2");
    static_assert(std::decay<MatV>::type::dimensionality == 3, " dimenionality == 3");
    if (Ai.size() == 0)
      return;
    assert(Ai.size() == std::get<0>(V.sizes()));
    int nbatch = Ai.size();
    int npol   = noncollinear ? 2 : 1;
    int NMO    = std::get<0>((*Ai[0]).sizes());
    int NAEA   = std::get<1>((*Ai[0]).sizes());
    int M      = NMO / npol;
    assert(NMO % npol == 0);
    assert(std::get<0>(P1.sizes()) == NMO);
    assert(std::get<1>(P1.sizes()) == NMO);
    assert(std::get<1>(V.sizes()) == M);
    assert(std::get<2>(V.sizes()) == M);
    TTensor TMN({nbatch, NMO, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    TTensor T1({nbatch, NMO, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    TTensor T2({nbatch, NMO, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    // could be batched when csrmm is batched
    if (TA == 'H' || TA == 'h')
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(ma::H(P1), *Ai[ib], TMN[ib]);
    }
    else if (TA == 'T' || TA == 't')
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(ma::T(P1), *Ai[ib], TMN[ib]);
    }
    else
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(P1, *(Ai[ib]), TMN[ib]);
    }

    // Apply V
    if (noncollinear)
    {
      // treat 2 polarizations as separate elements in the batch
      using TTensor_ref = boost::multi::array_ref<T, 3, decltype(TMN.origin())>;
      TTensor_ref TMN_(TMN.origin(), {npol * nbatch, M, NAEA});
      TTensor_ref T1_(T1.origin(), {npol * nbatch, M, NAEA});
      TTensor_ref T2_(T2.origin(), {npol * nbatch, M, NAEA});
      SlaterDeterminantOperations::batched::apply_expM_noncollinear(V, TMN_, T1_, T2_, order, TA);
    }
    else
    {
      SlaterDeterminantOperations::batched::apply_expM(V, TMN, T1, T2, order, TA);
    }

    if (TA == 'H' || TA == 'h')
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(ma::H(P1), TMN[ib], *Ai[ib]);
    }
    else if (TA == 'T' || TA == 't')
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(ma::T(P1), TMN[ib], *Ai[ib]);
    }
    else
    {
      for (int ib = 0; ib < nbatch; ib++)
        ma::product(P1, TMN[ib], *Ai[ib]);
    }
  }

  // C[nwalk, M, N]
  template<class MatA, class MatB, class MatC, class TVec>
  void BatchedMixedDensityMatrix(std::vector<MatA>& hermA,
                                 std::vector<MatB>& Bi,
                                 MatC&& C,
                                 T LogOverlapFactor,
                                 TVec&& ovlp,
                                 bool compact = false,
                                 bool herm    = true)
  {
    if (Bi.size() == 0)
      return;
    static_assert((pointedType<MatA>::dimensionality == 2 or pointedType<MatA>::dimensionality == -2),
                  "Wrong dimensionality");
    static_assert(pointedType<MatB>::dimensionality == 2, "Wrong dimensionality");
    static_assert(std::decay<MatC>::type::dimensionality == 3, "Wrong dimensionality");
    static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
    int NMO    = (herm ? std::get<1>((*hermA[0]).sizes()) : std::get<0>((*hermA[0]).sizes()));
    int NAEA   = (herm ? std::get<0>((*hermA[0]).sizes()) : std::get<1>((*hermA[0]).sizes()));
    int nbatch = Bi.size();
    assert(C.size() == nbatch);
    assert(ovlp.size() == nbatch);
    int n1 = nbatch, n2 = NAEA, n3 = NMO;
    if (compact)
    {
      n1 = n2 = n3 = 0;
    }
    TTensor TNN3D({nbatch, NAEA, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    TTensor TNM3D({n1, n2, n3}, buffer_manager.get_generator().template get_allocator<T>());
    IVector IWORK(iextensions<1u>{nbatch * (NMO + 1)}, buffer_manager.get_generator().template get_allocator<int>());
    SlaterDeterminantOperations::batched::MixedDensityMatrix(hermA, Bi, std::forward<MatC>(C), LogOverlapFactor,
                                                             std::forward<TVec>(ovlp), TNN3D, TNM3D, IWORK, compact,
                                                             herm);
  }

  template<class MatA, class MatB, class MatC, class TVec>
  void BatchedDensityMatrices(const std::vector<MatA>& Left,
                              const std::vector<MatB>& Right,
                              std::vector<MatC>& G,
                              T LogOverlapFactor,
                              TVec&& ovlp,
                              bool compact = false,
                              bool herm    = true)
  {
    if (Left.size() == 0)
      return;
    if (Right.size() == 0)
      return;
    static_assert((pointedType<MatA>::dimensionality == 2 or pointedType<MatA>::dimensionality == -2),
                  "Wrong dimensionality");
    static_assert(pointedType<MatB>::dimensionality == 2, "Wrong dimensionality");
    static_assert(pointedType<MatC>::dimensionality == 2, "Wrong dimensionality");
    static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
    int NMO    = (herm ? std::get<1>((*Left[0]).sizes()) : std::get<0>((*Left[0]).sizes()));
    int NAEA   = (herm ? std::get<0>((*Left[0]).sizes()) : std::get<1>((*Left[0]).sizes()));
    int nbatch = Left.size();
    assert(Right.size() == nbatch);
    assert(G.size() == nbatch);
    assert(ovlp.size() == nbatch);
    int n1 = nbatch, n2 = NAEA, n3 = NMO;
    if (compact)
    {
      n1 = n2 = n3 = 0;
    }
    TTensor TNN3D({nbatch, NAEA, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    TTensor TNM3D({n1, n2, n3}, buffer_manager.get_generator().template get_allocator<T>());
    IVector IWORK(iextensions<1u>{nbatch * (NMO + 1)}, buffer_manager.get_generator().template get_allocator<int>());
    SlaterDeterminantOperations::batched::DensityMatrices(Left, Right, G, LogOverlapFactor, std::forward<TVec>(ovlp),
                                                          TNN3D, TNM3D, IWORK, compact, herm);
  }

  template<class MatA, class MatB, class TVec>
  void BatchedOverlap(std::vector<MatA>& hermA,
                      std::vector<MatB>& Bi,
                      T LogOverlapFactor,
                      TVec&& ovlp,
                      bool herm = true)
  {
    static_assert((pointedType<MatA>::dimensionality == 2 or pointedType<MatA>::dimensionality == -2),
                  "Wrong dimensionality");
    static_assert(pointedType<MatB>::dimensionality == 2, "Wrong dimensionality");
    if (Bi.size() == 0)
      return;
    assert(hermA.size() > 0);
    static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
    int NMO    = (herm ? std::get<1>((*hermA[0]).sizes()) : std::get<0>((*hermA[0]).sizes()));
    int NAEA   = (herm ? std::get<0>((*hermA[0]).sizes()) : std::get<1>((*hermA[0]).sizes()));
    int nbatch = Bi.size();
    assert(ovlp.size() == nbatch);
    TTensor TNN3D({nbatch, NAEA, NAEA}, buffer_manager.get_generator().template get_allocator<T>());
    IVector IWORK(iextensions<1u>{nbatch * (NMO + 1)}, buffer_manager.get_generator().template get_allocator<int>());
    SlaterDeterminantOperations::batched::Overlap(hermA, Bi, LogOverlapFactor, std::forward<TVec>(ovlp), TNN3D, IWORK,
                                                  herm);
  }

  template<class MatA, class PTR>
  void BatchedOrthogonalize(std::vector<MatA>& Ai, T LogOverlapFactor, PTR detR)
  {
    static_assert(pointedType<MatA>::dimensionality == 2, "Wrong dimensionality");
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    // QR on the transpose
    if (Ai.size() == 0)
      return;
    int NMO    = std::get<0>((*Ai[0]).sizes());
    int NAEA   = std::get<1>((*Ai[0]).sizes());
    int nbatch = Ai.size();
    TTensor AT({nbatch, NAEA, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    TMatrix T_({nbatch, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    TMatrix scl({nbatch, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    int sz = ma::gqr_optimal_workspace_size(AT[0]);
    TVector WORK(iextensions<1u>{nbatch * sz}, buffer_manager.get_generator().template get_allocator<T>());
    IVector IWORK(iextensions<1u>{nbatch * (NMO + 1)}, buffer_manager.get_generator().template get_allocator<int>());
    for (int i = 0; i < nbatch; i++)
      ma::transpose(*Ai[i], AT[i]);
    // careful, expects fortran order
    geqrfStrided(NMO, NAEA, AT.origin(), NMO, NMO * NAEA, T_.origin(), NMO, IWORK.origin(), nbatch);
    using ma::determinant_from_geqrf;
    using ma::scale_columns;
    for (int i = 0; i < nbatch; i++)
      *(detR + i) = determinant_from_geqrf(NAEA, AT[i].origin(), NMO, scl[i].origin(), LogOverlapFactor);
    gqrStrided(NMO, NAEA, NAEA, AT.origin(), NMO, NMO * NAEA, T_.origin(), NMO, WORK.origin(), sz, IWORK.origin(),
               nbatch);
    for (int i = 0; i < nbatch; i++)
    {
      ma::transpose(AT[i], *Ai[i]);
      scale_columns(NMO, NAEA, (*Ai[i]).origin(), (*Ai[i]).stride(0), scl[i].origin());
    }
#else
    int nw = Ai.size();
    for (int i = 0; i < nw; i++)
      *(detR + i) = Orthogonalize(*Ai[i], LogOverlapFactor);
#endif
  }

  template<class MatA>
  void BatchedOrthogonalize(std::vector<MatA>& Ai, T LogOverlapFactor)
  {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    // QR on the transpose
    if (Ai.size() == 0)
      return;
    int NMO    = std::get<0>((*Ai[0]).sizes());
    int NAEA   = std::get<1>((*Ai[0]).sizes());
    int nbatch = Ai.size();
    TTensor AT({nbatch, NAEA, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    TMatrix T_({nbatch, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    TMatrix scl({nbatch, NMO}, buffer_manager.get_generator().template get_allocator<T>());
    int sz = ma::gqr_optimal_workspace_size(AT[0]);
    TVector WORK(iextensions<1u>{nbatch * sz}, buffer_manager.get_generator().template get_allocator<T>());
    IVector IWORK(iextensions<1u>{nbatch * (NMO + 1)}, buffer_manager.get_generator().template get_allocator<int>());
    for (int i = 0; i < nbatch; i++)
      ma::transpose(*Ai[i], AT[i]);
    // careful, expects fortran order
    geqrfStrided(NMO, NAEA, AT.origin(), NMO, NMO * NAEA, T_.origin(), NMO, IWORK.origin(), nbatch);
    using ma::determinant_from_geqrf;
    using ma::scale_columns;
    for (int i = 0; i < nbatch; i++)
      determinant_from_geqrf(NAEA, AT[i].origin(), NMO, scl[i].origin());
    gqrStrided(NMO, NAEA, NAEA, AT.origin(), NMO, NMO * NAEA, T_.origin(), NMO, WORK.origin(), sz, IWORK.origin(),
               nbatch);
    for (int i = 0; i < nbatch; i++)
    {
      ma::transpose(AT[i], *Ai[i]);
      scale_columns(NMO, NAEA, (*Ai[i]).origin(), (*Ai[i]).stride(0), scl[i].origin());
    }
#else
    int nw = Ai.size();
    for (int i = 0; i < nw; i++)
      Orthogonalize(*Ai[i], LogOverlapFactor);
#endif
  }

protected:
  using Base::buffer_manager;
  using Base::work_size;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
