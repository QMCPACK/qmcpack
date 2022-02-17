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

#ifndef QMCPLUSPLUS_AFQMC_GENERALIZEDFOCKMATRIX_HPP
#define QMCPLUSPLUS_AFQMC_GENERALIZEDFOCKMATRIX_HPP

#include "AFQMC/config.h"
#include <vector>
#include <string>
#include <iostream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/Timer.h"

#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
/* 
 * Observable class that calculates the walker averaged "full" 1 RDM.
 * In this context, "full" means that no contraction over the RDM is
 * being performed. The resulting RDM will be [spin][x*NMO][x*NMO],
 * where x:2 for NONCOLLINEAR and 1 for everything else.
 */
class generalizedFockMatrix : public AFQMCInfo
{
  // allocators
  using Allocator        = device_allocator<ComplexType>;
  using Allocator_shared = localTG_allocator<ComplexType>;

  // type defs
  using pointer              = typename Allocator::pointer;
  using const_pointer        = typename Allocator::const_pointer;
  using pointer_shared       = typename Allocator_shared::pointer;
  using const_pointer_shared = typename Allocator_shared::const_pointer;

  using CVector_ref    = boost::multi::array_ref<ComplexType, 1, pointer>;
  using CMatrix_ref    = boost::multi::array_ref<ComplexType, 2, pointer>;
  using CVector        = boost::multi::array<ComplexType, 1, Allocator>;
  using CMatrix        = boost::multi::array<ComplexType, 2, Allocator>;
  using sharedCMatrix  = boost::multi::array<ComplexType, 2, Allocator_shared>;
  using sharedCTensor  = boost::multi::array<ComplexType, 3, Allocator_shared>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType, 1>;
  using stdCMatrix_ref = boost::multi::array_ref<ComplexType, 2>;
  using stdCVector     = boost::multi::array<ComplexType, 1>;
  using stdCMatrix     = boost::multi::array<ComplexType, 2>;
  using stdIMatrix     = boost::multi::array<int, 2>;
  using mpi3CVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using mpi3IMatrix    = boost::multi::array<int, 2, shared_allocator<int>>;
  using mpi3CMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using mpi3CTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using mpi3C4Tensor   = boost::multi::array<ComplexType, 4, shared_allocator<ComplexType>>;

  using stack_alloc_type = DeviceBufferManager::template allocator_t<ComplexType>;
  using Static3Tensor    = boost::multi::static_array<ComplexType, 3, stack_alloc_type>;

public:
  generalizedFockMatrix(afqmc::TaskGroup_& tg_,
                        AFQMCInfo& info,
                        xmlNodePtr cur,
                        WALKER_TYPES wlk,
                        HamiltonianOperations* HOps_,
                        int nave_ = 1,
                        int bsize = 1)
      : AFQMCInfo(info),
        TG(tg_),
        walker_type(wlk),
        HamOp(HOps_),
        writer(false),
        block_size(bsize),
        nave(nave_),
        counter(0),
        denom(iextensions<1u>{0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMAverage({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    using std::copy_n;
    using std::fill_n;

    app_log() << "  --  Adding Back Propagated generalized Fock matrix (genFock) estimator. -- \n";

    dm_size = NMO * NMO;
    if (walker_type == COLLINEAR)
      dm_size *= 2;
    else if (walker_type == NONCOLLINEAR)
      dm_size *= 4;

    using std::fill_n;
    writer = (TG.getGlobalRank() == 0);

    DMAverage = mpi3CTensor({2, nave, dm_size}, shared_allocator<ComplexType>{TG.TG_local()});
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0, 0.0));
  }

  template<class MatG, class MatG_host, class HostCVec1, class HostCVec2, class HostCVec3>
  void accumulate_reference(int iav,
                            int iref,
                            MatG&& G,
                            MatG_host&& G_host,
                            HostCVec1&& wgt,
                            HostCVec2&& Xw,
                            HostCVec3&& ovlp,
                            bool impsamp)
  {
    static_assert(std::decay<MatG>::type::dimensionality == 4, "Wrong dimensionality");
    static_assert(std::decay<MatG_host>::type::dimensionality == 4, "Wrong dimensionality");
    using std::fill_n;
    // assumes G[nwalk][spin][M][M]
    int nw(G.size());
    assert(G.size() == wgt.size());
    assert(wgt.size() == nw);
    assert(Xw.size() == nw);
    assert(ovlp.size() >= nw);
    assert(G.num_elements() == G_host.num_elements());
    assert(G.extensions() == G_host.extensions());
    assert(G[0].num_elements() == dm_size);

    // check structure dimensions
    if (iref == 0)
    {
      if (denom.size() != nw)
      {
        denom = mpi3CVector(iextensions<1u>{nw}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(DMWork.sizes()) != 3 || std::get<1>(DMWork.sizes()) != nw || std::get<2>(DMWork.sizes()) != dm_size)
      {
        DMWork = mpi3CTensor({3, nw, dm_size}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      fill_n(denom.origin(), denom.num_elements(), ComplexType(0.0, 0.0));
      fill_n(DMWork.origin(), DMWork.num_elements(), ComplexType(0.0, 0.0));
    }
    else
    {
      if (std::get<0>(denom.sizes()) != nw || std::get<0>(DMWork.sizes()) != 2 || std::get<1>(DMWork.sizes()) != nw || std::get<2>(DMWork.sizes()) != dm_size ||
          std::get<0>(DMAverage.sizes()) != 2 || std::get<1>(DMAverage.sizes()) != nave || std::get<2>(DMAverage.sizes()) != dm_size)
        APP_ABORT(" Error: Invalid state in accumulate_reference. \n\n\n");
    }

    DeviceBufferManager buffer_manager;
    Static3Tensor gFock({2, nw, dm_size}, buffer_manager.get_generator().template get_allocator<ComplexType>());

    HamOp->generalizedFockMatrix(G, gFock[0], gFock[1]);

    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), dm_size, TG.TG_local().size());

    if (TG.TG_local().root())
      for (int iw = 0; iw < nw; iw++)
        denom[iw] += Xw[iw];
    for (int ic = 0; ic < 2; ic++)
    {
      for (int iw = 0; iw < nw; iw++)
      {
        copy_n(make_device_ptr(gFock[ic][iw].origin()) + i0, (iN - i0), to_address(DMWork[2][iw].origin()) + i0);
        ma::axpy(Xw[iw], DMWork[2][iw].sliced(i0, iN), DMWork[ic][iw].sliced(i0, iN));
      }
    }
    TG.TG_local().barrier();
  }

  template<class HostCVec>
  void accumulate_block(int iav, HostCVec&& wgt, bool impsamp)
  {
    int nw(denom.size());
    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), dm_size, TG.TG_local().size());

    if (TG.TG_local().root())
      for (int iw = 0; iw < nw; iw++)
        denom[iw] = wgt[iw] / denom[iw];
    TG.TG_local().barrier();

    // DMAverage[iav][ij] += sum_iw DMWork[iw][ij] * denom[iw] = T( DMWork ) * denom
    ma::product(ComplexType(1.0, 0.0), ma::T(DMWork[0]({0, nw}, {i0, iN})), denom, ComplexType(1.0, 0.0),
                DMAverage[0][iav].sliced(i0, iN));
    ma::product(ComplexType(1.0, 0.0), ma::T(DMWork[1]({0, nw}, {i0, iN})), denom, ComplexType(1.0, 0.0),
                DMAverage[1][iav].sliced(i0, iN));
    TG.TG_local().barrier();
  }

  template<class HostCVec>
  void print(int iblock, hdf_archive& dump, HostCVec&& Wsum)
  {
    using std::fill_n;
    const int n_zero = 9;

    if (TG.TG_local().root())
    {
      ma::scal(ComplexType(1.0 / block_size), DMAverage);
      TG.TG_heads().reduce_in_place_n(to_address(DMAverage.origin()), DMAverage.num_elements(), std::plus<>(), 0);
      if (writer)
      {
        dump.push(std::string("GenFockPlus"));
        for (int i = 0; i < nave; ++i)
        {
          dump.push(std::string("Average_") + std::to_string(i));
          std::string padded_iblock =
              std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
          stdCVector_ref DMAverage_(to_address(DMAverage[0][i].origin()), {dm_size});
          dump.write(DMAverage_, "gfockp_" + padded_iblock);
          dump.write(Wsum[i], "denominator_" + padded_iblock);
          dump.pop();
        }
        dump.pop();
        dump.push(std::string("GenFockMinus"));
        for (int i = 0; i < nave; ++i)
        {
          dump.push(std::string("Average_") + std::to_string(i));
          std::string padded_iblock =
              std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
          stdCVector_ref DMAverage_(to_address(DMAverage[1][i].origin()), {dm_size});
          dump.write(DMAverage_, "gfockm_" + padded_iblock);
          dump.write(Wsum[i], "denominator_" + padded_iblock);
          dump.pop();
        }
        dump.pop();
      }
    }
    TG.TG_local().barrier();
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0, 0.0));
  }

private:
  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  HamiltonianOperations* HamOp;

  bool writer;

  int block_size;

  int nave;

  int counter;

  int dm_size;

  mpi3CVector denom;

  // DMAverage (2, nave, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CTensor DMAverage;
  // DMWork (3, nwalk, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CTensor DMWork;
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
