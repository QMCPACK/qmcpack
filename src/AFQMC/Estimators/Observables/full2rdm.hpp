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

#ifndef QMCPLUSPLUS_AFQMC_FULL2RDM_HPP
#define QMCPLUSPLUS_AFQMC_FULL2RDM_HPP

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
 * Observable class that calculates the walker averaged 2 RDM.
 * The resulting RDM will be [3*spin][i][k][j][l]  
 * where x:2 for NONCOLLINEAR and 1 for everything else.
 * For collinear, the spin ordering is (a,a,a,a), (a,a,b,b), (b,b,b,b) 
 */
class full2rdm : public AFQMCInfo
{
  // allocators
  using Allocator        = device_allocator<ComplexType>;
  using Allocator_shared = node_allocator<ComplexType>;

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
  using StaticVector     = boost::multi::static_array<ComplexType, 1, stack_alloc_type>;
  using StaticMatrix     = boost::multi::static_array<ComplexType, 2, stack_alloc_type>;

public:
  full2rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk, int nave_ = 1, int bsize = 1)
      : AFQMCInfo(info),
        TG(tg_),
        walker_type(wlk),
        writer(false),
        block_size(bsize),
        nave(nave_),
        counter(0),
        apply_rotation(false),
        XRot({0, 0}, make_node_allocator<ComplexType>(TG)),
        denom(iextensions<1u>{0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMAverage({0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork({0, 0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    using std::copy_n;
    using std::fill_n;

    app_log() << "  --  Adding Back Propagated 2RDM (TwoRDM) estimator. -- \n";

    std::string rot_file("");
    std::string path("/");
    std::string str("false");

    if (cur != NULL)
    {
      ParameterSet m_param;
      m_param.add(rot_file, "rotation");
      m_param.add(path, "path");
      m_param.put(cur);
    }

    if (rot_file != "")
    {
      if (not file_exists(rot_file))
      {
        app_error() << " Error: File with rotation matrix does not exist: " << rot_file << std::endl;
        APP_ABORT("");
      }
      apply_rotation = true;
      int dim[2];

      hdf_archive dump;
      if (TG.Node().root())
      {
        if (!dump.open(rot_file, H5F_ACC_RDONLY))
          APP_ABORT("Error opening orbitals file for n2r estimator.\n");
        dump.push(path, false);

        stdCMatrix R;
        if (!dump.readEntry(R, "RotationMatrix"))
          APP_ABORT("Error reading RotationMatrix.\n");
        if (std::get<1>(R.sizes()) != NMO)
          APP_ABORT("Error Wrong dimensions in RotationMatrix.\n");
        dim[0] = R.size();
        dim[1] = 0;
        // conjugate rotation matrix
        std::transform(R.origin(), R.origin() + R.num_elements(), R.origin(),
                       [](const auto& c) { return std::conj(c); });
        TG.Node().broadcast_n(dim, 2, 0);
        XRot = sharedCMatrix({dim[0], NMO}, make_node_allocator<ComplexType>(TG));
        copy_n(R.origin(), R.num_elements(), make_device_ptr(XRot.origin()));
        if (TG.Node().root())
          TG.Cores().broadcast_n(to_address(XRot.origin()), XRot.num_elements(), 0);

        dump.pop();
        dump.close();
      }
      else
      {
        TG.Node().broadcast_n(dim, 2, 0);
        XRot = sharedCMatrix({dim[0], NMO}, make_node_allocator<ComplexType>(TG));
        if (TG.Node().root())
          TG.Cores().broadcast_n(to_address(XRot.origin()), XRot.num_elements(), 0);
      }
      TG.Node().barrier();

      dm_size = XRot.size() * XRot.size() * XRot.size() * XRot.size();
    }
    else
    {
      dm_size = NMO * NMO * NMO * NMO;
    }

    // (a,a,a,a), (a,a,b,b)
    nspinblocks = 2;
    if (walker_type == COLLINEAR)
    {
      nspinblocks = 3; // (a,a,a,a), (a,a,b,b), (b,b,b,b)
    }
    else if (walker_type == NONCOLLINEAR)
      APP_ABORT(" Error: NONCOLLINEAR not yet implemented. \n\n\n");

    dm_size *= nspinblocks;

    using std::fill_n;
    writer = (TG.getGlobalRank() == 0);

    DMAverage = mpi3CMatrix({nave, dm_size}, shared_allocator<ComplexType>{TG.TG_local()});
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

    // check structure dimensions
    if (iref == 0)
    {
      if (denom.size() != nw)
      {
        denom = mpi3CVector(iextensions<1u>{nw}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(DMWork.sizes()) != nw || std::get<1>(DMWork.sizes()) != dm_size)
      {
        DMWork = mpi3CMatrix({nw, dm_size}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      fill_n(denom.origin(), denom.num_elements(), ComplexType(0.0, 0.0));
      fill_n(DMWork.origin(), DMWork.num_elements(), ComplexType(0.0, 0.0));
    }
    else
    {
      if (std::get<0>(denom.sizes()) != nw || std::get<0>(DMWork.sizes()) != nw || std::get<1>(DMWork.sizes()) != dm_size || std::get<0>(DMAverage.sizes()) != nave ||
          std::get<1>(DMAverage.sizes()) != dm_size)
        APP_ABORT(" Error: Invalid state in accumulate_reference. \n\n\n");
    }

    if (apply_rotation)
      acc_with_rotation(G, Xw);
    else
      acc_no_rotation(G, Xw);
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
    ma::product(ComplexType(1.0, 0.0), ma::T(DMWork({0, nw}, {i0, iN})), denom, ComplexType(1.0, 0.0),
                DMAverage[iav].sliced(i0, iN));
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
        dump.push(std::string("FullTwoRDM"));
        for (int i = 0; i < nave; ++i)
        {
          dump.push(std::string("Average_") + std::to_string(i));
          std::string padded_iblock =
              std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
          stdCVector_ref DMAverage_(to_address(DMAverage[i].origin()), {dm_size});
          dump.write(DMAverage_, "two_rdm_" + padded_iblock);
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

  bool writer;

  int block_size;

  int nave;

  int counter;

  int nspinblocks;

  int dm_size;

  bool apply_rotation;

  sharedCMatrix XRot;
  stdCVector Grot;

  mpi3CVector denom;

  // DMAverage (nave, nspinblocks, x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMAverage;
  // DMWork (nwalk, nspinblocks, x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMWork;

  // expects host accesible memory
  template<class MatG, class CVec>
  void acc_no_rotation(MatG&& G, CVec&& Xw)
  {
    // doing this 1 walker at a time and not worrying about speed
    int nw(G.size());

    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), NMO * NMO, TG.TG_local().size());
    int dN           = iN - i0;

    size_t M2(NMO * NMO);
    size_t M4(M2 * M2);
    size_t N = size_t(dN) * M2;
    DeviceBufferManager buffer_manager;
    StaticMatrix R({dN, NMO * NMO}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    CMatrix_ref Q(R.origin(), {NMO * NMO, NMO});
    CVector_ref R1D(R.origin(), R.num_elements());
    CVector_ref Q1D(Q.origin(), Q.num_elements());

    // put this in shared memory!!!
    StaticMatrix Gt({NMO, NMO}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    CMatrix_ref GtC(Gt.origin(), {NMO * NMO, 1});
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
    if (Grot.size() < R.num_elements())
      Grot = stdCVector(iextensions<1u>(R.num_elements()));
#endif

    if (TG.TG_local().size() > 1)
      APP_ABORT(" ncores > 1 not yet allowed in full2rdm.\n\n\n");

    // (ikjl) = Gik * Gjl - (same spin) Gil Gjk
    if (walker_type == COLLINEAR)
    {
      for (int iw = 0; iw < nw; iw++)
      {
        if (TG.TG_local().root())
          denom[iw] += Xw[iw];

        CMatrix_ref Gup(make_device_ptr(G[iw][0].origin()), {NMO * NMO, 1});
        CMatrix_ref Gdn(make_device_ptr(G[iw][1].origin()), {NMO * NMO, 1});
        // use ger !!!!

        //  (a,a,a,a)
        ma::product(Gup.sliced(i0, iN), ma::T(Gup), R);
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        using std::copy_n;
        copy_n(R.origin(), R.num_elements(), Grot.origin());
        ma::axpy(Xw[iw], Grot, DMWork[iw].sliced(size_t(i0) * M2, size_t(iN) * M2));
#else
        ma::axpy(Xw[iw], R1D, DMWork[iw].sliced(size_t(i0) * M2, size_t(iN) * M2));
#endif
        TG.TG_local().barrier();
        ma::transpose(G[iw][0], Gt);
        // EXX: (ikjl) -= Gil Gjk = Gt_kj  x G[i]l  for each i
        // parallelize this!!!
        for (int i = 0; i < NMO; ++i)
        {
          ma::product(ComplexType(-1.0), GtC, G[iw][0].sliced(i, i + 1), ComplexType(0.0), Q);
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
          using std::copy_n;
          copy_n(Q.origin(), Q.num_elements(), Grot.origin());
          ma::axpy(Xw[iw], Grot.sliced(0, Q.num_elements()),
                   DMWork[iw].sliced(size_t(i * NMO) * M2, size_t((i + 1) * NMO) * M2));
#else
          ma::axpy(Xw[iw], Q1D, DMWork[iw].sliced(size_t(i * NMO) * M2, size_t((i + 1) * NMO) * M2));
#endif
        }

        //  (a,a,b,b)
        ma::product(Gup.sliced(i0, iN), ma::T(Gdn), R);
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        using std::copy_n;
        copy_n(R.origin(), R.num_elements(), Grot.origin());
        ma::axpy(Xw[iw], Grot, DMWork[iw].sliced(M4 + size_t(i0) * M2, M4 + size_t(iN) * M2));
#else
        ma::axpy(Xw[iw], R1D, DMWork[iw].sliced(M4 + size_t(i0) * M2, M4 + size_t(iN) * M2));
#endif

        //  (b,b,b,b)
        ma::product(Gdn.sliced(i0, iN), ma::T(Gdn), R);
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        using std::copy_n;
        copy_n(R.origin(), R.num_elements(), Grot.origin());
        ma::axpy(Xw[iw], Grot, DMWork[iw].sliced(2 * M4 + size_t(i0) * M2, 2 * M4 + size_t(iN) * M2));
#else
        ma::axpy(Xw[iw], R1D, DMWork[iw].sliced(2 * M4 + size_t(i0) * M2, 2 * M4 + size_t(iN) * M2));
#endif
        TG.TG_local().barrier();
        ma::transpose(G[iw][1], Gt);
        // EXX: (ikjl) -= Gil Gjk = Gt_kj  x G[i]l  for each i
        // parallelize this!!!
        for (int i = 0; i < NMO; ++i)
        {
          ma::product(ComplexType(-1.0), GtC, G[iw][1].sliced(i, i + 1), ComplexType(0.0), Q);
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
          using std::copy_n;
          copy_n(Q.origin(), Q.num_elements(), Grot.origin());
          ma::axpy(Xw[iw], Grot.sliced(0, Q.num_elements()),
                   DMWork[iw].sliced(2 * M4 + size_t(i * NMO) * M2, 2 * M4 + size_t((i + 1) * NMO) * M2));
#else
          ma::axpy(Xw[iw], Q1D, DMWork[iw].sliced(2 * M4 + size_t(i * NMO) * M2, 2 * M4 + size_t((i + 1) * NMO) * M2));
#endif
        }
      }
    }
    else
    {
      /*
      int N = dN*NMO*NMO;
      set_buffer(N);
      CMatrix_ref R( Buff.origin(), {dN,NMO*NMO});
      CVector_ref R1D( Buff.origin(), {dN*NMO*NMO});
#if defined(ENABLE_CUDA)
      if(Grot.size() < R.num_elements()) 
        Grot = stdCVector(iextensions<1u>(R.num_elements()));
#endif
        

      for(int iw=0; iw<nw; iw++) {

        if(TG.TG_local().root()) denom[iw] += Xw[iw];

        CMatrix_ref Gw( to_address(G[iw].origin()), {G[0].num_elements(),1});
        // use ger later
        ma::product( Gw.sliced(i0,iN), ma::T(Gw), R );

#if defined(ENABLE_CUDA)
        using std::copy_n;
        copy_n(R.origin(),R.num_elements(),Grot.origin());
        ma::axpy( Xw[iw], Grot, DMWork[iw].sliced(i0*NMO*NMO,iN*NMO*NMO) );
#else
        ma::axpy( Xw[iw], R.sliced(i0,iN), DMWork[iw].sliced(i0,iN) );
#endif
      }
*/
    }
    TG.TG_local().barrier();
  }

  // G should be device accesible memory
  // Xw is in host
  template<class MatG, class CVec>
  void acc_with_rotation(MatG&& G, CVec&& Xw)
  {
    /*
    int i0, iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),int(XRot.size(0)),TG.TG_local().size());

    // can batch in the future if too slow
    // Grot = Xc * G * H(Xc)
    int nX = XRot.size(0);
    int sz = nX * (NMO + (iN-i0)); 
    int npts = (iN-i0)*nX;
    set_buffer(sz);
    CMatrix_ref T1(Buff.origin(),{(iN-i0),NMO});
    CMatrix_ref T2(T1.origin()+T1.num_elements(),{(iN-i0),nX});
    if(Grot.size() != npts) 
      Grot = stdCVector(iextensions<1u>(npts));

    // round-robin for now
    int cnt=0;
    for(int iw=0; iw<nw; iw++) {
      if(i0==iN || i0==XRot.size(0)) break;
      if(TG.TG_local().root()) denom[iw] += Xw[iw];
      ma::product(XRot.sliced(i0,iN),G[iw][0],T1);
      ma::product(T1,ma::H(XRot),T2);
      copy_n(T2.origin(),T2.num_elements(),Grot.origin());      
      ma::axpy( Xw[iw], Grot, DMWork[iw].sliced(i0*nX,i0*nX+npts) );
      if(walker_type == COLLINEAR) {
        ma::product(XRot.sliced(i0,iN),G[iw][1],T1);
        ma::product(T1,ma::H(XRot),T2);
        copy_n(T2.origin(),T2.num_elements(),Grot.origin());
        ma::axpy( Xw[iw], Grot, DMWork[iw].sliced((nX+i0)*nX,(nX+i0)*nX+npts) );
      }
    }
    TG.TG_local().barrier();
*/
  }
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
