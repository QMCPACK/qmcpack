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

#ifndef QMCPLUSPLUS_AFQMC_ATOMCENTERED_CORRELATORS_HPP
#define QMCPLUSPLUS_AFQMC_ATOMCENTERED_CORRELATORS_HPP

#include "AFQMC/config.h"
#include <vector>
#include <string>
#include <iostream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/Timer.h"

#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/detail/utilities.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
/* 
 * Observable class that calculates the walker averaged on-top pair density 
 * Alloc defines the allocator type used to store the orbital and temporary tensors.
 * In a device compilation, this would need out-of-card gemm available.
 */
class atomcentered_correlators : public AFQMCInfo
{
  // allocators
  using Allocator  = localTG_allocator<ComplexType>;
  using IAllocator = device_allocator<int>;

  // type defs
  using pointer     = device_ptr<ComplexType>;
  using int_pointer = device_ptr<int>;

  using CVector_ref    = boost::multi::array_ref<ComplexType, 1, pointer>;
  using CMatrix_ref    = boost::multi::array_ref<ComplexType, 2, pointer>;
  using CTensor_ref    = boost::multi::array_ref<ComplexType, 3, pointer>;
  using C4Tensor_ref   = boost::multi::array_ref<ComplexType, 4, pointer>;
  using C5Tensor_ref   = boost::multi::array_ref<ComplexType, 5, pointer>;
  using IVector        = boost::multi::array<int, 1, IAllocator>;
  using CVector        = boost::multi::array<ComplexType, 1, Allocator>;
  using CMatrix        = boost::multi::array<ComplexType, 2, Allocator>;
  using CTensor        = boost::multi::array<ComplexType, 3, Allocator>;
  using C4Tensor       = boost::multi::array<ComplexType, 4, Allocator>;
  using stdCVector     = boost::multi::array<ComplexType, 1>;
  using stdCMatrix     = boost::multi::array<ComplexType, 2>;
  using stdIVector     = boost::multi::array<int, 1>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType, 1>;
  using stdCMatrix_ref = boost::multi::array_ref<ComplexType, 2>;
  using mpi3CVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using mpi3CMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using mpi3CTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using mpi3C4Tensor   = boost::multi::array<ComplexType, 4, shared_allocator<ComplexType>>;

  using shm_stack_alloc_type = LocalTGBufferManager::template allocator_t<ComplexType>;
  using StaticMatrix         = boost::multi::static_array<ComplexType, 2, shm_stack_alloc_type>;
  using Static3Tensor        = boost::multi::static_array<ComplexType, 3, shm_stack_alloc_type>;

  // MAM: Note -
  // This class uses lots of memory, but can be safely moved to single precision.
  // Make S, XY and all temporary calculations on the device in single precision,
  // then copy_n_cast to host in double precision and accumulate also in DP.
  // NwIJ, NwIJ, denom, DMWork, DMAverage in DP
  // S, XY, Buff in SP
  // use in_place_cast from DP to SP and back for input G

public:
  atomcentered_correlators(afqmc::TaskGroup_& tg_,
                           AFQMCInfo& info,
                           xmlNodePtr cur,
                           WALKER_TYPES wlk,
                           int nave_ = 1,
                           int bsize = 1)
      : AFQMCInfo(info),
        alloc(make_localTG_allocator<ComplexType>(tg_)),
        block_size(bsize),
        nave(nave_),
        counter(0),
        TG(tg_),
        walker_type(wlk),
        ns2(0),
        writer(false),
        S({0, 0, 0}, make_node_allocator<ComplexType>(TG)),
        XY({0, 0}, make_node_allocator<ComplexType>(TG)),
        shapes(iextensions<1u>{0}, IAllocator{}),
        DMAverage2D({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork2D({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMAverage1D({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork1D({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        denom(iextensions<1u>{0}, shared_allocator<ComplexType>{TG.TG_local()}),
        NwIJ({0, 0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        NwI({0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    app_log() << "  --  Adding Back Propagated atom-centered correlation functions (atomcentered_correlators). -- \n ";
    std::string orb_file("");
    std::string str("false");

    if (cur != NULL)
    {
      ParameterSet m_param;
      m_param.add(orb_file, "filename");
      m_param.put(cur);
    }

    if (not file_exists(orb_file))
    {
      app_error() << " Error: orbitals file does not exist: " << orb_file << std::endl;
      APP_ABORT("");
    }


    // read orbitals
    hdf_archive dump;
    if (TG.Node().root())
    {
      if (!dump.open(orb_file, H5F_ACC_RDONLY))
      {
        app_error() << " Error opening orbitals file for atomcentered_correlators estimator. \n";
        APP_ABORT("");
      }

      dump.push("AtomicOverlaps", false);

      if (!dump.readEntry(orbs, "index_of_orbital_per_site"))
      {
        app_error() << " Error in atomcentered_correlators: Problems reading index_of_orbital_per_site. " << std::endl;
        APP_ABORT("");
      }
      nsites = orbs.size() - 1;
      if (nsites < 1)
      {
        app_error() << " Error in atomcentered_correlators: nsites < 1. " << std::endl;
        APP_ABORT("");
      }
      if (orbs[0] != 0)
      {
        app_error() << " Error in atomcentered_correlators: orbs[0] != 0. " << std::endl;
        APP_ABORT("");
      }
      for (int i = 0; i < nsites; i++)
        if (orbs[i] > orbs[i + 1])
        {
          app_error() << " Error in atomcentered_correlators: orbs[i] > orbs[i+1]. " << std::endl;
          APP_ABORT("");
        }
      TG.Node().broadcast_n(&nsites, 1, 0);
      TG.Node().broadcast_n(orbs.begin(), orbs.size(), 0);
      NAO = orbs.back();
      S   = CTensor({2, NAO, NMO}, make_node_allocator<ComplexType>(TG));
      XY  = CMatrix({NAO, NAO}, make_node_allocator<ComplexType>(TG));
      CMatrix_ref S_(make_device_ptr(S[0].origin()), {NAO, NMO});
      if (!dump.readEntry(S_, "Left"))
      {
        app_error() << " Error in atomcentered_correlators: Problems reading Left. " << std::endl;
        APP_ABORT("");
      }
      stdCMatrix_ref S2_(to_address(S[1].origin()), {NAO, NMO});
      if (!dump.readEntry(S2_, "Right"))
      {
        app_error() << " Error in atomcentered_correlators: Problems reading Right. " << std::endl;
        APP_ABORT("");
      }
      dump.pop();
      dump.close();
      // use faster approach
      // slowwwww!!!!!
      for (int i = 0; i < NAO; i++)
        for (int j = 0; j < NAO; j++)
          S[0][i][j] = std::conj(ComplexType(S[0][i][j]));
      ma::product(S[0], ma::T(S[1]), XY);
    }
    else
    {
      TG.Node().broadcast_n(&nsites, 1, 0);
      orbs.resize(nsites + 1);
      TG.Node().broadcast_n(orbs.begin(), orbs.size(), 0);
      NAO = orbs.back();
      S   = CTensor({2, NAO, NMO}, make_node_allocator<ComplexType>(TG));
      XY  = CMatrix({NAO, NAO}, make_node_allocator<ComplexType>(TG));
    }
    ns2 = nsites * nsites;
    TG.Node().barrier();
    app_log() << " Number of atomic sites: " << nsites << std::endl;
    app_log() << " Number of (pseudo) localized atomic orbitals: " << NAO << std::endl;

    using std::fill_n;
    writer = (TG.getGlobalRank() == 0);

    if (writer)
    {
      type_id1D.reserve(3);
      type_id1D.emplace_back("C");
      type_id1D.emplace_back("S");
      type_id1D.emplace_back("M");
      type_id2D.reserve(3);
      type_id2D.emplace_back("CC");
      type_id2D.emplace_back("SS");
      type_id2D.emplace_back("CS");
    }

    DMAverage1D = mpi3CTensor({nave, 3, nsites}, shared_allocator<ComplexType>{TG.TG_local()});
    fill_n(DMAverage1D.origin(), DMAverage1D.num_elements(), ComplexType(0.0, 0.0));
    DMAverage2D = mpi3CTensor({nave, 3, ns2}, shared_allocator<ComplexType>{TG.TG_local()});
    fill_n(DMAverage2D.origin(), DMAverage2D.num_elements(), ComplexType(0.0, 0.0));
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
    using ma::gemmBatched;
    using std::copy_n;
    using std::fill_n;
    // assumes G[nwalk][spin][M][M]
    int nw(G.size());
    assert(G.size() == wgt.size());
    assert(wgt.size() == nw);
    assert(Xw.size() == nw);
    assert(ovlp.size() >= nw);
    assert(G.num_elements() == G_host.num_elements());
    assert(G.extensions() == G_host.extensions());

    int nsp;
    if (walker_type == CLOSED)
      nsp = 1;
    else
      nsp = 2;

    // check structure dimensions
    if (iref == 0)
    {
      if (denom.size() != nw)
      {
        denom = mpi3CVector(iextensions<1u>{nw}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(DMWork1D.sizes()) != nw || std::get<1>(DMWork1D.sizes()) != 3 || std::get<2>(DMWork1D.sizes()) != nsites)
      {
        DMWork1D = mpi3CTensor({nw, 3, nsites}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(DMWork2D.sizes()) != nw || std::get<1>(DMWork2D.sizes()) != 3 || std::get<2>(DMWork2D.sizes()) != ns2)
      {
        DMWork2D = mpi3CTensor({nw, 3, ns2}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(NwIJ.sizes()) != nsp || std::get<1>(NwIJ.sizes()) != nw || std::get<2>(NwIJ.sizes()) != nsites || std::get<3>(NwIJ.sizes()) != nsites)
      {
        NwIJ = mpi3C4Tensor({nsp, nw, nsites, nsites}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (std::get<0>(NwI.sizes()) != nsp || std::get<1>(NwI.sizes()) != nw || std::get<2>(NwI.sizes()) != nsites)
      {
        NwI = mpi3CTensor({nsp, nw, nsites}, shared_allocator<ComplexType>{TG.TG_local()});
      }
      if (shapes.size() < 2 * nw * nsites * nsites)
        shapes = IVector(iextensions<1u>{2 * nw * nsites * nsites}, IAllocator{});
      fill_n(denom.origin(), denom.num_elements(), ComplexType(0.0, 0.0));
      fill_n(DMWork1D.origin(), DMWork1D.num_elements(), ComplexType(0.0, 0.0));
      fill_n(DMWork2D.origin(), DMWork2D.num_elements(), ComplexType(0.0, 0.0));
    }
    else
    {
      if (std::get<0>(denom.sizes()) != nw || std::get<0>(DMWork1D.sizes()) != nw || std::get<1>(DMWork1D.sizes()) != 3 || std::get<2>(DMWork1D.sizes()) != nsites ||
          std::get<0>(DMWork2D.sizes()) != nw || std::get<1>(DMWork2D.sizes()) != 3 || std::get<2>(DMWork2D.sizes()) != ns2 || std::get<0>(NwI.sizes()) != nsp ||
          std::get<1>(NwI.sizes()) != nw || std::get<2>(NwI.sizes()) != nsites || std::get<0>(NwIJ.sizes()) != nsp || std::get<1>(NwIJ.sizes()) != nw ||
          std::get<2>(NwIJ.sizes()) != nsites || std::get<3>(NwIJ.sizes()) != nsites || std::get<0>(DMAverage1D.sizes()) != nave || std::get<1>(DMAverage1D.sizes()) != 3 ||
          std::get<2>(DMAverage1D.sizes()) != nsites || std::get<0>(DMAverage2D.sizes()) != nave || std::get<1>(DMAverage2D.sizes()) != 3 ||
          std::get<2>(DMAverage2D.sizes()) != ns2)
        APP_ABORT(" Error: Invalid state in accumulate_reference. \n\n\n");
    }

    // calculate green functions in atom basis and send to host for processing/accumulation
    // if memory becomes a problem, then batch over walkers
    // calculate nwbatch based on size of available memory on global buffer
    int iw0(0);
    int i0, iN;
    int nwbatch(nw);
    std::vector<decltype(ma::pointer_dispatch(S.origin()))> Aarray;
    Aarray.reserve(nwbatch * nsites * nsites);
    std::vector<decltype(ma::pointer_dispatch(S.origin()))> Barray;
    Barray.reserve(nwbatch * nsites * nsites);
    std::vector<decltype(ma::pointer_dispatch(S.origin()))> Carray;
    Carray.reserve(nwbatch * nsites * nsites);
    LocalTGBufferManager buffer_manager;
    while (iw0 < nw)
    {
      int nwlk = std::min(nwbatch, nw - iw0);

      Static3Tensor QwI({nwlk, NAO, NMO}, buffer_manager.get_generator().template get_allocator<ComplexType>());
      Static3Tensor MwIJ({nwlk, NAO, NAO}, buffer_manager.get_generator().template get_allocator<ComplexType>());
      Static3Tensor devNwIJ({nwlk, nsites, nsites},
                            buffer_manager.get_generator().template get_allocator<ComplexType>());
      StaticMatrix devNwI({nwlk, nsites}, buffer_manager.get_generator().template get_allocator<ComplexType>());

      for (int is = 0; is < nsp; ++is)
      {
        std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), int(NAO), TG.TG_local().size());
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        // QwI[iw] = X * G[iw][is] =  S[0] * G[iw][is]
        for (int iw = 0; iw < nwlk; ++iw)
        {
          Aarray.emplace_back(ma::pointer_dispatch(S[0][i0].origin()));
          Barray.emplace_back(ma::pointer_dispatch(G[iw0 + iw][is].origin()));
          Carray.emplace_back(ma::pointer_dispatch(QwI[iw][i0].origin()));
        }
        // careful with fortran ordering
        gemmBatched('N', 'N', NMO, int(iN - i0), NMO, ComplexType(1.0), Barray.data(), NMO, Aarray.data(), NMO,
                    ComplexType(0.0), Carray.data(), NMO, int(Aarray.size()));

        // no need to wait, since I'm working on same segment I'm producing
        // Carray already contains QwI
        Aarray.clear();
        Barray.clear();
        // MwIJ[iw] = QwI * Y =  QwI * T(S[1])
        for (int iw = 0; iw < nwlk; ++iw)
        {
          Barray.emplace_back(ma::pointer_dispatch(S[1].origin()));
          Aarray.emplace_back(ma::pointer_dispatch(MwIJ[iw][i0].origin()));
        }
        // careful with fortran ordering
        gemmBatched('T', 'N', NAO, int(iN - i0), NMO, ComplexType(1.0), Barray.data(), NMO, Carray.data(), NMO,
                    ComplexType(0.0), Aarray.data(), NAO, int(Aarray.size()));
        TG.TG_local().barrier(); // now wait!

        // NwI[iw][I] = sum_p MwIJ[iw][Ip][Ip]
        Aarray.clear();
        Barray.clear();
        for (int iw = 0, p = 0, nt = 0; iw < nwlk; ++iw)
          for (int I = 0; I < nsites; ++I, ++p)
            if (p % TG.TG_local().size() == TG.TG_local().rank())
            {
              Aarray.emplace_back(ma::pointer_dispatch(devNwI[iw].origin() + I));
              Barray.emplace_back(ma::pointer_dispatch(MwIJ[iw][orbs[I]].origin() + orbs[I]));
              shapes[nt++] = int(orbs[I + 1] - orbs[I]);
            }
        using ma::batched_diagonal_sum;
        batched_diagonal_sum(shapes.data(), Barray.data(), NAO, ComplexType(1.0), Aarray.data(), int(Aarray.size()));

        fill_n(devNwIJ.origin(), devNwIJ.num_elements(), ComplexType(0.0));

        // NwIJ[iw][I][J] = sum_ij MwIJ[iw][Ij][Jj] * XY[Jj][Ii]
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int iw = 0, p = 0, nt = 0; iw < nwlk; ++iw)
          for (int I = 0; I < nsites; ++I)
            for (int J = 0; J < nsites; ++J, ++p)
              if (p % TG.TG_local().size() == TG.TG_local().rank())
              {
                Aarray.emplace_back(ma::pointer_dispatch(MwIJ[iw][orbs[I]].origin() + orbs[J]));
                Barray.emplace_back(ma::pointer_dispatch(XY[orbs[J]].origin() + orbs[I]));
                Carray.emplace_back(ma::pointer_dispatch(devNwIJ[iw][I].origin() + J));
                shapes[nt++] = int(orbs[I + 1] - orbs[I]);
                shapes[nt++] = int(orbs[J + 1] - orbs[J]);
              }
        using ma::batched_ab_ba;
        batched_ab_ba(shapes.data(), Aarray.data(), NAO, Barray.data(), NAO, ComplexType(1.0), Carray.data(),
                      int(Aarray.size()));
        TG.TG_local().barrier();

        // testing !!!
        //        fill_n(devNwIJ.origin(),devNwIJ.num_elements(),ComplexType(0.0));
        //        TG.TG_local().barrier();

        // NwIJ[iw][I][J] = sum_ij MwIJ[iw][Ij][Jj] * MwJI[iw][Jj][Ii]
        Aarray.clear();
        Barray.clear();
        Carray.clear();
        for (int iw = 0, p = 0, nt = 0; iw < nwlk; ++iw)
          for (int I = 0; I < nsites; ++I)
            for (int J = 0; J < nsites; ++J, ++p)
              if (p % TG.TG_local().size() == TG.TG_local().rank())
              {
                Aarray.emplace_back(ma::pointer_dispatch(MwIJ[iw][orbs[I]].origin() + orbs[J]));
                Barray.emplace_back(ma::pointer_dispatch(MwIJ[iw][orbs[J]].origin() + orbs[I]));
                Carray.emplace_back(ma::pointer_dispatch(devNwIJ[iw][I].origin() + J));
                shapes[nt++] = int(orbs[I + 1] - orbs[I]);
                shapes[nt++] = int(orbs[J + 1] - orbs[J]);
              }
        using ma::batched_ab_ba;
        batched_ab_ba(shapes.data(), Aarray.data(), NAO, Barray.data(), NAO, ComplexType(-1.0), Carray.data(),
                      int(Aarray.size()));
        TG.TG_local().barrier();

        std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), int(devNwIJ.num_elements()), TG.TG_local().size());
        copy_n(devNwIJ.origin() + i0, (iN - i0), NwIJ[is][iw0].origin() + i0);
        std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), int(devNwI.num_elements()), TG.TG_local().size());
        copy_n(devNwI.origin() + i0, (iN - i0), NwI[is][iw0].origin() + i0);
      }
      iw0 += nwlk;
    }
    TG.TG_local().barrier();

    // : Z(type,r,r') = < (nup_r OP1 ndn_r) (nup_r' OP2 ndn_r') >
    // where type= 0: charge-charge --> ++
    //             1: spin-spin     --> --
    //             2: charge-spin   --> +-
    // THIS COULD BE WRITTEN WITH OPENMP TO ENABLE GPU EASILY!!!
    for (int iw = 0; iw < nw; iw++)
    {
      if (TG.TG_local().root())
        denom[iw] += Xw[iw];
      if (iw % TG.TG_local().size() != TG.TG_local().rank())
        continue;
      if (walker_type == CLOSED)
      {
        stdCMatrix_ref DMcc(to_address(DMWork2D[iw][0].origin()), {nsites, nsites});
        stdCMatrix_ref DMss(to_address(DMWork2D[iw][1].origin()), {nsites, nsites});
        stdCMatrix_ref DMcs(to_address(DMWork2D[iw][2].origin()), {nsites, nsites});
        auto X_(2.0 * Xw[iw]);
        for (int i = 0; i < nsites; i++)
        {
          DMWork1D[iw][0][i] += X_ * NwI[0][iw][i];
          for (int j = 0; j < nsites; j++)
          {
            auto dJuu(NwI[0][iw][i] * NwI[0][iw][j]);
            auto dXuu(NwIJ[0][iw][i][j]);
            DMcc[i][j] += X_ * (2.0 * dJuu + dXuu);
            DMss[i][j] += X_ * (dXuu);
            // no cs in RHF
          }
        }
      }
      else if (walker_type == COLLINEAR)
      {
        stdCMatrix_ref DMcc(to_address(DMWork2D[iw][0].origin()), {nsites, nsites});
        stdCMatrix_ref DMss(to_address(DMWork2D[iw][1].origin()), {nsites, nsites});
        stdCMatrix_ref DMcs(to_address(DMWork2D[iw][2].origin()), {nsites, nsites});
        auto X_(Xw[iw]);
        for (int i = 0; i < nsites; i++)
        {
          DMWork1D[iw][0][i] += X_ * (NwI[0][iw][i] + NwI[1][iw][i]);
          DMWork1D[iw][1][i] += X_ * (NwI[0][iw][i] - NwI[1][iw][i]);
          DMWork1D[iw][2][i] += std::real(ComplexType(X_ * (NwI[0][iw][i] - NwI[1][iw][i])));
          for (int j = 0; j < nsites; j++)
          {
            auto uu(NwI[0][iw][i] * NwI[0][iw][j] + NwIJ[0][iw][i][j]);
            auto dd(NwI[1][iw][i] * NwI[1][iw][j] + NwIJ[1][iw][i][j]);
            auto ud(NwI[0][iw][i] * NwI[1][iw][j]);
            auto du(NwI[1][iw][i] * NwI[0][iw][j]);
            DMcc[i][j] += X_ * (uu + dd + ud + du);
            DMss[i][j] += X_ * (uu + dd - ud - du);
            DMcs[i][j] += X_ * (uu - dd - ud + du);
            // no cs in RHF
          }
        }
      }
      else if (walker_type == NONCOLLINEAR)
      {
        APP_ABORT(" Noncollinear not implemented \n\n\n");
      }
    }
    TG.TG_local().barrier();
  }

  template<class HostCVec>
  void accumulate_block(int iav, HostCVec&& wgt, bool impsamp)
  {
    int nw(denom.size());
    TG.TG_local().barrier();
    // this is meant to be small, so serializing
    if (TG.TG_local().root())
    {
      for (int iw = 0; iw < nw; iw++)
        denom[iw] = wgt[iw] / denom[iw];
      {
        stdCVector_ref DMAv1D(to_address(DMAverage1D[iav].origin()), {3 * nsites});
        stdCMatrix_ref DMWork2D_(to_address(DMWork1D.origin()), {nw, 3 * nsites});
        // DMAverage[iav][t][ij] += sum_iw DMWork[iw][t][ij] * denom[iw] = T( DMWork ) * denom
        ma::product(ComplexType(1.0, 0.0), ma::T(DMWork2D_), denom, ComplexType(1.0, 0.0), DMAv1D);
      }
      {
        stdCVector_ref DMAv1D(to_address(DMAverage2D[iav].origin()), {3 * ns2});
        stdCMatrix_ref DMWork2D_(to_address(DMWork2D.origin()), {nw, 3 * ns2});
        // DMAverage[iav][t][ij] += sum_iw DMWork[iw][t][ij] * denom[iw] = T( DMWork ) * denom
        ma::product(ComplexType(1.0, 0.0), ma::T(DMWork2D_), denom, ComplexType(1.0, 0.0), DMAv1D);
      }
    }
    TG.TG_local().barrier();
  }

  template<class HostCVec>
  void print(int iblock, hdf_archive& dump, HostCVec&& Wsum)
  {
    using std::fill_n;
    const int n_zero = 9;

    if (TG.TG_local().root())
    {
      ma::scal(ComplexType(1.0 / block_size), DMAverage1D);
      TG.TG_heads().reduce_in_place_n(to_address(DMAverage1D.origin()), DMAverage1D.num_elements(), std::plus<>(), 0);
      ma::scal(ComplexType(1.0 / block_size), DMAverage2D);
      TG.TG_heads().reduce_in_place_n(to_address(DMAverage2D.origin()), DMAverage2D.num_elements(), std::plus<>(), 0);
      if (writer)
      {
        dump.push(std::string("ATOM_CORRELATORS"));
        for (int t = 0; t < 3; ++t)
        {
          dump.push(type_id2D[t]);
          for (int i = 0; i < nave; ++i)
          {
            dump.push(std::string("Average_") + std::to_string(i));
            std::string padded_iblock =
                std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
            stdCVector_ref DMAverage2D_(to_address(DMAverage2D[i][t].origin()), {ns2});
            dump.write(DMAverage2D_, "correlator2D_" + type_id2D[t] + padded_iblock);
            dump.write(Wsum[i], "denominator_" + padded_iblock);
            dump.pop();
          }
          dump.pop();
        }
        for (int t = 0; t < 3; ++t)
        {
          dump.push(type_id1D[t]);
          for (int i = 0; i < nave; ++i)
          {
            dump.push(std::string("Average_") + std::to_string(i));
            std::string padded_iblock =
                std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
            stdCVector_ref DMAverage1D_(to_address(DMAverage1D[i][t].origin()), {nsites});
            dump.write(DMAverage1D_, "correlator1D_" + type_id1D[t] + padded_iblock);
            dump.write(Wsum[i], "denominator_" + padded_iblock);
            dump.pop();
          }
          dump.pop();
        }
        dump.pop();
      }
    }
    TG.TG_local().barrier();
    fill_n(DMAverage1D.origin(), DMAverage1D.num_elements(), ComplexType(0.0, 0.0));
    fill_n(DMAverage2D.origin(), DMAverage2D.num_elements(), ComplexType(0.0, 0.0));
  }

private:
  Allocator alloc;

  int block_size;

  int nave;

  int counter;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  int ns2;

  bool writer;

  std::vector<std::string> type_id1D;
  std::vector<std::string> type_id2D;

  int NAO;
  int nsites;
  std::vector<int> orbs;

  // S_IJ = H(X_I) * Y_J,
  // S[0] = X // stored as conj(X)
  // S[1] = Y
  CTensor S;
  //XY  = conj(X) * T(Y)
  CMatrix XY;
  IVector shapes;

  // type:0-cc, 1-ss, 2-cs (c=charge, s=spin)
  // np = # of atomic sites
  // DMAverage (nave, type, np*np )
  mpi3CTensor DMAverage2D;
  // DMWork (nwalk, type, np*np )
  mpi3CTensor DMWork2D;
  // DMAverage (nave, type, np )
  mpi3CTensor DMAverage1D;
  // DMWork (nwalk, type, np )
  mpi3CTensor DMWork1D;

  mpi3CVector denom;

  mpi3C4Tensor NwIJ;
  mpi3CTensor NwI;
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
