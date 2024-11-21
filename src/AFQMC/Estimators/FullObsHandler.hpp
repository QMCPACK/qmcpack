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

#ifndef QMCPLUSPLUS_AFQMC_FULLOBSHANDLER_HPP
#define QMCPLUSPLUS_AFQMC_FULLOBSHANDLER_HPP

#include <vector>
#include <string>
#include <iostream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Estimators/Observables/Observable.hpp"
#include "AFQMC/config.h"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
/*
 * This class manages a list of "full" observables.
 * Full observables are those that have a walker dependent left-hand side, 
 * which result from back propagation. 
 * This implementation of the class assumes a multi-determinant trial wavefunction,  
 * resulting in the loop over "references" (each determinant in the trial wavefunciton
 * being back-propagated). 
 * Given a walker set and an array of (back propagated) slater matrices, 
 * this routine will calculate and accumulate all requested observables.
 * To make the implementation of the BackPropagated class cleaner, 
 * this class also handles all the hdf5 I/O (given a hdf archive).
 */
class FullObsHandler : public AFQMCInfo
{
  // allocators
  using sharedAllocator = localTG_allocator<ComplexType>;

  using shared_pointer       = typename sharedAllocator::pointer;
  using const_shared_pointer = typename sharedAllocator::const_pointer;

  using devCMatrix_ptr = boost::multi::array_ptr<ComplexType, 2, device_ptr<ComplexType>>;

  using sharedCVector      = boost::multi::array<ComplexType, 1, sharedAllocator>;
  using sharedCVector_ref  = boost::multi::array_ref<ComplexType, 1, shared_pointer>;
  using sharedCMatrix_ref  = boost::multi::array_ref<ComplexType, 2, shared_pointer>;
  using sharedC4Tensor_ref = boost::multi::array_ref<ComplexType, 4, shared_pointer>;

  using mpi3C4Tensor = boost::multi::array<ComplexType, 4, shared_allocator<ComplexType>>;

  using stdCVector     = boost::multi::array<ComplexType, 1>;
  using stdCMatrix     = boost::multi::array<ComplexType, 2>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType, 1>;

  using shm_stack_alloc_type = LocalTGBufferManager::template allocator_t<ComplexType>;
  using StaticSHMVector      = boost::multi::static_array<ComplexType, 1, shm_stack_alloc_type>;
  using StaticSHM4Tensor     = boost::multi::static_array<ComplexType, 4, shm_stack_alloc_type>;

public:
  FullObsHandler(afqmc::TaskGroup_& tg_,
                 AFQMCInfo& info,
                 std::string name_,
                 xmlNodePtr cur,
                 WALKER_TYPES wlk,
                 Wavefunction& wfn)
      : AFQMCInfo(info),
        TG(tg_),
        walker_type(wlk),
        wfn0(wfn),
        writer(false),
        block_size(1),
        nave(1),
        name(name_),
        nspins((walker_type == COLLINEAR) ? 2 : 1),
        G4D_host({0, 0, 0, 0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    using std::fill_n;

    xmlNodePtr curRoot = cur;
    if (cur != NULL)
    {
      ParameterSet m_param;
      m_param.add(nave, "naverages");
      m_param.add(block_size, "block_size");
      m_param.put(cur);
    }

    if (nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    cur = curRoot->children;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      std::transform(cname.begin(), cname.end(), cname.begin(), (int (*)(int))tolower);
      if (cname == "onerdm")
      {
        properties.emplace_back(Observable(full1rdm(TG, info, cur, walker_type, nave, block_size)));
      }
      else if (cname == "gfock" || cname == "genfock" || cname == "ekt")
      {
        properties.emplace_back(Observable(
            generalizedFockMatrix(TG, info, cur, walker_type, wfn0.getHamiltonianOperations(), nave, block_size)));
      }
      else if (cname == "diag2rdm")
      {
        properties.emplace_back(Observable(diagonal2rdm(TG, info, cur, walker_type, nave, block_size)));
      }
      else if (cname == "twordm")
      {
        properties.emplace_back(Observable(full2rdm(TG, info, cur, walker_type, nave, block_size)));
      }
      else if (cname == "n2r" || cname == "ontop2rdm")
      {
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
        std::string str("false");
        ParameterSet m_param;
        m_param.add(str, "use_host_memory");
        m_param.put(cur);
        std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);
        if (str == "false" || str == "no")
        {
          properties.emplace_back(Observable(
              n2r<device_allocator<ComplexType>>(TG, info, cur, walker_type, false, device_allocator<ComplexType>{},
                                                 device_allocator<ComplexType>{}, nave, block_size)));
        }
        else
#endif
        {
          properties.emplace_back(Observable(
              n2r<shared_allocator<ComplexType>>(TG, info, cur, walker_type, true,
                                                 shared_allocator<ComplexType>{TG.TG_local()},
                                                 shared_allocator<ComplexType>{TG.Node()}, nave, block_size)));
        }
      }
      else if (cname == "realspace_correlators")
      {
        properties.emplace_back(Observable(realspace_correlators(TG, info, cur, walker_type, nave, block_size)));
      }
      else if (cname == "correlators")
      {
        properties.emplace_back(Observable(atomcentered_correlators(TG, info, cur, walker_type, nave, block_size)));
      }
      cur = cur->next;
    }

    if (properties.size() == 0)
      APP_ABORT("empty observables list is not allowed.\n");

    Gdims = std::make_tuple(NMO, NMO);
    if (walker_type == NONCOLLINEAR)
      Gdims = std::make_tuple(2 * NMO, 2 * NMO);
    dm_size = nspins * std::get<0>(Gdims) * std::get<1>(Gdims);

    writer = (TG.getGlobalRank() == 0);

    denominator = stdCVector(iextensions<1u>{nave});
    fill_n(denominator.begin(), denominator.num_elements(), ComplexType(0.0, 0.0));
  }

  void print(int iblock, hdf_archive& dump)
  {
    using std::fill_n;

    if (TG.TG_local().root())
    {
      ma::scal(ComplexType(1.0 / block_size), denominator);
      TG.TG_heads().reduce_in_place_n(to_address(denominator.origin()), denominator.num_elements(), std::plus<>(), 0);
    }

    for (auto& v : properties)
      v.print(iblock, dump, denominator);
    fill_n(denominator.origin(), denominator.num_elements(), ComplexType(0.0, 0.0));
  }

  template<class WlkSet, class MatR, class CVec, class MatD>
  void accumulate(int iav, WlkSet& wset, MatR&& Refs, CVec&& wgt, MatD&& DevdetR, bool impsamp)
  {
    if (iav < 0 || iav >= nave)
      APP_ABORT("Runtime Error: iav out of range in full1rdm::accumulate. \n\n\n");

    int nw(wset.size());
    int nrefs(std::get<1>(Refs.sizes()));
    double LogOverlapFactor(wset.getLogOverlapFactor());
    LocalTGBufferManager shm_buffer_manager;
    StaticSHM4Tensor G4D({nw, nspins, std::get<0>(Gdims), std::get<1>(Gdims)},
                         shm_buffer_manager.get_generator().template get_allocator<ComplexType>());
    StaticSHMVector DevOv(iextensions<1u>{2 * nw},
                          shm_buffer_manager.get_generator().template get_allocator<ComplexType>());
    sharedCMatrix_ref G2D(G4D.origin(), {nw, dm_size});

    if (G4D_host.num_elements() != G4D.num_elements())
    {
      G4D_host = mpi3C4Tensor(G4D.extensions(), shared_allocator<ComplexType>{TG.TG_local()});
      TG.TG_local().barrier();
    }

    stdCVector Xw(iextensions<1u>{nw});
    std::fill_n(Xw.origin(), Xw.num_elements(), ComplexType(1.0, 0.0));
    stdCVector Ov(iextensions<1u>{2 * nw});
    stdCMatrix detR(DevdetR);

    using SMType = typename WlkSet::reference::SMType;
    // MAM: The pointer type of GA/GB needs to be device_ptr, it can not be
    //      one of the shared_memory types. The dispatching in DensityMatrices is done
    //      through the pointer type of the result matrix (GA/GB).
    std::vector<devCMatrix_ptr> GA;
    std::vector<devCMatrix_ptr> GB;
    std::vector<SMType> RefsA;
    std::vector<SMType> RefsB;
    std::vector<SMType> SMA;
    std::vector<SMType> SMB;
    GA.reserve(nw);
    SMA.reserve(nw);
    RefsA.reserve(nw);
    if (walker_type == COLLINEAR)
      RefsB.reserve(nw);
    if (walker_type == COLLINEAR)
      SMB.reserve(nw);
    if (walker_type == COLLINEAR)
      GB.reserve(nw);

    if (impsamp)
      denominator[iav] += std::accumulate(wgt.begin(), wgt.end(), ComplexType(0.0));
    else
    {
      APP_ABORT(" Finish implementation of free projection. \n\n\n");
    }

    for (int iref = 0, is = 0; iref < nrefs; iref++, is += nspins)
    {
      // conjugated here!
      ComplexType CIcoeff(std::conj(wfn0.getReferenceWeight(iref)));

      //1. Calculate Green functions
      // Refs({wset.size(),nrefs,ref_size}
      RefsA.clear();
      RefsB.clear();
      SMA.clear();
      SMB.clear();
      GA.clear();
      GB.clear();
      // using SlaterMatrixAux to store References in device memory
      if (walker_type == COLLINEAR)
      {
        for (int iw = 0; iw < nw; iw++)
        {
          SMA.emplace_back(wset[iw].SlaterMatrixN(Alpha));
          SMB.emplace_back(wset[iw].SlaterMatrixN(Beta));
          GA.emplace_back(make_device_ptr(G2D[iw].origin()), iextensions<2u>{NMO, NMO});
          GB.emplace_back(make_device_ptr(G2D[iw].origin()) + NMO * NMO, iextensions<2u>{NMO, NMO});
          RefsA.emplace_back(wset[iw].SlaterMatrixAux(Alpha));
          RefsB.emplace_back(wset[iw].SlaterMatrixAux(Beta));
          copy_n(Refs[iw][iref].origin(), (*RefsA.back()).num_elements(), (*RefsA.back()).origin());
          copy_n(Refs[iw][iref].origin() + (*RefsA.back()).num_elements(), (*RefsB.back()).num_elements(),
                 (*RefsB.back()).origin());
        }
        wfn0.DensityMatrix(RefsA, SMA, GA, DevOv.sliced(0, nw), LogOverlapFactor, false, false);
        wfn0.DensityMatrix(RefsB, SMB, GB, DevOv.sliced(nw, 2 * nw), LogOverlapFactor, false, false);
      }
      else
      {
        for (int iw = 0; iw < nw; iw++)
        {
          SMA.emplace_back(wset[iw].SlaterMatrixN(Alpha));
          GA.emplace_back(make_device_ptr(G2D[iw].origin()), iextensions<2u>{NMO, NMO});
          RefsA.emplace_back(wset[iw].SlaterMatrixAux(Alpha));
          copy_n(Refs[iw][iref].origin(), (*RefsA.back()).num_elements(), (*RefsA.back()).origin());
        }
        wfn0.DensityMatrix(RefsA, SMA, GA, DevOv.sliced(0, nw), LogOverlapFactor, false, false);
      }

      //2. calculate and accumulate appropriate weights
      copy_n(DevOv.origin(), 2 * nw, Ov.origin());
      if (nrefs > 1)
      {
        if (walker_type == CLOSED)
        {
          for (int iw = 0; iw < nw; iw++)
            Xw[iw] = CIcoeff * Ov[iw] * Ov[iw] * std::conj(detR[iw][iref] * detR[iw][iref]);
        }
        else if (walker_type == COLLINEAR)
        {
          for (int iw = 0; iw < nw; iw++)
            Xw[iw] = CIcoeff * Ov[iw] * Ov[iw + nw] * std::conj(detR[iw][2 * iref] * detR[iw][2 * iref + 1]);
        }
        else if (walker_type == NONCOLLINEAR)
        {
          for (int iw = 0; iw < nw; iw++)
            Xw[iw] = CIcoeff * Ov[iw] * std::conj(detR[iw][iref]);
        }
      }
      if (nrefs == 1)
        for (int iw = 0; iw < nw; iw++)
          Xw[iw] = ComplexType(1.0);

      // MAM: Since most of the simpler estimators need G4D in host memory,
      //      I'm providing a copy of the structure there already
      TG.TG_local().barrier();
      int i0, iN;
      std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), int(G4D_host.num_elements()), TG.TG_local().size());
      copy_n(make_device_ptr(G4D.origin()) + i0, iN - i0, to_address(G4D_host.origin()) + i0);
      TG.TG_local().barrier();

      //3. accumulate references
      for (auto& v : properties)
        v.accumulate_reference(iav, iref, G4D, G4D_host, wgt, Xw, Ov, impsamp);
    }
    //4. accumulate block (normalize and accumulate sum over references)
    for (auto& v : properties)
      v.accumulate_block(iav, wgt, impsamp);
  }

private:
  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  Wavefunction& wfn0;

  bool writer;

  int block_size;

  int nave;

  std::string name;

  int nspins;
  int dm_size;
  std::tuple<int, int> Gdims;

  std::vector<Observable> properties;

  // denominator (nave, ...)
  stdCVector denominator;

  // space for G in host space
  mpi3C4Tensor G4D_host;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
