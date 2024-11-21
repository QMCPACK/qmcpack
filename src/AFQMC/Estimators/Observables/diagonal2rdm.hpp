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

#ifndef QMCPLUSPLUS_AFQMC_DIAGONAL2RDM_HPP
#define QMCPLUSPLUS_AFQMC_DIAGONAL2RDM_HPP

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

namespace qmcplusplus
{
namespace afqmc
{
/* 
 * Observable class that calculates the walker averaged diagonal of the 2 RDM 
 */
class diagonal2rdm : public AFQMCInfo
{
  // allocators
  using Allocator = device_allocator<ComplexType>;

  // type defs
  using pointer       = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CVector_ref    = boost::multi::array_ref<ComplexType, 1, pointer>;
  using CMatrix_ref    = boost::multi::array_ref<ComplexType, 2, pointer>;
  using CVector        = boost::multi::array<ComplexType, 1, Allocator>;
  using CMatrix        = boost::multi::array<ComplexType, 2, Allocator>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType, 1>;
  using stdCMatrix_ref = boost::multi::array_ref<ComplexType, 2>;
  using mpi3CVector    = boost::multi::array<ComplexType, 1, shared_allocator<ComplexType>>;
  using mpi3CMatrix    = boost::multi::array<ComplexType, 2, shared_allocator<ComplexType>>;
  using mpi3CTensor    = boost::multi::array<ComplexType, 3, shared_allocator<ComplexType>>;
  using mpi3C4Tensor   = boost::multi::array<ComplexType, 4, shared_allocator<ComplexType>>;

public:
  diagonal2rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk, int nave_ = 1, int bsize = 1)
      : AFQMCInfo(info),
        block_size(bsize),
        nave(nave_),
        counter(0),
        TG(tg_),
        walker_type(wlk),
        writer(false),
        hdf_walker_output(""),
        DMAverage({0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork({0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        denom(iextensions<1u>{0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    app_log() << "  --  Adding Back Propagated Diagonal 2RDM (Diag2RDM) estimator. -- \n";

    if (cur != NULL)
    {
      ParameterSet m_param;
      m_param.add(hdf_walker_output, "walker_output");
      m_param.put(cur);
    }

    if (hdf_walker_output != std::string(""))
    {
      hdf_walker_output = "G" + std::to_string(TG.TG_heads().rank()) + "_" + hdf_walker_output;
      hdf_archive dump;
      if (not dump.create(hdf_walker_output))
      {
        app_log() << "Problems creating walker output hdf5 file: " << hdf_walker_output << std::endl;
        APP_ABORT("Problems creating walker output hdf5 file.\n");
      }
      dump.push("DiagTwoRDM");
      dump.push("Metadata");
      dump.write(NMO, "NMO");
      dump.write(NAEA, "NUP");
      dump.write(NAEB, "NDOWN");
      int wlk_t_copy = walker_type; // the actual data type of enum is implementation-defined. convert to int for file
      dump.write(wlk_t_copy, "WalkerType");
      dump.pop();
      dump.pop();
      dump.close();
    }

    using std::fill_n;
    writer  = (TG.getGlobalRank() == 0);
    dm_size = NMO * (2 * NMO - 1);
    if (walker_type == CLOSED)
      dm_size -= NMO * (NMO - 1) / 2;

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

    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), dm_size, TG.TG_local().size());

    for (int iw = 0; iw < nw; iw++)
    {
      if (iw % TG.TG_local().size() == TG.TG_local().rank())
      {
        denom[iw] += Xw[iw];
        ComplexType* ptr(DMWork[iw].origin());
        if (walker_type == CLOSED)
        {
          auto&& Gu_(G_host[iw][0]);
          for (int i = 0, ij = 0; i < NMO; i++)
          {
            for (int j = i + 1; j < NMO; j++, ptr++)
              *ptr += Xw[iw] * (Gu_[i][i] * Gu_[j][j] - Gu_[i][j] * Gu_[j][i]);
            for (int j = NMO, j0 = 0; j < 2 * NMO; j++, j0++, ptr++)
              *ptr += Xw[iw] * (Gu_[i][i] * Gu_[j0][j0]);
          }
        }
        else if (walker_type == COLLINEAR)
        {
          auto&& Gu_(G_host[iw][0]);
          auto&& Gd_(G_host[iw][1]);
          for (int i = 0, ij = 0; i < 2 * NMO; i++)
          {
            if (i < NMO)
            {
              for (int j = i + 1; j < NMO; j++, ptr++)
                *ptr += Xw[iw] * (Gu_[i][i] * Gu_[j][j] - Gu_[i][j] * Gu_[j][i]);
              for (int j = NMO, j0 = 0; j < 2 * NMO; j++, j0++, ptr++)
                *ptr += Xw[iw] * (Gu_[i][i] * Gd_[j0][j0]);
            }
            else
            {
              int i0 = i - NMO;
              for (int j = i + 1, j0 = i + 1 - NMO; j < 2 * NMO; j++, j0++, ptr++)
                *ptr += Xw[iw] * (Gd_[i0][i0] * Gd_[j0][j0] - Gd_[i0][j0] * Gd_[j0][i0]);
              ;
            }
          }
        }
        else
        {
          auto&& G_(G_host[iw][0]);
          for (int i = 0, ij = 0; i < 2 * NMO; i++)
            for (int j = i + 1; j < 2 * NMO; j++, ptr++)
              *ptr += Xw[iw] * (G_[i][i] * G_[j][j] - G_[i][j] * G_[j][i]);
        }
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

    if (hdf_walker_output != std::string(""))
    {
      const int n_zero = 9;
      for (int iw = 0; iw < nw; iw++)
        ma::scal(ComplexType(1.0, 0.0) / denom[iw], DMWork[iw].sliced(i0, iN));
      TG.TG_local().barrier();

      if (TG.TG_local().root())
      {
        hdf_archive dump;
        if (iav == 0)
          counter++;
        if (not dump.open(hdf_walker_output))
        {
          app_log() << "Problems opening walker output hdf5 file: " << hdf_walker_output << std::endl;
          APP_ABORT("Problems opening walker output hdf5 file.\n");
        }
        dump.push("DiagTwoRDM");
        dump.push(std::string("Group") + std::to_string(TG.TG_heads().rank()));
        dump.push(std::string("Average_") + std::to_string(iav));
        std::string padded_num = std::string(n_zero - std::to_string(counter).length(), '0') + std::to_string(counter);
        dump.write(wgt, "weights_" + padded_num);
        stdCMatrix_ref DM(to_address(DMWork.origin()), {nw, dm_size});
        dump.write(DM, "diag_two_rdm_" + padded_num);
        dump.pop();
        dump.pop();
        dump.pop();
        dump.close();

        // adjust denom
        for (int iw = 0; iw < nw; iw++)
          denom[iw] = wgt[iw];
      }
      TG.TG_local().barrier();
    }
    else
    {
      if (TG.TG_local().root())
        for (int iw = 0; iw < nw; iw++)
          denom[iw] = wgt[iw] / denom[iw];
    }
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
        dump.push(std::string("DiagTwoRDM"));
        for (int i = 0; i < nave; ++i)
        {
          dump.push(std::string("Average_") + std::to_string(i));
          std::string padded_iblock =
              std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
          stdCVector_ref DMAverage_(to_address(DMAverage[i].origin()), {dm_size});
          dump.write(DMAverage_, "diag_two_rdm_" + padded_iblock);
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
  int block_size;

  int nave;

  int counter;

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  int dm_size;

  bool writer;

  std::string hdf_walker_output;

  // DMAverage (nave, spin*spin*x*NMO*(x*NMO-1)/2 ), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMAverage;
  // DMWork (nwalk, spin*spin*x*NMO*(x*NMO-1)/2 ), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMWork;

  mpi3CVector denom;
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
