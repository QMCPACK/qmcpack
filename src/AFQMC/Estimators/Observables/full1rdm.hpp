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

#ifndef QMCPLUSPLUS_AFQMC_FULL1RDM_HPP
#define QMCPLUSPLUS_AFQMC_FULL1RDM_HPP

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
 * Observable class that calculates the walker averaged 1 RDM.
 * The resulting RDM will be [spin][x*NMO][x*NMO],
 * where x:2 for NONCOLLINEAR and 1 for everything else.
 */
class full1rdm : public AFQMCInfo
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
  using StaticMatrix     = boost::multi::static_array<ComplexType, 2, stack_alloc_type>;

public:
  full1rdm(afqmc::TaskGroup_& tg_, AFQMCInfo& info, xmlNodePtr cur, WALKER_TYPES wlk, int nave_ = 1, int bsize = 1)
      : AFQMCInfo(info),
        TG(tg_),
        walker_type(wlk),
        writer(false),
        block_size(bsize),
        nave(nave_),
        counter(0),
        hdf_walker_output(""),
        nskip_walker_output(0),
        ncnt_walker_output(0),
        apply_rotation(false),
        XRot({0, 0}, make_node_allocator<ComplexType>(TG)),
        print_from_list(false),
        index_list({0, 0}, shared_allocator<int>{TG.Node()}),
        denom(iextensions<1u>{0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMAverage({0, 0}, shared_allocator<ComplexType>{TG.TG_local()}),
        DMWork({0, 0}, shared_allocator<ComplexType>{TG.TG_local()})
  {
    using std::copy_n;
    using std::fill_n;

    app_log() << "  --  Adding Back Propagated Full 1RDM (OneRDM) estimator. -- \n";

    std::string rot_file("");
    std::string path("/");
    std::string str("false");

    if (cur != NULL)
    {
      ParameterSet m_param;
      m_param.add(hdf_walker_output, "walker_output");
      m_param.add(nskip_walker_output, "nskip_output");
      m_param.add(rot_file, "rotation");
      m_param.add(path, "path");
      m_param.add(str, "with_index_list");
      m_param.put(cur);
    }

    if (rot_file != "")
    {
      if (not file_exists(rot_file))
      {
        app_error() << " Error: File with rotation matrix does not exist: " << rot_file << std::endl;
        APP_ABORT("");
      }
      apply_rotation  = true;
      print_from_list = (str == "true" || str == "yes");
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
        stdIMatrix I;
        if (print_from_list)
        {
          if (!dump.readEntry(I, "Indices"))
            APP_ABORT("Error reading Indices.\n");
          if (std::get<1>(I.sizes()) != 2)
            APP_ABORT("Error Wrong dimensions in Indices.\n");
          dim[1] = std::get<0>(I.sizes());
        }
        TG.Node().broadcast_n(dim, 2, 0);
        XRot = sharedCMatrix({dim[0], NMO}, make_node_allocator<ComplexType>(TG));
        copy_n(R.origin(), R.num_elements(), make_device_ptr(XRot.origin()));
        if (TG.Node().root())
          TG.Cores().broadcast_n(to_address(XRot.origin()), XRot.num_elements(), 0);
        if (print_from_list)
        {
          index_list = mpi3IMatrix({dim[1], 2}, shared_allocator<int>{TG.Node()});
          copy_n(I.origin(), I.num_elements(), make_device_ptr(index_list.origin()));
          if (TG.Node().root())
            TG.Cores().broadcast_n(to_address(index_list.origin()), index_list.num_elements(), 0);
        }

        dump.pop();
        dump.close();
      }
      else
      {
        TG.Node().broadcast_n(dim, 2, 0);
        XRot = sharedCMatrix({dim[0], NMO}, make_node_allocator<ComplexType>(TG));
        if (TG.Node().root())
          TG.Cores().broadcast_n(to_address(XRot.origin()), XRot.num_elements(), 0);
        if (print_from_list)
        {
          index_list = mpi3IMatrix({dim[1], 2}, shared_allocator<int>{TG.Node()});
          if (TG.Node().root())
            TG.Cores().broadcast_n(to_address(index_list.origin()), index_list.num_elements(), 0);
        }
      }
      TG.Node().barrier();

      if (print_from_list)
        dm_size = index_list.size();
      else
        dm_size = XRot.size() * XRot.size();
    }
    else
    {
      // can also add print_from_list option without rotation later on
      dm_size = NMO * NMO;
    }

    if (walker_type == COLLINEAR)
      dm_size *= 2;
    else if (walker_type == NONCOLLINEAR)
      dm_size *= 4;

    if (hdf_walker_output != std::string(""))
    {
      hdf_walker_output = "G" + std::to_string(TG.TG_heads().rank()) + "_" + hdf_walker_output;
      hdf_archive dump;
      if (not dump.create(hdf_walker_output))
      {
        app_log() << "Problems creating walker output hdf5 file: " << hdf_walker_output << std::endl;
        APP_ABORT("Problems creating walker output hdf5 file.\n");
      }
      dump.push("FullOneRDM");
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
      acc_no_rotation(G_host, Xw);
  }


  template<class HostCVec>
  void accumulate_block(int iav, HostCVec&& wgt, bool impsamp)
  {
    int nw(denom.size());
    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), dm_size, TG.TG_local().size());

    if (iav == 0)
      ncnt_walker_output++;
    if (hdf_walker_output != std::string("") && ncnt_walker_output % (nskip_walker_output + 1) == 0)
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
        dump.push("FullOneRDM");
        dump.push(std::string("Group") + std::to_string(TG.TG_heads().rank()));
        dump.push(std::string("Average_") + std::to_string(iav));
        std::string padded_num = std::string(n_zero - std::to_string(counter).length(), '0') + std::to_string(counter);
        dump.write(wgt, "weights_" + padded_num);
        stdCMatrix_ref DM(to_address(DMWork.origin()), {nw, dm_size});
        dump.write(DM, "one_rdm_" + padded_num);
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
    //app_log()<<" DMAv: " <<iav <<" " <<ma::sum(DMAverage[iav]) <<"\n";
    //TG.Global().barrier();
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
        dump.push(std::string("FullOneRDM"));
        for (int i = 0; i < nave; ++i)
        {
          dump.push(std::string("Average_") + std::to_string(i));
          std::string padded_iblock =
              std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
          stdCVector_ref DMAverage_(to_address(DMAverage[i].origin()), {dm_size});
          dump.write(DMAverage_, "one_rdm_" + padded_iblock);
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

  int dm_size;

  std::string hdf_walker_output;

  int nskip_walker_output;

  int ncnt_walker_output;

  bool apply_rotation;

  sharedCMatrix XRot;
  stdCVector Grot;

  bool print_from_list;

  mpi3IMatrix index_list;

  mpi3CVector denom;

  // DMAverage (nave, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMAverage;
  // DMWork (nwalk, spin*x*NMO*x*NMO), x=(1:CLOSED/COLLINEAR, 2:NONCOLLINEAR)
  mpi3CMatrix DMWork;

  // expects host accesible memory
  template<class MatG, class CVec>
  void acc_no_rotation(MatG&& G, CVec&& Xw)
  {
    int nw(G.size());
    assert(G[0].num_elements() == dm_size);

    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), dm_size, TG.TG_local().size());

    stdCMatrix_ref G2D(to_address(G.origin()), {nw, dm_size});

    for (int iw = 0; iw < nw; iw++)
    {
      if (TG.TG_local().root())
        denom[iw] += Xw[iw];
      ma::axpy(Xw[iw], G2D[iw].sliced(i0, iN), DMWork[iw].sliced(i0, iN));
    }
    TG.TG_local().barrier();
    //app_log()<<" DMWORK: " <<ma::sum(DMWork) <<"\n";
    //TG.Global().barrier();
  }

  // G should be device accesible memory
  // Xw is in host
  template<class MatG, class CVec>
  void acc_with_rotation(MatG&& G, CVec&& Xw)
  {
    int nw(G.size());
    assert(std::get<2>(G.sizes()) == std::get<3>(G.sizes()));
    assert(std::get<2>(G.sizes()) == std::get<1>(XRot.sizes()));

    if (walker_type == NONCOLLINEAR)
      APP_ABORT("Error: Not yet implemented: acc_with_rotation && noncollinear.\n");

    int i0, iN;
    std::tie(i0, iN) = FairDivideBoundary(TG.TG_local().rank(), int(XRot.size()), TG.TG_local().size());

    // can batch in the future if too slow
    // Grot = Xc * G * H(Xc)
    int nX   = XRot.size();
    int npts = (iN - i0) * nX;
    DeviceBufferManager buffer_manager;
    StaticMatrix T1({(iN - i0), NMO}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    StaticMatrix T2({(iN - i0), nX}, buffer_manager.get_generator().template get_allocator<ComplexType>());
    if (Grot.size() != npts)
      Grot = stdCVector(iextensions<1u>(npts));

    // round-robin for now
    int cnt = 0;
    for (int iw = 0; iw < nw; iw++)
    {
      if (i0 == iN || i0 == XRot.size())
        break;
      if (TG.TG_local().root())
        denom[iw] += Xw[iw];
      ma::product(XRot.sliced(i0, iN), G[iw][0], T1);
      ma::product(T1, ma::H(XRot), T2);
      copy_n(T2.origin(), T2.num_elements(), Grot.origin());
      if (print_from_list)
      {
        for (int i = 0; i < index_list.size(); i++)
        {
          if (index_list[i][0] >= i0 && index_list[i][0] < iN)
          {
            int ij = (index_list[i][0] - i0) * nX + index_list[i][1];
            DMWork[iw][i] += Xw[iw] * Grot[ij];
          }
        }
      }
      else
        ma::axpy(Xw[iw], Grot, DMWork[iw].sliced(i0 * nX, i0 * nX + npts));
      if (walker_type == COLLINEAR)
      {
        ma::product(XRot.sliced(i0, iN), G[iw][1], T1);
        ma::product(T1, ma::H(XRot), T2);
        copy_n(T2.origin(), T2.num_elements(), Grot.origin());
        if (print_from_list)
        {
          for (int i = 0, ie = index_list.size(); i < ie; i++)
          {
            if (index_list[i][0] >= i0 && index_list[i][0] < iN)
            {
              int ij = (index_list[i][0] - i0) * nX + index_list[i][1];
              DMWork[iw][i + ie] += Xw[iw] * Grot[ij];
            }
          }
        }
        else
          ma::axpy(Xw[iw], Grot, DMWork[iw].sliced((nX + i0) * nX, (nX + i0) * nX + npts));
      }
    }
    TG.TG_local().barrier();
  }
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
