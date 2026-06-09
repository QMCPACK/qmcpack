//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SplineSetReader.h"
#include "OneSplineOrbData.hpp"
#include "Utilities/FairDivide.h"
#include <Timer.h>
#if defined(QMC_COMPLEX)
#include "SplineC2C.h"
#include "SplineC2COMPTarget.h"
#else
#include "SplineR2R.h"
#include "SplineC2R.h"
#include "SplineC2ROMPTarget.h"
#endif
#include "Message/CommOperators.h"
#include "spline2/SplineUtils.h"
#include "spline2/MultiBspline.hpp"
#include "spline2/MultiBsplineOffload.hpp"
#if defined(HAVE_MPI)
#include "spline2/MultiBsplineMPIShared.hpp"
#endif


namespace qmcplusplus
{
template<typename ST>
SplineSetReader<ST>::SplineSetReader(EinsplineSetBuilder* e, bool use_duplex_splines)
    : BsplineReader(e, use_duplex_splines)
{}

template<typename ST>
std::unique_ptr<SPOSet> SplineSetReader<ST>::create_spline_set(const std::string& my_name,
                                                               int spin,
                                                               const std::pair<int, int>& distributed_and_shared_ranks,
                                                               const BandInfoGroup& bandgroup)
{
  const int N = bandgroup.getNumDistinctOrbitals();

  auto [distributed_ranks, shared_ranks] = distributed_and_shared_ranks;

  if (use_offload)
  {
    if (distributed_ranks > 1 || shared_ranks > 1)
      app_warning() << "Offload implemenation doesn't support distributing or sharing the memory of spline "
                       "coefficients. Overriding distributed_ranks and shared_ranks to 1."
                    << std::endl;
    distributed_ranks = 1;
    shared_ranks      = 1;
  }

  auto dist_comm_ptr = std::make_unique<Communicate>(*myComm, myComm->size() / (distributed_ranks * shared_ranks));

  app_log() << "  Using " << (use_duplex_splines_ ? "complex" : "real") << " einspline table." << std::endl;
  if (distributed_ranks > 1)
    app_log() << "  Distributed across " << distributed_ranks << " MPI ranks." << std::endl;
  if (shared_ranks > 1)
    app_log() << "  Shared across " << shared_ranks << " MPI ranks." << std::endl;

  const TinyVector<int, 3> half_g = use_duplex_splines_
      ? TinyVector<int, 3>(0, 0, 0)
      : computeHalfG(mybuilder->TargetPtcl.getLattice().BoxBConds, mybuilder->primcell_kpoints,
                     bandgroup.myBands[0].TwistIndex);

  Ugrid xyz_grid[3];
  typename bspline_traits<ST, 3>::BCType xyz_bc[3];
  set_grid(mybuilder->MeshSize, half_g, xyz_grid, xyz_bc);

  const size_t num_splines = getAlignedSize<ST>(use_duplex_splines_ ? N * 2 : N);
  std::unique_ptr<MultiBsplineBase<ST>> multi_splines_ptr;
  if (use_offload)
    multi_splines_ptr = std::make_unique<MultiBsplineOffload<ST>>(xyz_grid, xyz_bc, num_splines);
#if defined(HAVE_MPI)
  else if (distributed_ranks * shared_ranks > 1)
    multi_splines_ptr = std::make_unique<MultiBsplineMPIShared<ST>>(xyz_grid, xyz_bc, num_splines,
                                                                    std::move(dist_comm_ptr), distributed_ranks);
#endif
  else
    multi_splines_ptr = std::make_unique<MultiBspline<ST>>(xyz_grid, xyz_bc, num_splines);

  auto& multi_splines(*multi_splines_ptr);
  app_log() << "MEMORY " << multi_splines.sizeInByte() / (1 << 20) << " MB allocated "
            << "for the coefficients in 3D spline orbital representation" << std::endl;

  std::unique_ptr<BsplineSet> bspline;
#if defined(QMC_COMPLEX)
  if (use_offload)
    bspline = std::make_unique<SplineC2COMPTarget<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                       std::move(multi_splines_ptr), use_offload);
  else
    bspline = std::make_unique<SplineC2C<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                              std::move(multi_splines_ptr), use_offload);
#else
  if (use_duplex_splines_)
  {
    if (use_offload)
      bspline = std::make_unique<SplineC2ROMPTarget<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                         std::move(multi_splines_ptr), use_offload);
    else
      bspline = std::make_unique<SplineC2R<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                std::move(multi_splines_ptr), use_offload);
  }
  else
    bspline = std::make_unique<SplineR2R<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                              std::move(multi_splines_ptr), use_offload);
#endif

  app_log() << "  ClassName = " << bspline->getClassName() << std::endl;

  bspline->resizeStorage(N);

  bspline->HalfG = half_g;

  if (use_duplex_splines_)
  {
    //baseclass handles twists
    check_twists(*bspline, bandgroup);
  }

  bool foundspline = lookforSplineDataDumpFile(bandgroup, bspline->getKeyword(), sizeof(ST));
  if (foundspline && myComm->rank() == 0)
  {
    Timer now;
    hdf_archive h5f(myComm);
    const auto splinefile = getSplineDumpFileName(bandgroup);
    h5f.open(splinefile, H5F_ACC_RDONLY);
    foundspline = SplineUtils<ST>::read(multi_splines, h5f);
    if (foundspline)
      app_log() << "  Successfully restored 3D B-spline coefficients from " << splinefile << ". The reading time is "
                << now.elapsed() << " sec." << std::endl;
  }

  /* create a sub communicator. spline table is shared across MPI ranks with identical subcomm rank id.
   *  myComm->size() == 12 ; shared_ranks = 2; distributed_ranks = 3;
   *          | group 0 | group 1 |
   *          |  0 1 2  |  3 4 5  |
   *          |  6 7 8  | 9 10 11 |
   */
  Communicate shared_comm(*myComm, shared_ranks, distributed_ranks);
  Communicate dist_comm(shared_comm, shared_comm.size() / distributed_ranks);

  if (!foundspline)
  {
    if (shared_comm.getGroupID() == 0)
    {
      multi_splines.flush_zero(dist_comm.rank());
      Timer now;
      initialize_spline_pio_gather(spin, bandgroup, half_g, bspline->BandIndexMap, multi_splines, dist_comm);
      app_log() << "  SplineSetReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;
    }

    if (saveSplineCoefs && myComm->rank() == 0)
    {
      Timer now;
      const std::string splinefile(getSplineDumpFileName(bandgroup));
      hdf_archive h5f;
      h5f.create(splinefile);
      std::string classname = bspline->getClassName();
      h5f.write(classname, "class_name");
      int sizeD = sizeof(ST);
      h5f.write(sizeD, "sizeof");
      SplineUtils<ST>::write(multi_splines, h5f);
      h5f.close();
      app_log() << "  Stored spline coefficients in " << splinefile << " for potential reuse. The writing time is "
                << now.elapsed() << " sec." << std::endl;
    }
  }

  {
    myComm->barrier();
    Timer now;
    if (shared_comm.getGroupID() == 0)
      SplineUtils<ST>::bcast(multi_splines, dist_comm.rank(), dist_comm.getInterGroupComm());
    myComm->barrier();
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
  }

  return bspline;
}

template<typename ST>
void SplineSetReader<ST>::initialize_spline_pio_gather(const int spin,
                                                       const BandInfoGroup& bandgroup,
                                                       const TinyVector<int, 3>& half_g,
                                                       const aligned_vector<int>& band_index_map,
                                                       MultiBsplineBase<ST>& multi_splines,
                                                       Communicate& dist_comm) const
{
  const auto& block_offsets = multi_splines.getBlockOffsets();
  const size_t iblock       = dist_comm.rank();
  auto& comm                = dist_comm.getInterGroupComm();
  //distribute bands over processor groups
  const int orb_block_start = use_duplex_splines_ ? block_offsets[iblock] / 2 : block_offsets[iblock];
  const int orb_block_end   = use_duplex_splines_ ? block_offsets[iblock + 1] / 2 : block_offsets[iblock + 1];
  const int Nbands          = std::min(orb_block_end, bandgroup.getNumDistinctOrbitals()) - orb_block_start;
  const int Nprocs          = comm.size();
  const int Nbandgroups     = Nbands > 0 ? std::min(Nbands, Nprocs) : 1;

  Communicate band_group_comm(comm, Nbandgroups);
  std::vector<int> band_groups(Nbandgroups + 1, 0);
  FairDivideLow(Nbands, Nbandgroups, band_groups);
  const int iorb_first = orb_block_start + band_groups[band_group_comm.getGroupID()];
  const int iorb_last  = orb_block_start + band_groups[band_group_comm.getGroupID() + 1];

  app_log() << "Start transforming plane waves to 3D B-Splines." << std::endl;
  OneSplineOrbData oneband(mybuilder->MeshSize, half_g, use_duplex_splines_);
  hdf_archive h5f(&band_group_comm, false);
  Vector<std::complex<double>> cG(mybuilder->Gvecs[0].size());
  const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
  if (band_group_comm.isGroupLeader())
  {
    auto& group_leader_comm = band_group_comm.getInterGroupComm();
    h5f.open(mybuilder->H5FileName, H5F_ACC_RDONLY);
    for (int iorb = iorb_first; iorb < iorb_last; iorb++)
    {
      const auto& cur_band = cur_bands[band_index_map[iorb]];
      const int ti         = cur_band.TwistIndex;
      readOneOrbitalCoefs(psi_g_path(ti, spin, cur_band.BandIndex), h5f, cG);
      oneband.fft_spline(cG, mybuilder->Gvecs[0], mybuilder->primcell_kpoints[ti], rotate);
      if (use_duplex_splines_)
      {
        multi_splines.set_spline(oneband.get_spline_r(), iorb * 2);
        multi_splines.set_spline(oneband.get_spline_i(), iorb * 2 + 1);
      }
      else
        multi_splines.set_spline(oneband.get_spline_r(), iorb);
    }

    {
      group_leader_comm.barrier();
      Timer now;
      if (use_duplex_splines_)
      {
        std::vector<int> offset(band_groups.size());
        for (int i = 0; i < offset.size(); i++)
          offset[i] = band_groups[i] * 2;
        SplineUtils<ST>::gatherv(multi_splines, iblock, offset, group_leader_comm);
      }
      else
        SplineUtils<ST>::gatherv(multi_splines, iblock, band_groups, group_leader_comm);
      app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
    }
  }
}

template class SplineSetReader<float>;
#if !defined(QMC_MIXED_PRECISION)
template class SplineSetReader<double>;
#endif
} // namespace qmcplusplus
