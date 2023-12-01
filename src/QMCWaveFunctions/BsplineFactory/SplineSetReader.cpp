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
#include "einspline_helper.hpp"
#include <Timer.h>
#if defined(QMC_COMPLEX)
#include "SplineC2C.h"
#include "SplineC2COMPTarget.h"
#else
#include "SplineR2R.h"
#include "SplineC2R.h"
#include "SplineC2ROMPTarget.h"
#endif

namespace qmcplusplus
{
template<typename SA>
SplineSetReader<SA>::SplineSetReader(EinsplineSetBuilder* e) : BsplineReader(e)
{}

template<typename SA>
std::unique_ptr<SPOSet> SplineSetReader<SA>::create_spline_set(const std::string& my_name,
                                                               int spin,
                                                               const BandInfoGroup& bandgroup)
{
  auto bspline = std::make_unique<SA>(my_name);
  app_log() << "  ClassName = " << bspline->getClassName() << std::endl;
  bool foundspline = createSplineDataSpaceLookforDumpFile(bandgroup, *bspline);
  if (foundspline)
  {
    Timer now;
    hdf_archive h5f(myComm);
    const auto splinefile = getSplineDumpFileName(bandgroup);
    h5f.open(splinefile, H5F_ACC_RDONLY);
    foundspline = bspline->read_splines(h5f);
    if (foundspline)
      app_log() << "  Successfully restored 3D B-spline coefficients from " << splinefile << ". The reading time is "
                << now.elapsed() << " sec." << std::endl;
  }

  if (!foundspline)
  {
    bspline->flush_zero();

    Timer now;
    initialize_spline_pio_gather(spin, bandgroup, *bspline);
    app_log() << "  SplineSetReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;

    if (saveSplineCoefs && myComm->rank() == 0)
    {
      Timer now;
      const std::string splinefile(getSplineDumpFileName(bandgroup));
      hdf_archive h5f;
      h5f.create(splinefile);
      std::string classname = bspline->getClassName();
      h5f.write(classname, "class_name");
      int sizeD = sizeof(typename SA::DataType);
      h5f.write(sizeD, "sizeof");
      bspline->write_splines(h5f);
      h5f.close();
      app_log() << "  Stored spline coefficients in " << splinefile << " for potential reuse. The writing time is "
                << now.elapsed() << " sec." << std::endl;
    }
  }

  {
    Timer now;
    bspline->bcast_tables(myComm);
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
  }

  return bspline;
}

template<typename SA>
bool SplineSetReader<SA>::createSplineDataSpaceLookforDumpFile(const BandInfoGroup& bandgroup, SA& bspline) const
{
  if (bspline.isComplex())
    app_log() << "  Using complex einspline table" << std::endl;
  else
    app_log() << "  Using real einspline table" << std::endl;

  //baseclass handles twists
  check_twists(bspline, bandgroup);

  Ugrid xyz_grid[3];

  typename SA::BCType xyz_bc[3];
  bool havePsig = set_grid(bspline.HalfG, xyz_grid, xyz_bc);
  if (!havePsig)
    myComm->barrier_and_abort("SplineSetReader needs psi_g. Set precision=\"double\".");
  bspline.create_spline(xyz_grid, xyz_bc);

  int foundspline = 0;
  Timer now;
  if (myComm->rank() == 0)
  {
    now.restart();
    hdf_archive h5f(myComm);
    foundspline = h5f.open(getSplineDumpFileName(bandgroup), H5F_ACC_RDONLY);
    if (foundspline)
    {
      std::string aname("none");
      foundspline = h5f.readEntry(aname, "class_name");
      foundspline = (aname.find(bspline.getKeyword()) != std::string::npos);
    }
    if (foundspline)
    {
      int sizeD   = 0;
      foundspline = h5f.readEntry(sizeD, "sizeof");
      foundspline = (sizeD == sizeof(typename SA::DataType));
    }
    h5f.close();
  }
  myComm->bcast(foundspline);
  return foundspline;
}

template<typename SA>
void SplineSetReader<SA>::readOneOrbitalCoefs(const std::string& s,
                                              hdf_archive& h5f,
                                              Vector<std::complex<double>>& cG) const
{
  if (!h5f.readEntry(cG, s))
  {
    std::ostringstream msg;
    msg << "SplineSetReader Failed to read band(s) from h5 file. "
        << "Attempted dataset " << s << " with " << cG.size() << " complex numbers." << std::endl;
    throw std::runtime_error(msg.str());
  }
  double total_norm = compute_norm(cG);
  if ((checkNorm) && (std::abs(total_norm - 1.0) > PW_COEFF_NORM_TOLERANCE))
  {
    std::ostringstream msg;
    msg << "SplineSetReader The orbital dataset " << s << " has a wrong norm " << total_norm
        << ", computed from plane wave coefficients!" << std::endl
        << "This may indicate a problem with the HDF5 library versions used "
        << "during wavefunction conversion or read." << std::endl;
    throw std::runtime_error(msg.str());
  }
}

template<typename SA>
void SplineSetReader<SA>::initialize_spline_pio_gather(const int spin,
                                                       const BandInfoGroup& bandgroup,
                                                       SA& bspline) const
{
  //distribute bands over processor groups
  int Nbands            = bandgroup.getNumDistinctOrbitals();
  const int Nprocs      = myComm->size();
  const int Nbandgroups = std::min(Nbands, Nprocs);
  Communicate band_group_comm(*myComm, Nbandgroups);
  std::vector<int> band_groups(Nbandgroups + 1, 0);
  FairDivideLow(Nbands, Nbandgroups, band_groups);
  int iorb_first = band_groups[band_group_comm.getGroupID()];
  int iorb_last  = band_groups[band_group_comm.getGroupID() + 1];

  app_log() << "Start transforming plane waves to 3D B-Splines." << std::endl;
  OneSplineOrbData oneband(mybuilder->MeshSize, bspline.HalfG, bspline.isComplex());
  hdf_archive h5f(&band_group_comm, false);
  Vector<std::complex<double>> cG(mybuilder->Gvecs[0].size());
  const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
  if (band_group_comm.isGroupLeader())
  {
    h5f.open(mybuilder->H5FileName, H5F_ACC_RDONLY);
    for (int iorb = iorb_first; iorb < iorb_last; iorb++)
    {
      const auto& cur_band = cur_bands[bspline.BandIndexMap[iorb]];
      const int ti         = cur_band.TwistIndex;
      readOneOrbitalCoefs(psi_g_path(ti, spin, cur_band.BandIndex), h5f, cG);
      oneband.fft_spline(cG, mybuilder->Gvecs[0], mybuilder->primcell_kpoints[ti], rotate);
      bspline.set_spline(&oneband.get_spline_r(), &oneband.get_spline_i(), cur_band.TwistIndex, iorb, 0);
    }

    {
      band_group_comm.getGroupLeaderComm()->barrier();
      Timer now;
      bspline.gather_tables(band_group_comm.getGroupLeaderComm());
      app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
    }
  }
}

#if defined(QMC_COMPLEX)
template class SplineSetReader<SplineC2C<float>>;
template class SplineSetReader<SplineC2COMPTarget<float>>;
#if !defined(QMC_MIXED_PRECISION)
template class SplineSetReader<SplineC2C<double>>;
template class SplineSetReader<SplineC2COMPTarget<double>>;
#endif
#else
template class SplineSetReader<SplineC2R<float>>;
template class SplineSetReader<SplineR2R<float>>;
template class SplineSetReader<SplineC2ROMPTarget<float>>;
#if !defined(QMC_MIXED_PRECISION)
template class SplineSetReader<SplineC2R<double>>;
template class SplineSetReader<SplineR2R<double>>;
template class SplineSetReader<SplineC2ROMPTarget<double>>;
#endif
#endif
} // namespace qmcplusplus
