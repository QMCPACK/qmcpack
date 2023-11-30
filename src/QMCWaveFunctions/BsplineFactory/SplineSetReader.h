//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * The most general reader class for the following classes using the full single grid for the supercell
 * - SplineR2R
 * - SplineC2C
 * - SplineC2R
 * Each band is initialized with UBspline_3d_d and both real and imaginary parts are passed to the objects
 * which will convert the data type to their internal precision. 
 */
#ifndef QMCPLUSPLUS_SPLINESET_READER_H
#define QMCPLUSPLUS_SPLINESET_READER_H
#include "mpi/collectives.h"
#include "mpi/point2point.h"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
class OneSplineOrbData
{
  Array<std::complex<double>, 3> FFTbox;
  Array<double, 3> splineData_r, splineData_i;
  double rotate_phase_r, rotate_phase_i;
  UBspline_3d_d* spline_r = nullptr;
  UBspline_3d_d* spline_i = nullptr;
  fftw_plan FFTplan       = nullptr;

  const TinyVector<int, 3>& mesh_size_;
  const bool isComplex_;

  void create(const TinyVector<int, 3>& halfG)
  {
    const int nx = mesh_size_[0];
    const int ny = mesh_size_[1];
    const int nz = mesh_size_[2];
    //perform FFT using FFTW
    FFTbox.resize(nx, ny, nz);
    FFTplan = fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(FFTbox.data()),
                               reinterpret_cast<fftw_complex*>(FFTbox.data()), +1, FFTW_ESTIMATE);
    splineData_r.resize(nx, ny, nz);
    if (isComplex_)
      splineData_i.resize(nx, ny, nz);

    TinyVector<double, 3> start(0.0);
    TinyVector<double, 3> end(1.0);
    spline_r = einspline::create(spline_r, start, end, mesh_size_, halfG);
    if (isComplex_)
      spline_i = einspline::create(spline_i, start, end, mesh_size_, halfG);
  }

public:
  OneSplineOrbData(const TinyVector<int, 3>& mesh_size, const TinyVector<int, 3>& halfG, const bool isComplex)
      : mesh_size_(mesh_size), isComplex_(isComplex)
  {
    create(halfG);
  }

  ~OneSplineOrbData() { clear(); }

  auto getRotatePhase() const { return std::complex<double>(rotate_phase_r, rotate_phase_i); }
  auto& get_spline_r() { return *spline_r; }
  auto& get_spline_i() { return *spline_i; }

  void clear()
  {
    einspline::destroy(spline_r);
    einspline::destroy(spline_i);
    if (FFTplan != nullptr)
      fftw_destroy_plan(FFTplan);
    FFTplan = nullptr;
  }

  /** fft and spline cG
   * @param cG psi_g to be processed
   * @param ti twist index
   * @param iorb orbital index
   *
   * Perform FFT and spline to spline_r and spline_i
   */
  void fft_spline(const Vector<std::complex<double>>& cG,
                  const std::vector<TinyVector<int, 3>>& gvecs,
                  const TinyVector<double, 3>& primcell_kpoint,
                  const bool rotate)
  {
    unpack4fftw(cG, gvecs, mesh_size_, FFTbox);
    fftw_execute(FFTplan);
    if (isComplex_)
    {
      if (rotate)
        fix_phase_rotate_c2c(FFTbox, splineData_r, splineData_i, primcell_kpoint, rotate_phase_r, rotate_phase_i);
      else
      {
        split_real_components_c2c(FFTbox, splineData_r, splineData_i);
        rotate_phase_r = 1.0;
        rotate_phase_i = 0.0;
      }
      einspline::set(spline_r, splineData_r.data());
      einspline::set(spline_i, splineData_i.data());
    }
    else
    {
      fix_phase_rotate_c2r(FFTbox, splineData_r, primcell_kpoint, rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r, splineData_r.data());
    }
  }
};

/** General SplineSetReader to handle any unitcell
 */
template<typename SA>
class SplineSetReader : public BsplineReader
{
public:
  using DataType   = typename SA::DataType;
  using SplineType = typename SA::SplineType;

  SplineSetReader(EinsplineSetBuilder* e) : BsplineReader(e) {}

  std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                            int spin,
                                            const BandInfoGroup& bandgroup) override
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
        int sizeD = sizeof(DataType);
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

  /** create data space in the spline object and try open spline dump files.
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   * @return true if dumpfile pass class name and data type size check
   */
  bool createSplineDataSpaceLookforDumpFile(const BandInfoGroup& bandgroup, SA& bspline) const
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
        foundspline = (sizeD == sizeof(DataType));
      }
      h5f.close();
    }
    myComm->bcast(foundspline);
    return foundspline;
  }

  /** read planewave coefficients from h5 file
   * @param s data set full path in h5
   * @param h5f hdf5 file handle
   * @param cG vector to store coefficients
   */
  void readOneOrbitalCoefs(const std::string& s, hdf_archive& h5f, Vector<std::complex<double>>& cG) const
  {
    if (!h5f.readEntry(cG, s))
    {
      std::ostringstream msg;
      msg << "SplineSetReader Failed to read band(s) from h5 file. " << "Attempted dataset " << s << " with "
          << cG.size() << " complex numbers." << std::endl;
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

  /** transforming planewave orbitals to 3D B-spline orbitals in real space.
   * @param spin orbital dataset spin index
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   */
  void initialize_spline_pio_gather(const int spin, const BandInfoGroup& bandgroup, SA& bspline) const
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
};
} // namespace qmcplusplus
#endif
