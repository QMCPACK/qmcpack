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
/** General SplineSetReader to handle any unitcell
 */
template<typename SA>
struct SplineSetReader : public BsplineReaderBase
{
  typedef SA splineset_t;
  typedef typename splineset_t::DataType DataType;
  typedef typename splineset_t::SplineType SplineType;

  Array<std::complex<double>, 3> FFTbox;
  Array<double, 3> splineData_r, splineData_i;
  double rotate_phase_r, rotate_phase_i;
  UBspline_3d_d* spline_r;
  UBspline_3d_d* spline_i;
  splineset_t* bspline;
  fftw_plan FFTplan;

  SplineSetReader(EinsplineSetBuilder* e)
      : BsplineReaderBase(e), spline_r(nullptr), spline_i(nullptr), bspline(nullptr), FFTplan(nullptr)
  {}

  ~SplineSetReader() override { clear(); }

  void clear()
  {
    einspline::destroy(spline_r);
    einspline::destroy(spline_i);
    if (FFTplan != nullptr)
      fftw_destroy_plan(FFTplan);
    FFTplan = nullptr;
  }

  // set info for Hybrid
  virtual void initialize_hybridrep_atomic_centers() {}
  // transform cG to radial functions
  virtual void create_atomic_centers_Gspace(Vector<std::complex<double>>& cG, Communicate& band_group_comm, int iorb) {}

  /** for exporting data from multi_UBspline_3d_d to multi_UBspline_3d_z
   *  This is only used by the legacy EinsplineSet class. To be deleted together with EinsplineSet.
   */
  void export_MultiSpline(multi_UBspline_3d_z** target)
  {
    *target                        = new multi_UBspline_3d_z;
    const auto* source_MultiSpline = (multi_UBspline_3d_d*)bspline->SplineInst->getSplinePtr();

    (*target)->spcode   = MULTI_U3D;
    (*target)->tcode    = DOUBLE_COMPLEX;
    (*target)->coefs    = (complex_double*)source_MultiSpline->coefs;
    (*target)->x_stride = source_MultiSpline->x_stride / 2;
    (*target)->y_stride = source_MultiSpline->y_stride / 2;
    (*target)->z_stride = source_MultiSpline->z_stride / 2;

    (*target)->x_grid.start     = source_MultiSpline->x_grid.start;
    (*target)->x_grid.end       = source_MultiSpline->x_grid.end;
    (*target)->x_grid.num       = source_MultiSpline->x_grid.num;
    (*target)->x_grid.delta     = source_MultiSpline->x_grid.delta;
    (*target)->x_grid.delta_inv = source_MultiSpline->x_grid.delta_inv;
    (*target)->y_grid.start     = source_MultiSpline->y_grid.start;
    (*target)->y_grid.end       = source_MultiSpline->y_grid.end;
    (*target)->y_grid.num       = source_MultiSpline->y_grid.num;
    (*target)->y_grid.delta     = source_MultiSpline->y_grid.delta;
    (*target)->y_grid.delta_inv = source_MultiSpline->y_grid.delta_inv;
    (*target)->z_grid.start     = source_MultiSpline->z_grid.start;
    (*target)->z_grid.end       = source_MultiSpline->z_grid.end;
    (*target)->z_grid.num       = source_MultiSpline->z_grid.num;
    (*target)->z_grid.delta     = source_MultiSpline->z_grid.delta;
    (*target)->z_grid.delta_inv = source_MultiSpline->z_grid.delta_inv;

    (*target)->xBC.lCode  = source_MultiSpline->xBC.lCode;
    (*target)->xBC.rCode  = source_MultiSpline->xBC.rCode;
    (*target)->xBC.lVal_r = source_MultiSpline->xBC.lVal;
    (*target)->xBC.lVal_i = source_MultiSpline->xBC.lVal;
    (*target)->xBC.rVal_r = source_MultiSpline->xBC.rVal;
    (*target)->xBC.rVal_i = source_MultiSpline->xBC.rVal;
    (*target)->yBC.lCode  = source_MultiSpline->yBC.lCode;
    (*target)->yBC.rCode  = source_MultiSpline->yBC.rCode;
    (*target)->yBC.lVal_r = source_MultiSpline->yBC.lVal;
    (*target)->yBC.lVal_i = source_MultiSpline->yBC.lVal;
    (*target)->yBC.rVal_r = source_MultiSpline->yBC.rVal;
    (*target)->yBC.rVal_i = source_MultiSpline->yBC.rVal;
    (*target)->zBC.lCode  = source_MultiSpline->zBC.lCode;
    (*target)->zBC.rCode  = source_MultiSpline->zBC.rCode;
    (*target)->zBC.lVal_r = source_MultiSpline->zBC.lVal;
    (*target)->zBC.lVal_i = source_MultiSpline->zBC.lVal;
    (*target)->zBC.rVal_r = source_MultiSpline->zBC.rVal;
    (*target)->zBC.rVal_i = source_MultiSpline->zBC.rVal;

    (*target)->num_splines = source_MultiSpline->num_splines / 2;
    (*target)->coefs_size  = source_MultiSpline->coefs_size / 2;
    // (*target)->lapl3 = (complex_double*) malloc (6*sizeof(double)*(*target)->z_stride);
  }

  /** for exporting data from multi_UBspline_3d_d to multi_UBspline_3d_z
   *  This is only used by the legacy EinsplineSet class. To be deleted together with EinsplineSet.
   */
  void export_MultiSpline(multi_UBspline_3d_d** target)
  {
    *target = (multi_UBspline_3d_d*)bspline->SplineInst->getSplinePtr();
  }

  std::unique_ptr<SPOSet> create_spline_set(int spin, const BandInfoGroup& bandgroup)
  {
    ReportEngine PRE("SplineSetReader", "create_spline_set(spin,SPE*)");
    //Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
    //double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;
    bspline = new splineset_t;
    app_log() << "  ClassName = " << bspline->getClassName() << std::endl;
    if (bspline->is_complex)
      app_log() << "  Using complex einspline table" << std::endl;
    else
      app_log() << "  Using real einspline table" << std::endl;

    // set info for Hybrid
    this->initialize_hybridrep_atomic_centers();

    //baseclass handles twists
    check_twists(bspline, bandgroup);

    Ugrid xyz_grid[3];

    typename splineset_t::BCType xyz_bc[3];
    bool havePsig = set_grid(bspline->HalfG, xyz_grid, xyz_bc);
    if (!havePsig)
      myComm->barrier_and_abort("SplineSetReader needs psi_g. Set precision=\"double\".");
    bspline->create_spline(xyz_grid, xyz_bc);

    std::ostringstream oo;
    oo << bandgroup.myName << ".g" << MeshSize[0] << "x" << MeshSize[1] << "x" << MeshSize[2] << ".h5";

    const std::string splinefile(oo.str());
    bool root       = (myComm->rank() == 0);
    int foundspline = 0;
    Timer now;
    if (root)
    {
      now.restart();
      hdf_archive h5f(myComm);
      foundspline = h5f.open(splinefile, H5F_ACC_RDONLY);
      if (foundspline)
      {
        std::string aname("none");
        foundspline = h5f.readEntry(aname, "class_name");
        foundspline = (aname.find(bspline->KeyWord) != std::string::npos);
      }
      if (foundspline)
      {
        int sizeD   = 0;
        foundspline = h5f.readEntry(sizeD, "sizeof");
        foundspline = (sizeD == sizeof(typename splineset_t::DataType));
      }
      if (foundspline)
      {
        foundspline = bspline->read_splines(h5f);
        if (foundspline)
          app_log() << "  Successfully restored coefficients from " << splinefile << ". The reading time is "
                    << now.elapsed() << " sec." << std::endl;
      }
      h5f.close();
    }
    myComm->bcast(foundspline);
    if (foundspline)
    {
      now.restart();
      bspline->bcast_tables(myComm);
      app_log() << "  SplineSetReader bcast the full table " << now.elapsed() << " sec." << std::endl;
      app_log().flush();
    }
    else
    {
      bspline->flush_zero();

      int nx = MeshSize[0];
      int ny = MeshSize[1];
      int nz = MeshSize[2];
      if (havePsig) //perform FFT using FFTW
      {
        FFTbox.resize(nx, ny, nz);
        FFTplan = fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(FFTbox.data()),
                                   reinterpret_cast<fftw_complex*>(FFTbox.data()), +1, FFTW_ESTIMATE);
        splineData_r.resize(nx, ny, nz);
        if (bspline->is_complex)
          splineData_i.resize(nx, ny, nz);

        TinyVector<double, 3> start(0.0);
        TinyVector<double, 3> end(1.0);
        spline_r = einspline::create(spline_r, start, end, MeshSize, bspline->HalfG);
        if (bspline->is_complex)
          spline_i = einspline::create(spline_i, start, end, MeshSize, bspline->HalfG);

        now.restart();
        initialize_spline_pio_gather(spin, bandgroup);
        app_log() << "  SplineSetReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;

        fftw_destroy_plan(FFTplan);
        FFTplan = NULL;
      }
      else //why, don't know
        initialize_spline_psi_r(spin, bandgroup);
      if (saveSplineCoefs && root)
      {
        now.restart();
        hdf_archive h5f;
        h5f.create(splinefile);
        std::string classname = bspline->getClassName();
        h5f.write(classname, "class_name");
        int sizeD = sizeof(typename splineset_t::DataType);
        h5f.write(sizeD, "sizeof");
        bspline->write_splines(h5f);
        h5f.close();
        app_log() << "  Stored spline coefficients in " << splinefile << " for potential reuse. The writing time is "
                  << now.elapsed() << " sec." << std::endl;
      }
    }

    clear();
    return std::unique_ptr<SPOSet>{bspline};
  }

  /** fft and spline cG
   * @param cG psi_g to be processed
   * @param ti twist index
   * @param iorb orbital index
   *
   * Perform FFT and spline to spline_r and spline_i
   */
  inline void fft_spline(Vector<std::complex<double>>& cG, int ti)
  {
    unpack4fftw(cG, mybuilder->Gvecs[0], MeshSize, FFTbox);
    fftw_execute(FFTplan);
    if (bspline->is_complex)
    {
      fix_phase_rotate_c2c(FFTbox, splineData_r, splineData_i, mybuilder->TwistAngles[ti], rotate_phase_r,
                           rotate_phase_i);
      einspline::set(spline_r, splineData_r.data());
      einspline::set(spline_i, splineData_i.data());
    }
    else
    {
      fix_phase_rotate_c2r(FFTbox, splineData_r, mybuilder->TwistAngles[ti], rotate_phase_r, rotate_phase_i);
      einspline::set(spline_r, splineData_r.data());
    }
  }


  /** initialize the splines
   */
  void initialize_spline_pio_gather(int spin, const BandInfoGroup& bandgroup)
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
    hdf_archive h5f(&band_group_comm, false);
    Vector<std::complex<double>> cG(mybuilder->Gvecs[0].size());
    const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
    if (band_group_comm.isGroupLeader())
      h5f.open(mybuilder->H5FileName, H5F_ACC_RDONLY);
    for (int iorb = iorb_first; iorb < iorb_last; iorb++)
    {
      if (band_group_comm.isGroupLeader())
      {
        int iorb_h5   = bspline->BandIndexMap[iorb];
        int ti        = cur_bands[iorb_h5].TwistIndex;
        std::string s = psi_g_path(ti, spin, cur_bands[iorb_h5].BandIndex);
        if (!h5f.readEntry(cG, s))
        {
          std::ostringstream msg;
          msg << "SplineSetReader Failed to read band(s) from h5 file. "
              << "Attemped dataset " << s << " with " << cG.size() << " complex numbers." << std::endl;
          throw std::runtime_error(msg.str());
        }
        double total_norm = compute_norm(cG);
        if ((checkNorm) && (std::abs(total_norm - 1.0) > PW_COEFF_NORM_TOLERANCE))
        {
          std::ostringstream msg;
          msg << "SplineSetReader The orbital " << iorb_h5 << " has a wrong norm " << total_norm
              << ", computed from plane wave coefficients!" << std::endl
              << "This may indicate a problem with the HDF5 library versions used "
              << "during wavefunction conversion or read." << std::endl;
          throw std::runtime_error(msg.str());
        }
        fft_spline(cG, ti);
        bspline->set_spline(spline_r, spline_i, cur_bands[iorb_h5].TwistIndex, iorb, 0);
      }
      this->create_atomic_centers_Gspace(cG, band_group_comm, iorb);
    }

    myComm->barrier();
    Timer now;
    if (band_group_comm.isGroupLeader())
    {
      now.restart();
      bspline->gather_tables(band_group_comm.GroupLeaderComm);
      app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
    }
    now.restart();
    bspline->bcast_tables(myComm);
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
  }

  void initialize_spline_psi_r(int spin, const BandInfoGroup& bandgroup)
  {
    //not used by may be enabled later
    myComm->barrier_and_abort("SplineSetReaderP initialize_spline_psi_r implementation not finished.");
    int nx = MeshSize[0];
    int ny = MeshSize[1];
    int nz = MeshSize[2];
    splineData_r.resize(nx, ny, nz);
    splineData_i.resize(ny, ny, nz);
    Array<std::complex<double>, 3> rawData(nx, ny, nz);
    const std::vector<BandInfo>& cur_bands = bandgroup.myBands;

    // create single splines as an intermediate storage
    TinyVector<double, 3> start(0.0);
    TinyVector<double, 3> end(1.0);
    spline_r = einspline::create(spline_r, start, end, MeshSize, bspline->HalfG);
    if (bspline->is_complex)
      spline_i = einspline::create(spline_i, start, end, MeshSize, bspline->HalfG);

    //this will be parallelized with OpenMP
    int N = bandgroup.getNumDistinctOrbitals();
    for (int iorb = 0; iorb < N; ++iorb)
    {
      //check dimension
      if (myComm->rank() == 0)
      {
        std::string path = psi_r_path(cur_bands[iorb].TwistIndex, spin, cur_bands[iorb].BandIndex);
        HDFAttribIO<Array<std::complex<double>, 3>> h_splineData(rawData);
        h_splineData.read(mybuilder->H5FileID, path.c_str());
        simd::copy(splineData_r.data(), splineData_i.data(), rawData.data(), rawData.size());
      }
      mpi::bcast(*myComm, splineData_r);
      einspline::set(spline_r, splineData_r.data());
      if (bspline->is_complex)
      {
        mpi::bcast(*myComm, splineData_i);
        einspline::set(spline_i, splineData_i.data());
      }
      bspline->set_spline(spline_r, spline_i, cur_bands[iorb].TwistIndex, iorb, 0);
    }
  }
};
} // namespace qmcplusplus
#endif
