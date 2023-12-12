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
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "OneSplineOrbData.hpp"
#include <spline/einspline_engine.hpp>
#include "einspline_helper.hpp"

namespace qmcplusplus
{
void OneSplineOrbData::create(const TinyVector<int, 3>& halfG)
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

  OneSplineOrbData::OneSplineOrbData(const TinyVector<int, 3>& mesh_size, const TinyVector<int, 3>& halfG, const bool isComplex)
      : mesh_size_(mesh_size), isComplex_(isComplex)
  {
    create(halfG);
  }

  OneSplineOrbData::~OneSplineOrbData() { clear(); }

  void OneSplineOrbData::clear()
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
  void OneSplineOrbData::fft_spline(const Vector<std::complex<double>>& cG,
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
} // namespace qmcplusplus
