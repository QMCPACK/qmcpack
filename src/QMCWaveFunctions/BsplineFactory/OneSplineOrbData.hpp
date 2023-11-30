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


#ifndef QMCPLUSPLUS_ONESPLINEORBDATA_H
#define QMCPLUSPLUS_ONESPLINEORBDATA_H

#include <complex>
#include <fftw3.h>
#include <einspline/bspline_structs.h>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsArray.h>

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

  void create(const TinyVector<int, 3>& halfG);

public:
  OneSplineOrbData(const TinyVector<int, 3>& mesh_size, const TinyVector<int, 3>& halfG, const bool isComplex);

  ~OneSplineOrbData();

  auto getRotatePhase() const { return std::complex<double>(rotate_phase_r, rotate_phase_i); }
  auto& get_spline_r() { return *spline_r; }
  auto& get_spline_i() { return *spline_i; }

  void clear();

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
                  const bool rotate);
};
} // namespace qmcplusplus
#endif
