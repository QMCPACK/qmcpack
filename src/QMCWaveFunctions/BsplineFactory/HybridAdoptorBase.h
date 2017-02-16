//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory 
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory

//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file HybridAdoptorBase.h
 *
 * Hybrid adoptor base class
 */
#ifndef QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H
#define QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H

//#include <QMCWaveFunctions/lcao/SoaSphericalTensor.h>

namespace qmcplusplus
{

template<typename ST>
struct AtomicOrbitalSoA
{
  static const int D=3;
  using AtomicSplineType=typename bspline_traits<ST,1>::SplineType;
  using AtomicBCType=typename bspline_traits<ST,1>::BCType;
  using AtomicSingleSplineType=UBspline_1d_d;
  using PointType=TinyVector<ST,D>;

  ST cutoff;
  //SoaSphericalTensor<ST> Ylm;
  AtomicSplineType* MultiSpline;

  AtomicOrbitalSoA():MultiSpline(nullptr)
  {
  }

  //~AtomicOrbitalSoA();

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
  }

  inline void set_spline(AtomicSingleSplineType* spline_r, AtomicSingleSplineType* spline_i, int twist, int ispline, int level)
  {
  }

  void set_spline(ST* restrict psi_r, ST* restrict psi_i, int twist, int ispline, int level)
  {
  }

  bool read_splines(hdf_archive& h5f)
  {
    //load spline coefficients
  }

  bool write_splines(hdf_archive& h5f)
  {
    //dump spline coefficients
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    //evaluate only V
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    //missing
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const PointType& r, VGL& vgl)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, HT& hess)
  {
    //Needed to do tensor product here
  }
};

/** adoptor class to match 
 *
 */
template<typename ST>
struct HybridAdoptorBase
{
  static const int D=3;
  using PointType=TinyVector<ST,D>;

  std::vector<AtomicOrbitalSoA<ST> > AtomicCenters;
  //mapping supercell to primitive cell

  HybridAdoptorBase() { }

  bool read_splines(hdf_archive& h5f)
  {
    //load spline coefficients
    //for each center
  }

  bool write_splines(hdf_archive& h5f)
  {
    //dump spline coefficients
    //for each center
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& psi)
  {
    //evaluate only V
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    //missing
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const PointType& r, VGL& vgl)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, HT& hess)
  {
    //Needed to do tensor product here
  }
};

}
#endif
