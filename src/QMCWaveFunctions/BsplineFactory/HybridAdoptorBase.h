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

#include <Particle/DistanceTableData.h>
#include <QMCWaveFunctions/lcao/SoaSphericalTensor.h>

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
  using value_type=ST;

  ST cutoff;
  int lmax;
  SoaSphericalTensor<ST> Ylm;
  AtomicSplineType* MultiSpline;

  AtomicOrbitalSoA(int Lmax):Ylm(Lmax), MultiSpline(nullptr), lmax(Lmax) { }

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
    //cutoff, spline radius, spline points, check coeff size
  }

  bool write_splines(hdf_archive& h5f)
  {
    //dump center info, including cutoff, spline radius, spline points, lm
    //name, position, consistency
    //dump spline coefficients    //dump spline coefficients
  }

  template<typename VV>
  inline void evaluate_v(const PointType& r, VV& myV)
  {
    //evaluate only V
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const PointType& r, VV& myV, GV& myG, VV& myH)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const PointType& r, VV& myV, GV& myG, HT& myH)
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
  using PointType=typename AtomicOrbitalSoA<ST>::PointType;

  std::vector<AtomicOrbitalSoA<ST> > AtomicCenters;
  // I-e distance table
  DistanceTableData* ei_dist;
  //mapping supercell to primitive cell
  std::vector<int> Super2Prim;

  HybridAdoptorBase() { }

  bool read_splines(hdf_archive& h5f)
  {
    // for each center
    // loop load center info, including name, position, check consistency
    // Initial class with index
    // call read_splines
    // push to the vector
  }

  bool write_splines(hdf_archive& h5f)
  {
    //loop for each center
    // write_splines
  }

  template<typename VV>
  inline bool evaluate_v(VV& myV)
  {
    //evaluate only V
    bool inAtom=false;
    const int center_idx=ei_dist->get_first_neighbor_temporal();
    if(center_idx<0) abort();
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( ei_dist->Temp_r[center_idx] < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_v(ei_dist->Temp_dr[center_idx], myV);
    }
    return inAtom;
  }

  template<typename VV, typename GV>
  inline bool evaluate_vgl(VV& myV, GV& myG, VV& myH)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  inline bool evaluate_vgh(VV& myV, GV& myG, HT& myH)
  {
    bool inAtom=false;
    const int center_idx=ei_dist->get_first_neighbor_temporal();
    if(center_idx<0) abort();
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( ei_dist->Temp_r[center_idx] < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgh(ei_dist->Temp_dr[center_idx], myV, myG, myH);
    }
    return inAtom;
  }
};

}
#endif
