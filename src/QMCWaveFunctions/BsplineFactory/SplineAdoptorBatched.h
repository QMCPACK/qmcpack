//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** \file SplineAdoptorBatched.h
 * \class SplineAdoptorBatched
 * Base class for SplineAdoptor's used for BsplineSet<SplineAdoptor>
 * Specifies that a SplineXXXAdopter provides these functions
 * - evaluate_v    value only
 * - evaluate_vgl  vgl
 * - evaluate_vgh  vgh
 * Specializations are implemented  in Spline*Adoptor.h and include
 * - SplineC2RAdoptor<ST,TT,D> : real wavefunction using complex einspline, tiling
 * - SplineC2CAdoptor<ST,TT,D> : complex wavefunction using complex einspline, tiling
 * - SplineR2RAdoptor<ST,TT,D> : real wavefunction using real einspline, a single twist
 * where ST (TT) is the precision of the einspline (SPOSetBase).
 *
 * typedefs and data members are duplicated for each adoptor class.
 * @todo Specalization and optimization for orthorhombic cells to use vgl not vgh
 */
#ifndef QMCPLUSPLUS_SPLINEADOPTORBASEVECTORIZED_H
#define QMCPLUSPLUS_SPLINEADOPTORBASEVECTORIZED_H

//#include "QMCWaveFunctions/BsplineFactory/BsplineDevice.h"
//#include "QMCWaveFunctions/BsplineFactory/BsplineDeviceCUDA.h"
#include "Lattice/CrystalLattice.h"
#include "simd/allocator.hpp"
#include "SplineAdoptor.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

/** base class any SplineAdoptor
 *
 * This handles SC and twist and declare storage for einspline
 */
template<template<typename, unsigned> class DEVICE, typename ST, unsigned D>
class SplineAdoptorBatched : public SplineAdoptor<ST, D>
{
public:
  using BaseType = SplineAdoptor<ST, D>;
  using PointType = TinyVector<ST,D>;
  using SingleSplineType=UBspline_3d_d;
  using DataType = ST; 

  DEVICE<ST,D> bspline_dev;


  //static_assert(std::is_base_of<BsplineDevice<DEVICE, ST, D>, DEVICE>, "DEVICE must inherit from BsplineDevice");

  ///true if the computed values are complex
  SplineAdoptorBatched() { };
  SplineAdoptorBatched(const SplineAdoptorBatched& rhs)=default;

  // inline void resizeStorage(size_t n, size_t nvals)
  // {
  //   init_base(n);
  // }
  
  inline void init_base(int n)
  {
    BaseType::init_base(n);
  }

  // template<typename VV>
  // inline void evaluate_v(const std::vector<PointType>& P, const int iat, VV& psi)
  // {
  //   PointType ru(PrimLattice.toUnit_floor(P));
  //   bspline_dev.evaluate(ru,myV);
  //   bspline_dev.assign_v(r,myV,psi);
  // }

  ///remap kpoints to group general kpoints & special kpoints
  int remap_kpoints()
  {
    return BaseType::remap_kpoints();
  }

  void evaluate(const ParticleSet& P, int iat, typename SPOSet::ValueVector_t& psi)
  {
    APP_ABORT("SplineAdoptorBatched doesn't implement single walker evaluates");
  }
};

}
#endif
