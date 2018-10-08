//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineDeviceCPU.h
 *
 * \class BsplineDeviceCPU
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
#ifndef QMCPLUSPLUS_BSPLINEDEVICECPU_H
#define QMCPLUSPLUS_BSPLINEDEVICECPU_H

#include <iostream>
#include "QMCWaveFunctions/BsplineFactory/BsplineDevice.h"
#include "einspline/bspline_base.h"
#include "einspline/bspline_structs.h"
namespace qmcplusplus
{

/** base class any SplineAdoptor
 *
 * This handles SC and twist and declare storage for einspline
 */

template<typename ST, unsigned D>
class BsplineDeviceCPU : BsplineDevice<BsplineDeviceCPU<ST,D>, ST, D>
{
public:
  using SingleBsplineType=UBspline_3d_d;
  void implementation()
  {
    std::cout<< "implemented for CUDA\n";
  }
};

}
#endif
