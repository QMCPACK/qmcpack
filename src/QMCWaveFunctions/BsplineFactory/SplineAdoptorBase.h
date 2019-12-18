//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file SplineAdoptorBase.h
 *
 * Base class for SplineAdoptor's used for BsplineSet<SplineAdoptor>
 * Specifies that a SplineXXXAdopter provides these functions
 * - evaluate_v    value only
 * - evaluate_vgl  vgl
 * - evaluate_vgh  vgh
 * Specializations are implemented  in Spline*Adoptor.h and include
 * - SplineC2RAdoptor<ST,TT,D> : real wavefunction using complex einspline, tiling
 * - SplineC2CAdoptor<ST,TT,D> : complex wavefunction using complex einspline, tiling
 * - SplineR2RAdoptor<ST,TT,D> : real wavefunction using real einspline, a single twist
 * where ST (TT) is the precision of the einspline (SPOSet).
 *
 * typedefs and data members are duplicated for each adoptor class.
 * @todo Specalization and optimization for orthorhombic cells to use vgl not vgh
 */
#ifndef QMCPLUSPLUS_SPLINEADOPTORBASE_H
#define QMCPLUSPLUS_SPLINEADOPTORBASE_H

namespace qmcplusplus
{
/** base class any SplineAdoptor
 *
 * This handles SC and twist and declare storage for einspline
 */
template<typename ST, unsigned D>
struct SplineAdoptorBase
{
  using PointType        = TinyVector<ST, D>;
  using SingleSplineType = UBspline_3d_d;
  using DataType         = ST;

  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, D> GGt;
};

} // namespace qmcplusplus
#endif
