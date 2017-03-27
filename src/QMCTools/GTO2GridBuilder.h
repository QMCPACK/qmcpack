//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIAN_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIAN_H

#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCTools/MolecularOrbitalBasis.h"
#include "QMCTools/RGFBuilderBase.h"

namespace qmcplusplus
{

/**Class to convert GaussianTypeOrbital to a radial orbital on a log grid.
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 */
struct GTO2GridBuilder: public RGFBuilderBase
{
  ///Boolean
  bool Normalized;
  ///constructor
  GTO2GridBuilder(bool normalized=false):Normalized(normalized) {}
  //bool addGrid(xmlNodePtr cur);
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
  bool putCommon(xmlNodePtr cur);

  bool addGrid(xmlNodePtr cur);
};
}
#endif
