//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
