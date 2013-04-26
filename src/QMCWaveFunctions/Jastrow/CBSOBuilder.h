//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SPHERICALBASISBUILDER_H
#define QMCPLUSPLUS_SPHERICALBASISBUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"

namespace qmcplusplus
{

/** Builder of cubic-bspline radial orbitals associated with a center
 *
 * CBSO stands for CubicBSplineOrbitals
 * For a center,
 *   - only one linear grid is used
 *   - any number of radial orbitals  of CubicBsplineSingle
 * Will test the group function later using CubicBsplineGroup<T,GRIDTYPE>
 */
class CBSOBuilder: public QMCTraits
{
  //typedef CubicBsplineGroup<RealType,LINEAR_1DGRID> RadialOrbitalType;
public:
  ///spline engine
  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
  ///numerical functor
  typedef CubicSplineSingle<RealType,SplineEngineType> RadialOrbitalType;
  ///spherical basis set using splines
  typedef SphericalBasisSet<RadialOrbitalType>         CenteredOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///number of grid points [0,Npts]
  int NumGridPoints;
  ///maximum cutoff
  RealType m_rcut;
  ///the quantum number of this node
  QuantumNumberType m_nlms;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///the species
  std::string m_species;

  ///constructor
  CBSOBuilder(xmlNodePtr cur=NULL);

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter);

  ///add a grid
  bool addGrid(xmlNodePtr cur);

  /** add a radial functor
   * @param cur xml element
   * @param nlms quantum number
   */
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

  /** put common element
   * @param cur xml element
   */
  bool putCommon(xmlNodePtr cur);

private:
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1252 $   $Date: 2006-08-13 10:15:01 -0500 (Sun, 13 Aug 2006) $
 * $Id: CBSOBuilder.h 1252 2006-08-13 15:15:01Z jnkim $
 ***************************************************************************/
