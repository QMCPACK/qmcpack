//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BSPLINE_AO_BUILDER_H
#define QMCPLUSPLUS_BSPLINE_AO_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"

namespace qmcplusplus
{

class BsplineAOBuilder: public QMCTraits
{

public:

  typedef BsplineFunctor<RealType>                 RadialOrbitalType;
  typedef SphericalBasisSet<RadialOrbitalType>     CenteredOrbitalType;

  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///the species
  std::string m_species;
  ///constructor
  BsplineAOBuilder(xmlNodePtr cur=NULL);

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
  {
    m_orbitals = oset;
    m_species = acenter;
  }

  bool addGrid(xmlNodePtr cur)
  {
    return true;
  }

  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms, bool useSphericalHarmonicsNormalization=true);

  bool putCommon(xmlNodePtr cur);

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1253 $   $Date: 2006-08-13 10:15:48 -0500 (Sun, 13 Aug 2006) $
 * $Id: BsplineAOBuilder.h 1253 2006-08-13 15:15:48Z jnkim $
 ***************************************************************************/
