//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DISTANCETABLE_H
#define QMCPLUSPLUS_DISTANCETABLE_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus
{
/** Class to manage multiple DistanceTableData objects.
 *
 * \date  2008-09-19
 * static data members are removed. DistanceTable::add functions
 * are kept for compatibility only. New codes should use a member function
 * of ParticleSet to add a distance table
 * int ParticleSet::addTable(const ParticleSet& source)
 *
 * \deprecated There is only one instance of the data memebers of
 * DistanceTable in an application and the data are shared by many objects.
 * Note that static data members and functions are used
 * (based on singleton and factory patterns).
 *\todo DistanceTable should work as a factory, as well, to instantiate DistanceTableData
 * subject to different boundary conditions.
 * Lattice/CrystalLattice.h and Lattice/CrystalLattice.cpp can be owned by DistanceTable
 * to generically control the crystalline structure.
 */

///free function to create a distable table of s-s
DistanceTableData* createDistanceTable(ParticleSet& s, std::ostream& description);

///free function create a distable table of s-t
DistanceTableData* createDistanceTableAB(const ParticleSet& s, ParticleSet& t, std::ostream& description);
DistanceTableData* createDistanceTableABOMP(const ParticleSet& s, ParticleSet& t, std::ostream& description);

inline DistanceTableData* createDistanceTable(const ParticleSet& s, ParticleSet& t, std::ostream& description)
{
  // during P-by-P move, the cost of single particle evaluation of distance tables
  // is determined by the number of source particles.
  // Thus the implementation selection is determined by the source particle set.
#if defined(ENABLE_OFFLOAD)
  if (s.getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD)
    return createDistanceTableABOMP(s, t, description);
  else
#endif
    return createDistanceTableAB(s, t, description);
}

} // namespace qmcplusplus
#endif
