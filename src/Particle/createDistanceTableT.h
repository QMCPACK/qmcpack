//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of
// Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DISTANCETABLET_H
#define QMCPLUSPLUS_DISTANCETABLET_H

#include "Particle/ParticleSetT.h"

namespace qmcplusplus
{
/** Class to manage multiple DistanceTable objects.
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
 *\todo DistanceTable should work as a factory, as well, to instantiate
 *DistanceTable subject to different boundary conditions.
 * Lattice/CrystalLattice.h and Lattice/CrystalLattice.cpp can be owned by
 *DistanceTable to generically control the crystalline structure.
 */

/// free function to create a distable table of s-s
template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableAAT(ParticleSetT<T>& s, std::ostream& description);

template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableAATOMPTarget(ParticleSetT<T>& s, std::ostream& description);

template <typename T>
inline std::unique_ptr<DistanceTableT<T>>
createDistanceTableT(ParticleSetT<T>& s, std::ostream& description)
{
    // during P-by-P move, the cost of single particle evaluation of distance
    // tables is determined by the number of source particles. Thus the
    // implementation selection is determined by the source particle set.
    if (s.getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD)
        return createDistanceTableAATOMPTarget(s, description);
    else
        return createDistanceTableAAT(s, description);
}

/// free function create a distable table of s-t
template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableABT(
    const ParticleSetT<T>& s, ParticleSetT<T>& t, std::ostream& description);

template <typename T>
std::unique_ptr<DistanceTableT<T>>
createDistanceTableABTOMPTarget(
    const ParticleSetT<T>& s, ParticleSetT<T>& t, std::ostream& description);

template <typename T>
inline std::unique_ptr<DistanceTableT<T>>
createDistanceTableT(
    const ParticleSetT<T>& s, ParticleSetT<T>& t, std::ostream& description)
{
    // during P-by-P move, the cost of single particle evaluation of distance
    // tables is determined by the number of source particles. Thus the
    // implementation selection is determined by the source particle set.
    if (s.getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD)
        return createDistanceTableABTOMPTarget(s, t, description);
    else
        return createDistanceTableABT(s, t, description);
}

} // namespace qmcplusplus
#endif
