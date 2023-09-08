//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/** @file DynamicCoordinatesT.h
 */
#ifndef QMCPLUSPLUS_DYNAMICCOORDINATEST_H
#define QMCPLUSPLUS_DYNAMICCOORDINATEST_H

#include <memory>

#include "OhmmsSoA/VectorSoaContainer.h"
#include "ParticleSetTraits.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
class ResourceCollection;

/** enumerator for DynamicCoordinates kinds
 */
enum class DynamicCoordinateKind
{
  DC_POS,         // SoA positions
  DC_POS_OFFLOAD, // SoA positions with OpenMP offload
};

/** quantum variables of all the particles
 */
template <typename T>
class DynamicCoordinatesT
{
public:
    using RealType = typename ParticleSetTraits<T>::RealType;
    using PosType = typename ParticleSetTraits<T>::PosType;
    using ParticlePos = typename LatticeParticleTraits<T>::ParticlePos;
    using PosVectorSoa =
        VectorSoaContainer<RealType, ParticleSetTraits<T>::DIM>;

    DynamicCoordinatesT(const DynamicCoordinateKind kind_in) :
        variable_kind_(kind_in)
    {
    }

    DynamicCoordinatesT(const DynamicCoordinatesT&) = default;
    DynamicCoordinatesT&
    operator=(const DynamicCoordinatesT&) = delete;

    DynamicCoordinateKind
    getKind() const
    {
        return variable_kind_;
    }

    virtual ~DynamicCoordinatesT() = default;

    virtual std::unique_ptr<DynamicCoordinatesT>
    makeClone() = 0;

    /** resize internal storages based on the number of particles
     *  @param n the number of particles
     */
    virtual void
    resize(size_t n) = 0;
    /// return the number of particles
    virtual size_t
    size() const = 0;

    /// overwrite the positions of all the particles.
    virtual void
    setAllParticlePos(const ParticlePos& R) = 0;
    /// overwrite the position of one the particle.
    virtual void
    setOneParticlePos(const PosType& pos, size_t iat) = 0;
    /** copy the active positions of particles with a uniform id in all the
     * walkers to a single internal buffer.
     *  @param coords_list a batch of DynamicCoordinates
     *  @param iat paricle id, uniform across coords_list
     *  @param new_positions proposed positions
     */
    virtual void
    mw_copyActivePos(
        const RefVectorWithLeader<DynamicCoordinatesT>& coords_list, size_t iat,
        const std::vector<PosType>& new_positions) const
    {
        assert(this == &coords_list.getLeader());
    }

    /** overwrite the positions of particles with a uniform id in all the
     * walkers upon acceptance.
     *  @param coords_list a batch of DynamicCoordinates
     *  @param iat paricle id, uniform across coords_list
     *  @param new_positions proposed positions
     *  @param isAccepted accept/reject info
     */
    virtual void
    mw_acceptParticlePos(
        const RefVectorWithLeader<DynamicCoordinatesT>& coords_list, size_t iat,
        const std::vector<PosType>& new_positions,
        const std::vector<bool>& isAccepted) const = 0;

    /// all particle position accessor
    virtual const PosVectorSoa&
    getAllParticlePos() const = 0;
    /// one particle position accessor
    virtual PosType
    getOneParticlePos(size_t iat) const = 0;

    /// secure internal data consistency after p-by-p moves
    virtual void
    donePbyP()
    {
    }

    /// initialize a shared resource and hand it to a collection
    virtual void
    createResource(ResourceCollection& collection) const
    {
    }

    /// acquire a shared resource from a collection
    virtual void
    acquireResource(ResourceCollection& collection,
        const RefVectorWithLeader<DynamicCoordinatesT>& coords_list) const
    {
    }

    /// return a shared resource to a collection
    virtual void
    releaseResource(ResourceCollection& collection,
        const RefVectorWithLeader<DynamicCoordinatesT>& coords_list) const
    {
    }

protected:
    /// type of dynamic coordinates
    const DynamicCoordinateKind variable_kind_;
};

/** create DynamicCoordinates based on kind
 */
template <typename T>
std::unique_ptr<DynamicCoordinatesT<T>> createDynamicCoordinatesT(
    const DynamicCoordinateKind kind = DynamicCoordinateKind::DC_POS);
} // namespace qmcplusplus
#endif
