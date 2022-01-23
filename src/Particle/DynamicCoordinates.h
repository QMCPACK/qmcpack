//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file DynamicCoordinates.h
 */
#ifndef QMCPLUSPLUS_DYNAMICCOORDINATES_H
#define QMCPLUSPLUS_DYNAMICCOORDINATES_H

#include <memory>
#include "Configuration.h"
#include "OhmmsSoA/VectorSoaContainer.h"
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
class DynamicCoordinates
{
public:
  using RealType     = QMCTraits::RealType;
  using PosType      = QMCTraits::PosType;
  using ParticlePos  = PtclOnLatticeTraits::ParticlePos;
  using PosVectorSoa = VectorSoaContainer<RealType, QMCTraits::DIM>;

  DynamicCoordinates(const DynamicCoordinateKind kind_in) : variable_kind_(kind_in) {}

  DynamicCoordinates(const DynamicCoordinates&)            = default;
  DynamicCoordinates& operator=(const DynamicCoordinates&) = delete;

  DynamicCoordinateKind getKind() const { return variable_kind_; }

  virtual ~DynamicCoordinates() = default;

  virtual std::unique_ptr<DynamicCoordinates> makeClone() = 0;

  /** resize internal storages based on the number of particles
   *  @param n the number of particles
   */
  virtual void resize(size_t n) = 0;
  /// return the number of particles
  virtual size_t size() const = 0;

  /// overwrite the positions of all the particles.
  virtual void setAllParticlePos(const ParticlePos& R) = 0;
  /// overwrite the position of one the particle.
  virtual void setOneParticlePos(const PosType& pos, size_t iat) = 0;
  /** copy the active positions of particles with a uniform id in all the walkers to a single internal buffer.
   *  @param coords_list a batch of DynamicCoordinates
   *  @param iat paricle id, uniform across coords_list
   *  @param new_positions proposed positions
   */
  virtual void mw_copyActivePos(const RefVectorWithLeader<DynamicCoordinates>& coords_list,
                                size_t iat,
                                const std::vector<PosType>& new_positions) const
  {
    assert(this == &coords_list.getLeader());
  }

  /** overwrite the positions of particles with a uniform id in all the walkers upon acceptance.
   *  @param coords_list a batch of DynamicCoordinates
   *  @param iat paricle id, uniform across coords_list
   *  @param new_positions proposed positions
   *  @param isAccepted accept/reject info
   */
  virtual void mw_acceptParticlePos(const RefVectorWithLeader<DynamicCoordinates>& coords_list,
                                    size_t iat,
                                    const std::vector<PosType>& new_positions,
                                    const std::vector<bool>& isAccepted) const = 0;

  /// all particle position accessor
  virtual const PosVectorSoa& getAllParticlePos() const = 0;
  /// one particle position accessor
  virtual PosType getOneParticlePos(size_t iat) const = 0;

  /// secure internal data consistency after p-by-p moves
  virtual void donePbyP() {}

  /// initialize a shared resource and hand it to a collection
  virtual void createResource(ResourceCollection& collection) const {}

  /// acquire a shared resource from a collection
  virtual void acquireResource(ResourceCollection& collection,
                               const RefVectorWithLeader<DynamicCoordinates>& coords_list) const
  {}

  /// return a shared resource to a collection
  virtual void releaseResource(ResourceCollection& collection,
                               const RefVectorWithLeader<DynamicCoordinates>& coords_list) const
  {}

protected:
  /// type of dynamic coordinates
  const DynamicCoordinateKind variable_kind_;
};
} // namespace qmcplusplus
#endif
