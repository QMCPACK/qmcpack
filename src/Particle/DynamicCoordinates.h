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
#include "OhmmsSoA/Container.h"

namespace qmcplusplus
{
/** enumerator for DynamicCoordinates kinds
 */
enum class DynamicCoordinateKind
{
  DC_POS, // SoA positions
  DC_POS_OFFLOAD, // SoA positions with OpenMP offload
};

/** quantum variables of all the particles
 */
class DynamicCoordinates
{
public:
  using RealType = QMCTraits::RealType;
  using PosType = QMCTraits::PosType;
  using ParticlePos_t = PtclOnLatticeTraits::ParticlePos_t;
  using PosVectorSoa = VectorSoaContainer<RealType, QMCTraits::DIM>;

  DynamicCoordinates(const DynamicCoordinateKind kind_in) : variable_kind_(kind_in) {}

  DynamicCoordinates(const DynamicCoordinates&) = default;
  DynamicCoordinates& operator=(const DynamicCoordinates&) = delete;

  DynamicCoordinateKind getKind() const { return variable_kind_; }

  virtual ~DynamicCoordinates() = default;

  virtual std::unique_ptr<DynamicCoordinates> makeClone() = 0;

  virtual void resize(size_t n) = 0;
  virtual size_t size() = 0;

  virtual void setAllParticlePos(const ParticlePos_t& R) = 0;
  virtual void setOneParticlePos(const PosType& pos, size_t iat) = 0;

  virtual const PosVectorSoa& getAllParticlePos() const = 0;
  virtual PosType getOneParticlePos(size_t iat) const = 0;

  virtual void donePbyP() { }
protected:
  const DynamicCoordinateKind variable_kind_;
};
} // namespace qmcplusplus
#endif
