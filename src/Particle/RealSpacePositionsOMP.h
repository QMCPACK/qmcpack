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


/** @file RealSpacePostionsOffload.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_OFFLOAD_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_OFFLOAD_H

#include "Particle/QuantumVariables.h"
#include "OhmmsSoA/Container.h"
#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
class RealSpacePositionsOMP : public QuantumVariables
{
public:
  RealSpacePositionsOMP() : QuantumVariables(QuantumVariableKind::QV_POS_OFFLOAD), RSoA_device_ptr(nullptr) {}
  RealSpacePositionsOMP(const RealSpacePositionsOMP& in)
    : QuantumVariables(QuantumVariableKind::QV_POS_OFFLOAD), RSoA(in.RSoA)
  {
    RSoA_hostview.attachReference(RSoA.size(), RSoA.capacity(), RSoA.data());
    auto* pos_ptr = RSoA.data();
    PRAGMA_OFFLOAD("#pragma omp target data use_device_ptr(pos_ptr)")
    {
      RSoA_device_ptr = pos_ptr;
    }
    updateH2D();
  }

  std::unique_ptr<QuantumVariables> makeClone() override { return std::make_unique<RealSpacePositionsOMP>(*this); }

  void resize(size_t n) override
  {
    if (RSoA.size() != n)
    {
      RSoA.resize(n);
      RSoA_hostview.attachReference(RSoA.size(), RSoA.capacity(), RSoA.data());
      auto* pos_ptr = RSoA.data();
      PRAGMA_OFFLOAD("#pragma omp target data use_device_ptr(pos_ptr)")
      {
        RSoA_device_ptr = pos_ptr;
      }
    }
  }

  size_t size() override { return RSoA_hostview.size(); }

  void setAllParticlePos(const ParticlePos_t& R) override
  {
    resize(R.size());
    RSoA_hostview.copyIn(R);
    updateH2D();
  }

  void setOneParticlePos(const PosType& pos, size_t iat) override
  {
    RSoA_hostview(iat) = pos;
    /* This was too slow due to overhead.
    RealType x     = pos[0];
    RealType y     = pos[1];
    RealType z     = pos[2];
    RealType* data = RSoA.data();
    size_t offset  = RSoA.capacity();

    PRAGMA_OFFLOAD("omp target map(to : x, y, z, iat)")
    {
      data[iat]              = x;
      data[iat + offset]     = y;
      data[iat + offset * 2] = z;
    }
    */
  }

  const PosVectorSoa& getAllParticlePos() override { return RSoA_hostview; }
  PosType getOneParticlePos(size_t iat) const override { return RSoA_hostview[iat]; }

  void donePbyP() override { updateH2D(); }

  RealType* getDevicePtr() const { return RSoA_device_ptr; }
private:
  ///particle positions in SoA layout
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> RSoA;

  /// the pointer of RSoA memory on the device
  RealType* RSoA_device_ptr;

  ///host view of RSoA
  PosVectorSoa RSoA_hostview;

  void updateH2D()
  {
    RealType* data = RSoA.data();
    PRAGMA_OFFLOAD("omp target update to(data[0:RSoA.capacity()*QMCTraits::DIM])")
  }
};
} // namespace qmcplusplus
#endif
