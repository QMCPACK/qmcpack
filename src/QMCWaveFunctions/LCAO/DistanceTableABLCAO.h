//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Anouar Benali, abenali.sci@hotmail.com, Qubit Pharmaceuticals.
//
// File created by:  Anouar Benali, abenali.sci@hotmail.com, Qubit Pharmaceuticals.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DISTANCETABLE_AB_LCAO_H
#define QMCPLUSPLUS_DISTANCETABLE_AB_LCAO_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{

/** Specialized distance table for LCAO basis evaluation on GPU
 * 
 * This class provides GPU-optimized computation of ion-electron displacements
 * for LCAO basis set evaluation. It avoids CPU-GPU transfers by working directly
 * with device pointers and the fused new position buffer.
 */
class DistanceTableABLCAO : public DistanceTable
{
public:
  using RealType  = OHMMS_PRECISION;
  using PosType   = TinyVector<RealType, OHMMS_DIM>;
  using IndexType = int;

  template<typename DT>
  using OffloadPinnedVector = Vector<DT, OffloadPinnedAllocator<DT>>;

private:
  int num_ions_; ///< Number of ion centers
  int num_elec_; ///< Number of electrons


public:
  /** Constructor
   * @param ions ion particle set
   * @param elecs electron particle set
   */
  DistanceTableABLCAO(const ParticleSet& ions, const ParticleSet& elecs);
  ~DistanceTableABLCAO() = default;

  // Required pure virtual methods from DistanceTable (unused stubs for LCAO)
  void evaluate(ParticleSet& P) override;
  void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old = true) override;
  void update(IndexType jat) override;
  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override;

  /** Multi-walker evaluation for LCAO basis
   * @param elec_list multi-walker electron particle sets
   * @param ions ion particle set
   * @param iat target electron index
   * @param num_centers number of ion centers to evaluate
   * @param displ_list_tr output displacement vectors [dim + 3*(walker + center*nw)]
   * @param Tv_list output lattice translation vectors (for PBC)
   * 
   * Computes ion-electron displacements using fused new position buffer.
   * For PBC systems, applies minimum image convention and computes lattice translations.
   */
  void mw_evaluate(const RefVectorWithLeader<ParticleSet>& elec_list,
                   const ParticleSet& ions,
                   int iat,
                   int num_centers,
                   OffloadPinnedVector<RealType>& displ_list_tr,
                   OffloadPinnedVector<RealType>& Tv_list);
};

} // namespace qmcplusplus
#endif
