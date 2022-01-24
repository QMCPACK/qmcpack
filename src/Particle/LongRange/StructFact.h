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


#ifndef QMCPLUSPLUS_STRUCTFACT_H
#define QMCPLUSPLUS_STRUCTFACT_H

//#define USE_REAL_STRUCT_FACTOR
#include "ParticleSet.h"
#include "DynamicCoordinates.h"
#include "LongRange/KContainer.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{
/** @ingroup longrange
 *\brief Calculates the structure-factor for a particle set
 *
 * Structure factor per species
 *   Rhok[alpha][k] \f$ \equiv \rho_{k}^{\alpha} = \sum_{i} e^{i{\bf k}\cdot{\bf r_i}}\f$
 * Structure factor per particle
 *   eikr[i][k]
 */
class StructFact : public QMCTraits
{
public:
  //Typedef for the lattice-type
  using ParticleLayout = PtclOnLatticeTraits::ParticleLayout;
  /** false, if the structure factor is not actively updated
   *
   * Default is false. Particle-by-particle update functions, makeMove,
   * acceptMove and rejectMove are costly and do not need to be performed
   * unless Hamiltonian uses pbyp.
   */
  bool DoUpdate;
  /** enumeration for the methods to handle mixed bconds
   *
   * Allow overwriting lattice::SuperCellEnum to use D-dim k-point sets with mixed BC
   */
  int SuperCellEnum;
  ///1-D container for the phase
  Vector<RealType> phiV;
  ///2-D container for the phase
#if defined(USE_REAL_STRUCT_FACTOR)
  Matrix<RealType> rhok_r, rhok_i;
  Matrix<RealType> eikr_r, eikr_i;
  Vector<RealType> eikr_r_temp, eikr_i_temp;
#else
  Matrix<ComplexType> rhok;
  ///eikr[particle-index][K]
  Matrix<ComplexType> eikr;
  ///eikr[K] for a proposed move
  Vector<ComplexType> eikr_temp;
#endif
  /** Constructor - copy ParticleSet and init. k-shells
   * @param nptcls number of particles
   * @param ns number of species
   * @param lattice long range box
   * @param kc cutoff for k
   */
  StructFact(int nptcls, int ns, const ParticleLayout& lattice, RealType kc);
  /// desructor
  ~StructFact();

  /** Recompute Rhok if lattice changed
   * @param kc cut-off K
   */
  void updateNewCell(const ParticleLayout& lattice, RealType kc);
  /**  Update Rhok if all particles moved
   */
  void updateAllPart(const ParticleSet& P);

  static void mw_updateAllPart(const RefVectorWithLeader<StructFact>& sk_list,
                               const RefVectorWithLeader<ParticleSet>& p_list);

  /** evaluate eikr_temp for eikr for the proposed move
   * @param active index of the moved particle
   * @param pos proposed position
   */
  void makeMove(int active, const PosType& pos);
  /** update eikr and rhok with eikr_temp
   * @param active index of the moved particle
   * @param gid group id of the active particle
   */
  void acceptMove(int active, int gid, const PosType& rold);
  /** discard any temporary data
   * @param active index of the moved particle
   * @param gid group id of the active particle
   *
   * Do nothing
   */
  void rejectMove(int active, int gid);
  /** @brief switch on the storage per particle
   * if StorePerParticle was false, this function allocates memory and precompute data
   * if StorePerParticle was true, this function is no-op
   */
  void turnOnStorePerParticle(const ParticleSet& P);

  /// accessor of StorePerParticle
  bool isStorePerParticle() const { return StorePerParticle; }

  /// access k_lists_ read only
  const KContainer& getKLists() const { return *k_lists_; }

private:
  /// Compute all rhok elements from the start
  void computeRhok(const ParticleSet& P);
  /** resize the internal data
   * @param np number of species
   * @param nptcl number of particles
   * @param nkpts
   */
  void resize(int nkpts);

  /// K-Vector List.
  const std::shared_ptr<KContainer> k_lists_;
  /// number of particles
  const size_t num_ptcls;
  /// number of species
  const size_t num_species;
  /** Whether intermediate data is stored per particle. default false
   * storing data per particle needs significant amount of memory but some calculation may request it.
   * storing data per particle specie is more cost-effective
   */
  bool StorePerParticle;
  /// timer for updateAllPart
  NewTimer& update_all_timer_;
};
} // namespace qmcplusplus

#endif
