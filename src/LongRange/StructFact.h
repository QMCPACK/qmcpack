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
#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "LongRange/KContainer.h"
#include "OhmmsPETE/OhmmsVector.h"

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
class StructFact: public QMCTraits
{

public:
  typedef PooledData<RealType>  BufferType;

  /** false, if the structure factor is not actively updated
   *
   * Default is false. Particle-by-particle update functions, makeMove,
   * acceptMove and rejectMove are costly and do not need to be performed
   * unless Hamiltonian uses pbyp.
   */
  bool DoUpdate;
  /// default false, the per particle data is not saved
  bool StorePerParticle;
  /** enumeration for the methods to handle mixed bconds
   *
   * Allow overwriting lattice::SuperCellEnum to use D-dim k-point sets with mixed BC
   */
  int SuperCellEnum;
  /// K-Vector List.
  KContainer KLists;
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
   * @param P Reference particle set
   * @param kc cutoff for k
   */
  StructFact(ParticleSet& P, RealType kc);
  /// desructor 
  ~StructFact();

  /** Recompute Rhok if lattice changed
   * @param kc cut-off K
   */
  void UpdateNewCell(ParticleSet& P, RealType kc);
  /**  Update Rhok if all particles moved
   */
  void UpdateAllPart(ParticleSet& P);

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
  /// Update Rhok if 1 particle moved
  //void Update1Part(const PosType& rold, const PosType& rnew,int iat,int GroupID);

  //Buffer methods. For PbyP MC where data must be stored in an anonymous
  //buffer between iterations
  /** @brief register rhok data to buf so that it can copyToBuffer and copyFromBuffer
   *
   * This function is used for particle-by-particle MC methods to register structure factor
   * to an anonymous buffer.
   */
  inline void registerData(BufferType& buf)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    buf.add(rhok_r.first_address(),rhok_r.last_address());
    buf.add(rhok_i.first_address(),rhok_i.last_address());
    buf.add(eikr_r.first_address(),eikr_r.last_address());
    buf.add(eikr_i.first_address(),eikr_i.last_address());
#else
    buf.add(rhok.first_address(),rhok.last_address());
    buf.add(eikr.first_address(),eikr.last_address());
#endif
  }

  /** @brief put rhok data to buf
   *
   * This function is used for particle-by-particle MC methods
   */
  inline void updateBuffer(BufferType& buf)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    buf.put(rhok_r.first_address(),rhok_r.last_address());
    buf.put(rhok_i.first_address(),rhok_i.last_address());
    buf.put(eikr_r.first_address(),eikr_r.last_address());
    buf.put(eikr_i.first_address(),eikr_i.last_address());
#else
    buf.put(rhok.first_address(),rhok.last_address());
    buf.put(eikr.first_address(),eikr.last_address());
#endif
  }
  /** @brief copy the data to an anonymous buffer
   *
   * Any data that will be used by the next iteration will be copied to a buffer.
   */
  inline void copyToBuffer(BufferType& buf)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    buf.put(rhok_r.first_address(),rhok_r.last_address());
    buf.put(rhok_i.first_address(),rhok_i.last_address());
    buf.put(eikr_r.first_address(),eikr_r.last_address());
    buf.put(eikr_i.first_address(),eikr_i.last_address());
#else
    buf.put(rhok.first_address(),rhok.last_address());
    buf.put(eikr.first_address(),eikr.last_address());
#endif
  }
  /** @brief copy the data from an anonymous buffer
   *
   * Any data that was used by the previous iteration will be copied from a buffer.
   */
  inline void copyFromBuffer(BufferType& buf)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    buf.get(rhok_r.first_address(),rhok_r.last_address());
    buf.get(rhok_i.first_address(),rhok_i.last_address());
    buf.get(eikr_r.first_address(),eikr_r.last_address());
    buf.get(eikr_i.first_address(),eikr_i.last_address());
#else
    buf.get(rhok.first_address(),rhok.last_address());
    buf.get(eikr.first_address(),eikr.last_address());
#endif
  }

  /** @brief switch on the storage per particle
   *
   * allocate the memory and precompute the data
   */
  void turnOnStorePerParticle(ParticleSet& P);

private:
  ///data for recursive evaluation for a given position
  Matrix<ComplexType> C;
  ///Compute all rhok elements from the start
  void FillRhok(ParticleSet& P);
  ///Smart update of rhok for 1-particle move. Simply supply old+new position
  void UpdateRhok(const PosType& rold,
                  const PosType& rnew,int iat,int GroupID);
  /** resize the internal data
   * @param np number of species
   * @param nptcl number of particles
   * @param nkpts
   */
  void resize(int ns, int nptcl, int nkpts);
};
}

#endif
