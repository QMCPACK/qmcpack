//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H

#include "Particle/ParticleSet.h"
#include <limits>
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "DTModes.h"

namespace qmcplusplus
{
class ResourceCollection;

/** @ingroup nnlist
 * @brief Abstract class to manage pair data between two ParticleSets.
 *
 * Each DistanceTableData object is fined by Source and Target of ParticleSet types.
 *
 */
class DistanceTableData
{
public:
  static constexpr unsigned DIM = OHMMS_DIM;

  using IndexType = QMCTraits::IndexType;
  using RealType  = QMCTraits::RealType;
  using PosType   = QMCTraits::PosType;
  using DistRow   = Vector<RealType, aligned_allocator<RealType>>;
  using DisplRow  = VectorSoaContainer<RealType, DIM>;

protected:
  // FIXME. once DT takes only DynamicCoordinates, change this type as well.
  const ParticleSet& origin_;

  const size_t num_sources_;
  const size_t num_targets_;

  ///name of the table
  const std::string name_;

  /**defgroup SoA data */
  /*@{*/
  /** distances_[i][j] , [num_targets_][num_sources_]
   *  Note: Derived classes decide if it is a memory view or the actual storage
   *        For derived AA, only the lower triangle (j<i) data can be accessed safely.
   *            There is no bound check to protect j>=i terms as the nature of operator[].
   *            When the storage of the table is allocated as a single memory segment,
   *            out-of-bound access is still within the segment and
   *            thus doesn't trigger an alarm by the address sanitizer.
   *        For derived AB, the full table is up-to-date after pbyp move
   */
  std::vector<DistRow> distances_;

  /** displacements_[num_targets_]x[3][num_sources_]
   *  Note: Derived classes decide if it is a memory view or the actual storage
   *        displacements_[i][j] = r_A2[j] - r_A1[i], the opposite sign of AoS dr
   *        For derived AA, A1=A2=A, only the lower triangle (j<i) is defined. See the note of distances_
   *        For derived AB, A1=A, A2=B, the full table is allocated.
   */
  std::vector<DisplRow> displacements_;

  /** temp_r */
  DistRow temp_r_;

  /** temp_dr */
  DisplRow temp_dr_;
  /*@}*/

  ///operation modes defined by DTModes
  DTModes modes_;

public:
  ///constructor using source and target ParticleSet
  DistanceTableData(const ParticleSet& source, const ParticleSet& target, DTModes modes)
      : origin_(source),
        num_sources_(source.getTotalNum()),
        num_targets_(target.getTotalNum()),
        name_(source.getName() + "_" + target.getName()),
        modes_(modes)
  {}

  /// copy constructor. deleted
  DistanceTableData(const DistanceTableData&) = delete;

  ///virutal destructor
  virtual ~DistanceTableData() = default;

  ///get modes
  inline DTModes getModes() const { return modes_; }

  ///set modes
  inline void setModes(DTModes modes) { modes_ = modes; }

  ///return the name of table
  inline const std::string& getName() const { return name_; }

  ///returns the reference the origin particleset
  const ParticleSet& get_origin() const { return origin_; }

  ///returns the number of centers
  inline size_t centers() const { return origin_.getTotalNum(); }

  ///returns the number of centers
  inline size_t targets() const { return num_targets_; }

  ///returns the number of source particles
  inline size_t sources() const { return num_sources_; }

  /// return multi walker temporary pair distance table data pointer
  virtual const RealType* getMultiWalkerTempDataPtr() const
  {
    throw std::runtime_error(name_ + " multi walker data pointer for temp not supported");
    return nullptr;
  }

  /// return multi-walker full (all pairs) distance table data pointer
  virtual const RealType* getMultiWalkerDataPtr() const
  {
    throw std::runtime_error(name_ + " multi walker data pointer not supported");
    return nullptr;
  }

  /// return stride of per target pctl data. full table data = stride * num of target particles
  virtual size_t getPerTargetPctlStrideSize() const
  {
    throw std::runtime_error(name_ + " getPerTargetPctlStrideSize not supported");
    return 0;
  }

  /** return full table distances
   */
  const std::vector<DistRow>& getDistances() const { return distances_; }

  /** return full table displacements
   */
  const std::vector<DisplRow>& getDisplacements() const { return displacements_; }

  /** return a row of distances for a given target particle
   */
  const DistRow& getDistRow(int iel) const { return distances_[iel]; }

  /** return a row of displacements for a given target particle
   */
  const DisplRow& getDisplRow(int iel) const { return displacements_[iel]; }

  /** return old distances set up by move() for optimized distance table consumers
   */
  virtual const DistRow& getOldDists() const
  {
    throw std::runtime_error("DistanceTableData::getOldDists is used incorrectly! Contact developers on github.");
    return temp_r_; // dummy return to avoid compiler warning.
  }

  /** return old displacements set up by move() for optimized distance table consumers
   */
  virtual const DisplRow& getOldDispls() const
  {
    throw std::runtime_error("DistanceTableData::getOldDispls is used incorrectly! Contact developers on github.");
    return temp_dr_; // dummy return to avoid compiler warning.
  }

  /** return the temporary distances when a move is proposed
   */
  const DistRow& getTempDists() const { return temp_r_; }

  /** return the temporary displacements when a move is proposed
   */
  const DisplRow& getTempDispls() const { return temp_dr_; }

  /** evaluate the full Distance Table
   * @param P the target particle set
   */
  virtual void evaluate(ParticleSet& P) = 0;
  virtual void mw_evaluate(const RefVectorWithLeader<DistanceTableData>& dt_list,
                           const RefVectorWithLeader<ParticleSet>& p_list) const
  {
#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      dt_list[iw].evaluate(p_list[iw]);
  }

  /** recompute multi walker internal data, recompute
   * @param dt_list the distance table batch
   * @param p_list the target particle set batch
   * @param recompute if true, must recompute. Otherwise, implementation dependent.
   */
  virtual void mw_recompute(const RefVectorWithLeader<DistanceTableData>& dt_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            const std::vector<bool>& recompute) const
  {
#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      if (recompute[iw])
        dt_list[iw].evaluate(p_list[iw]);
  }

  /** evaluate the temporary pair relations when a move is proposed
   * @param P the target particle set
   * @param rnew proposed new position
   * @param iat the particle to be moved
   * @param prepare_old if true, prepare (temporary) old distances and displacements for using getOldDists and getOldDispls functions in acceptMove.
   *
   * Note: some distance table consumers (WaveFunctionComponent) have optimized code paths which require prepare_old = true for accepting a move.
   * Drivers/Hamiltonians know whether moves will be accepted or not and manage this flag when calling ParticleSet::makeMoveXXX functions.
   */
  virtual void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old = true) = 0;

  /** walker batched version of move. this function may be implemented asynchronously.
   * Additional synchroniziation for collecting results should be handled by the caller.
   * If DTModes::NEED_TEMP_DATA_ON_HOST, host data will be updated.
   * If no consumer requests data on the host, the transfer is skipped.
   */
  virtual void mw_move(const RefVectorWithLeader<DistanceTableData>& dt_list,
                       const RefVectorWithLeader<ParticleSet>& p_list,
                       const std::vector<PosType>& rnew_list,
                       const IndexType iat = 0,
                       bool prepare_old    = true) const
  {
#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      dt_list[iw].move(p_list[iw], rnew_list[iw], iat, prepare_old);
  }

  /** update the distance table by the pair relations from the temporal position.
   *  Used when a move is accepted in regular mode
   * @param iat the particle with an accepted move
   */
  virtual void update(IndexType jat) = 0;

  /** fill partially the distance table by the pair relations from the temporary or old particle position.
   *  Used in forward mode when a move is reject
   * @param iat the particle with an accepted move
   * @param from_temp if true, copy from temp. if false, copy from old
   */
  virtual void updatePartial(IndexType jat, bool from_temp)
  {
    if (from_temp)
      update(jat);
  }

  /** walker batched version of updatePartial.
   * If not DTModes::NEED_TEMP_DATA_ON_HOST, host data is not up-to-date and host distance table will not be updated.
   */
  virtual void mw_updatePartial(const RefVectorWithLeader<DistanceTableData>& dt_list,
                                IndexType jat,
                                const std::vector<bool>& from_temp)
  {
#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      dt_list[iw].updatePartial(jat, from_temp[iw]);
  }

  /** finalize distance table calculation after particle-by-particle moves
   * if update() doesn't make the table up-to-date during p-by-p moves
   * finalizePbyP takes action to bring the table up-to-date
   */
  virtual void finalizePbyP(const ParticleSet& P) {}

  /** walker batched version of finalizePbyP
   * If not DTModes::NEED_TEMP_DATA_ON_HOST, host distance table data is not updated at all during p-by-p
   * Thus, a recompute is necessary to update the whole host distance table for consumers like the Coulomb potential.
   */
  virtual void mw_finalizePbyP(const RefVectorWithLeader<DistanceTableData>& dt_list,
                               const RefVectorWithLeader<ParticleSet>& p_list) const
  {
#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      dt_list[iw].finalizePbyP(p_list[iw]);
  }

  /** build a compact list of a neighbor for the iat source
   * @param iat source particle id
   * @param rcut cutoff radius
   * @param jid compressed index
   * @param dist compressed distance
   * @param displ compressed displacement
   * @return number of target particles within rcut
   */
  virtual size_t get_neighbors(int iat,
                               RealType rcut,
                               int* restrict jid,
                               RealType* restrict dist,
                               PosType* restrict displ) const
  {
    return 0;
  }

  /** find the first nearest neighbor
   * @param iat source particle id
   * @param r distance
   * @param dr displacement
   * @param newpos if true, use the data in temp_r_ and temp_dr_ for the proposed move.
   *        if false, use the data in distance_[iat] and displacements_[iat]
   * @return the id of the nearest particle, -1 not found
   */
  virtual int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const = 0;

  inline void print(std::ostream& os) { throw std::runtime_error("DistanceTableData::print is not supported"); }

  /// initialize a shared resource and hand it to a collection
  virtual void createResource(ResourceCollection& collection) const {}

  /// acquire a shared resource from a collection
  virtual void acquireResource(ResourceCollection& collection,
                               const RefVectorWithLeader<DistanceTableData>& dt_list) const
  {}

  /// return a shared resource to a collection
  virtual void releaseResource(ResourceCollection& collection,
                               const RefVectorWithLeader<DistanceTableData>& dt_list) const
  {}
};
} // namespace qmcplusplus
#endif
