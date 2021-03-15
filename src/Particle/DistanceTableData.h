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
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "OhmmsSoA/VectorSoaContainer.h"
#include <limits>
#include <bitset>

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
  const ParticleSet* Origin;

  const int N_sources;
  const int N_targets;

  /**defgroup SoA data */
  /*@{*/
  /** distances_[i][j] , [N_targets][N_sources]
   *  Note: Derived classes decide if it is a memory view or the actual storage
   *        For derived AA, only the lower triangle (j<i) is defined and up-to-date after pbyp move.
   *          The upper triangle is symmetric to the lower one only when the full table is evaluated from scratch.
   *          Avoid using the upper triangle because we may change the code to only allocate the lower triangle part.
   *        For derived AB, the full table is up-to-date after pbyp move
   */
  std::vector<DistRow> distances_;

  /** displacements_[N_targets]x[3][N_sources]
   *  Note: Derived classes decide if it is a memory view or the actual storage
   *        displacements_[i][j] = r_A2[j] - r_A1[i], the opposite sign of AoS dr
   *        For derived AA, A1=A2=A, only the lower triangle (j<i) is defined.
   *        For derived AB, A1=A, A2=B, the full table is allocated.
   */
  std::vector<DisplRow> displacements_;

  /** temp_r */
  DistRow temp_r_;

  /** temp_dr */
  DisplRow temp_dr_;
  /*@}*/

  /** whether full table needs to be ready at anytime or not
   * Optimization can be implemented during forward PbyP move when the full table is not needed all the time.
   * DT consumers should know if full table is needed or not and request via addTable.
   */
  bool need_full_table_;

  /** set to particle id after move() with prepare_old = true. -1 means not prepared.
   * It is intended only for safety checks, not for codepath selection.
   */
  int old_prepared_elec_id;

  ///name of the table
  const std::string name_;

public:
  ///constructor using source and target ParticleSet
  DistanceTableData(const ParticleSet& source, const ParticleSet& target)
      : Origin(&source),
        N_sources(source.getTotalNum()),
        N_targets(target.getTotalNum()),
        need_full_table_(false),
        old_prepared_elec_id(-1),
        name_(source.getName() + "_" + target.getName())
  {}

  ///virutal destructor
  virtual ~DistanceTableData() = default;

  ///get need_full_table_
  inline bool getFullTableNeeds() const { return need_full_table_; }

  ///set need_full_table_
  inline void setFullTableNeeds(bool is_needed) { need_full_table_ = is_needed; }

  ///return the name of table
  inline const std::string& getName() const { return name_; }

  ///returns the reference the origin particleset
  const ParticleSet& origin() const { return *Origin; }

  ///returns the number of centers
  inline IndexType centers() const { return Origin->getTotalNum(); }

  ///returns the number of centers
  inline IndexType targets() const { return N_targets; }

  ///returns the number of source particles
  inline IndexType sources() const { return N_sources; }

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
    APP_ABORT("DistanceTableData::getOldDists is used incorrectly! Contact developers on github.");
    return temp_r_; // dummy return to avoid compiler warning.
  }

  /** return old displacements set up by move() for optimized distance table consumers
   */
  virtual const DisplRow& getOldDispls() const
  {
    APP_ABORT("DistanceTableData::getOldDispls is used incorrectly! Contact developers on github.");
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

  /** evaluate the temporary pair relations when a move is proposed
   * @param P the target particle set
   * @param rnew proposed new position
   * @param iat the particle to be moved
   * @param prepare_old if true, prepare (temporary) old distances and displacements for using getOldDists and getOldDispls functions in acceptMove.
   *
   * Note: some distance table consumers (WaveFunctionComponent) have optimized code paths which require prepare_old = true for accepting a move.
   * Drivers/Hamiltonians know whether moves will be accepted or not and manage this flag when calling ParticleSet::makeMoveXXX functions.
   */
  virtual void move(const ParticleSet& P, const PosType& rnew, const IndexType iat = 0, bool prepare_old = true) = 0;
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
  virtual int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
    APP_ABORT("DistanceTableData::get_first_neighbor is not implemented in calling base class");
    return 0;
  }

  inline void print(std::ostream& os)
  {
    APP_ABORT("DistanceTableData::print is not supported")
    //os << "Table " << Origin->getName() << std::endl;
    //for (int i = 0; i < r_m.size(); i++)
    //  os << r_m[i] << " ";
    //os << std::endl;
  }

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
