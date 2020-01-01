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
#include "simd/allocator.hpp"
#include <OhmmsSoA/VectorSoaContainer.h>
#include <limits>
#include <bitset>

namespace qmcplusplus
{
#ifndef ENABLE_SOA
/** @defgroup nnlist Distance-table group
 * @brief class to manage a set of data for distance relations between ParticleSet objects.
 */
template<class T, unsigned N>
struct TempDisplacement
{
  ///new distance
  T r1;
  ///inverse of the new distance
  T rinv1;
  ///new displacement
  TinyVector<T, N> dr1;
  inline TempDisplacement() : r1(0.0), rinv1(0.0) {}
  inline void reset()
  {
    r1    = 0.0;
    rinv1 = 0.0;
    dr1   = 0.0;
  }
};
#endif

/** enumerator for DistanceTableData::DTType
 *
 * - DT_AOS Use original AoS type
 * - DT_SOA Use SoA type
 * - DT_AOS_PREFERRED Create AoS type, if possible.
 * - DT_SOA_PREFERRED Create SoA type, if possible.
 * The first user of each pair will decide the type of distance table.
 * It is the responsibility of the user class to check DTType.
 */
enum DistTableType
{
  DT_AOS = 0,
  DT_SOA,
  DT_AOS_PREFERRED,
  DT_SOA_PREFERRED
};

/** @ingroup nnlist
 * @brief Abstract class to manage pair data between two ParticleSets.
 *
 * Each DistanceTableData object is fined by Source and Target of ParticleSet types.
 *
 */
struct DistanceTableData
{
  static constexpr unsigned DIM = OHMMS_DIM;

  using IndexType    = QMCTraits::IndexType;
  using RealType     = QMCTraits::RealType;
  using PosType      = QMCTraits::PosType;
  using DistRowType  = Vector<RealType, aligned_allocator<RealType>>;
  using DisplRowType = VectorSoaContainer<RealType, DIM>;
#ifndef ENABLE_SOA
  using IndexVectorType = aligned_vector<IndexType>;
  using TempDistType    = TempDisplacement<RealType, DIM>;
  using ripair          = std::pair<RealType, IndexType>;
#endif

  ///Type of DT
  int DTType;

  const ParticleSet* Origin;

  int N_sources;
  int N_targets;
  int N_walkers;

#ifndef ENABLE_SOA
  ///number of pairs
  int npairs_m;

  /** @brief M.size() = N_sources+1
   *
   * M[i+i] - M[i] = the number of connected points to the i-th source
   */
  IndexVectorType M;

  /** @brief J.size() = M[N_sources]]
   *
   * J[nn] = the index of the connected point for the i-th point
   * satisfying  \f$M[i] <= nn < M[i+i]\f$
   */
  IndexVectorType J;

  /** @brief PairID.size() = M[N[SourceIndex]]
   *
   * PairID[nn] = the index of the connected point for the i-th point
   * satisfying  \f$PairIDM[i] <= nn < PairID[i+i]\f$
   */
  IndexVectorType PairID;

  /** Locator of the pair  */
  IndexVectorType IJ;

  /** @brief A NN relation of all the source particles with respect to an activePtcl
   *
   * This data is for particle-by-particle move.
   * When a MC move is propsed to the activePtcl, the old and new distance relation
   * is stored in Temp. When the move is accepted, the new data replace the old.
   * If the move is rejected, nothing is done and new data will be overwritten.
   */
  std::vector<TempDistType> Temp;
#endif

protected:
  /**defgroup SoA data */
  /*@{*/
  /** Distances[i][j] , [N_targets][N_sources]
   *  Note: Derived classes decide if it is a memory view or the actaully storage
   *        For derived AA, only the lower triangle (j<i) is up-to-date after pbyp move
   *          The upper triangle is symmetric to the lower one only when the full table is evaluated from scratch.
   *          Avoid using the upper triangle because we may change the code to only allocate the lower triangle part.
   *        For derived BA, the full table is up-to-date after pbyp move
   */
  std::vector<DistRowType> Distances;

  /** Displacements[N_targets]x[3][N_sources]
   *  Note: Derived classes decide if it is a memory view or the actaully storage
   *        Displacements[i][j] = r_A2[j] - r_A1[i], the opposite sign of AoS dr
   *        For derived AA, A1=A2=A, only the lower triangle (j<i) is allocated in memoryPool_displs_
   *          For this reason, Displacements[i] and Displacements[i+1] overlap in memory
   *          and they must be updated in order during PbyP move.
   *        For derived BA, A1=A, A2=B, the full table is allocated.
   */
  std::vector<DisplRowType> Displacements;

  /** temp_r */
  DistRowType Temp_r;

  /** temp_dr */
  DisplRowType Temp_dr;
  /*@}*/

  /** whether full table needs to be ready at anytime or not
   * Optimization can be implemented during forward PbyP move when the full table is not needed all the time.
   * DT consumers should know if full table is needed or not and request via addTable.
   */
  bool need_full_table_;

  ///name of the table
  std::string Name;

public:
  ///constructor using source and target ParticleSet
  DistanceTableData(const ParticleSet& source, const ParticleSet& target)
      : Origin(&source), N_sources(0), N_targets(0), N_walkers(0), need_full_table_(false)
  {}

  ///virutal destructor
  virtual ~DistanceTableData() {}

  ///get need_full_table_
  inline bool getFullTableNeeds() const { return need_full_table_; }

  ///set need_full_table_
  inline void setFullTableNeeds(bool is_needed) { need_full_table_ = is_needed; }

  ///return the name of table
  inline const std::string& getName() const { return Name; }

  ///set the name of table
  inline void setName(const std::string& tname) { Name = tname; }

  ///returns the reference the origin particleset
  const ParticleSet& origin() const { return *Origin; }

  inline bool is_same_type(int dt_type) const { return DTType == dt_type; }

#ifndef ENABLE_SOA
  //@{access functions to the distance, inverse of the distance and directional consine vector
  inline PosType dr(int j) const { return dr_m[j]; }
  inline RealType r(int j) const { return r_m[j]; }
  inline RealType rinv(int j) const { return rinv_m[j]; }
  //@}
#endif

  ///returns the number of centers
  inline IndexType centers() const { return Origin->getTotalNum(); }

  ///returns the number of centers
  inline IndexType targets() const { return N_targets; }

  ///returns the number of source particles
  inline IndexType sources() const { return N_sources; }

#ifndef ENABLE_SOA
  inline IndexType getTotNadj() const { return npairs_m; }

  /// return the distance |R[iadj(i,nj)]-R[i]|
  inline RealType distance(int i, int nj) const { return r_m[M[i] + nj]; }

  /// return the displacement R[iadj(i,nj)]-R[i]
  inline PosType displacement(int i, int nj) const { return dr_m[M[i] + nj]; }

  //!< Returns a number of neighbors of the i-th ptcl.
  inline IndexType nadj(int i) const { return M[i + 1] - M[i]; }
  //!< Returns the id of nj-th neighbor for i-th ptcl
  inline IndexType iadj(int i, int nj) const { return J[M[i] + nj]; }
  //!< Returns the id of j-th neighbor for i-th ptcl
  inline IndexType loc(int i, int j) const { return M[i] + j; }

  /** search the closest source particle within rcut
   * @param rcut cutoff radius
   * @return the index of the closest particle
   *
   * Return -1 if none is within rcut
   * @note It searches the temporary list after a particle move is made
   */
  inline IndexType find_closest_source(RealType rcut) const
  {
    int i = 0;
    while (i < N_sources)
    {
      if (Temp[i].r1 < rcut)
        return i;
      i++;
    }
    return -1;
  }

  /** search the closest source particle of the iel-th particle within rcut 
   * @param rcut cutoff radius
   * @return the index of the first particle
   *
   * Return -1 if none is within rcut
   * @note Check the real distance table and only works for the AsymmetricDistanceTable
   */
  inline IndexType find_closest_source(int iel, RealType rcut) const
  {
    int i = 0, nn = iel;
    while (nn < r_m.size())
    {
      if (r_m[nn] < rcut)
        return i;
      nn += N_targets;
      i++;
    }
    return -1;
  }
#endif

  /** return full table distances
   */
  const std::vector<DistRowType>& getDistances() const { return Distances; }

  /** return full table displacements
   */
  const std::vector<DisplRowType>& getDisplacements() const { return Displacements; }

  /** return a row of distances for a given target particle
   */
  const DistRowType& getDistRow(int iel) const { return Distances[iel]; }

  /** return a row of displacements for a given target particle
   */
  const DisplRowType& getDisplRow(int iel) const { return Displacements[iel]; }

  /** return old distances set up by move() for optimized distance table consumers
   */
  virtual const DistRowType& getOldDists() const
  {
    APP_ABORT("DistanceTableData::getOldDists is used incorrectly! Contact developers on github.");
    return Temp_r; // dummy return to avoid compiler warning.
  }

  /** return old displacements set up by move() for optimized distance table consumers
   */
  virtual const DisplRowType& getOldDispls() const
  {
    APP_ABORT("DistanceTableData::getOldDispls is used incorrectly! Contact developers on github.");
    return Temp_dr; // dummy return to avoid compiler warning.
  }

  /** return the temporary distances when a move is proposed
   */
  const DistRowType& getTemporalDists() const { return Temp_r; }

  /** return the temporary displacements when a move is proposed
   */
  const DisplRowType& getTemporalDispls() const { return Temp_dr; }

  /** evaluate the full Distance Table
   * @param P the target particle set
   */
  virtual void evaluate(ParticleSet& P) = 0;

  /** evaluate the temporary pair relations when a move is proposed
   * @param P the target particle set
   * @param rnew proposed new position
   * @param iat the particle to be moved
   * @param prepare_old if true, prepare old distances and displacements for using getOldDists and getOldDispls functions later.
   *
   * Note: some distance table consumers (WaveFunctionComponent) have optimized code paths which require prepare_old = true for accepting a move.
   * Drivers/Hamiltonians know whether moves will be accepted or not and manager this flag when calling ParticleSet::makeMoveXXX functions.
   */
  virtual void move(const ParticleSet& P, const PosType& rnew, const IndexType iat = 0, bool prepare_old = true) = 0;

  /** update the distance table by the pair relations if a move is accepted
   * @param iat the particle with an accepted move
   * @param forward If true, rows after iat will not be updated. If false, upon accept a move, the full table should be up-to-date
   */
  virtual void update(IndexType jat, bool forward = false) = 0;

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
   * @param newpos if true, use the data in Temp_r and Temp_dr for the proposed move.
   *        if false, use the data in Distance[iat] and Displacements[iat]
   * @return the id of the nearest particle, -1 not found
   */
  virtual int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
    APP_ABORT("DistanceTableData::get_first_neighbor is not implemented in calling base class");
    return 0;
  }

#ifndef ENABLE_SOA
  /** build a compact list of a neighbor for the iat source
   * @param iat source particle id
   * @param rcut cutoff radius
   * @param dist compressed distance
   * @return number of target particles within rcut
   */
  virtual size_t get_neighbors(int iat, RealType rcut, RealType* restrict dist) const { return 0; }

  /// find index and distance of each nearest neighbor particle
  virtual void nearest_neighbor(std::vector<ripair>& ri, bool transposed = false) const
  {
    APP_ABORT("DistanceTableData::nearest_neighbor is not implemented in calling base class");
  }

  /// find indices and distances of nearest neighbors particles to particle n
  virtual void nearest_neighbors(int n, int neighbors, std::vector<ripair>& ri, bool transposed = false)
  {
    APP_ABORT("DistanceTableData::nearest_neighbors is not implemented in calling base class");
  }

  /// find species resolved indices and distances of nearest particles to particle n
  virtual void nearest_neighbors_by_spec(int n,
                                         int neighbors,
                                         int spec_start,
                                         std::vector<ripair>& ri,
                                         bool transposed = false)
  {
    APP_ABORT("DistanceTableData::nearest_neighbors is not implemented in calling base class");
  }

  inline void check_neighbor_size(std::vector<ripair>& ri, bool transposed = false) const
  {
    int m;
    if (transposed)
      m = N_sources;
    else
      m = N_targets;
    if (ri.size() != m)
      APP_ABORT("DistanceTableData::check_neighbor_size  distance/index vector length is not equal to the number of "
                "neighbor particles");
  }
#endif

  inline void print(std::ostream& os)
  {
    os << "Table " << Origin->getName() << std::endl;
#ifndef ENABLE_SOA
    for (int i = 0; i < r_m.size(); i++)
      os << r_m[i] << " ";
    os << std::endl;
#endif
  }

#ifndef ENABLE_SOA
  /**defgroup storage data for nearest-neighbor relations
   */
  /*@{*/
  /** Cartesian distance \f$r(i,j) = |R(j)-R(i)|\f$ */
  std::vector<RealType> r_m;
  /** Cartesian distance \f$rinv(i,j) = 1/r(i,j)\f$ */
  std::vector<RealType> rinv_m;
  /** displacement vectors \f$dr(i,j) = R(j)-R(i)\f$  */
  std::vector<PosType> dr_m;
  /*@}*/
#endif

  /**resize the storage
   *@param npairs number of pairs which is evaluated by a derived class
   *@param nw number of copies
   *
   * The data for the pair distances, displacements
   *and the distance inverses are stored in a linear storage.
   * The logical view of these storages is (ipair,iwalker),
   * where 0 <= ipair < M[N[SourceIndex]] and 0 <= iwalker < N[WalkerIndex]
   * This scheme can handle both dense and sparse distance tables,
   * and full or half of the pairs.
   * Note that this function is protected and the derived classes are
   * responsible to call this function for memory allocation and any
   * change in the indices N.
   */
  void resize(int npairs, int nw)
  {
    N_walkers = nw;
    //if(nw==1)
    {
#ifndef ENABLE_SOA
      dr_m.resize(npairs);
      r_m.resize(npairs);
      rinv_m.resize(npairs);
      Temp.resize(N_sources);
#endif
    }
  }
};
} // namespace qmcplusplus
#endif
