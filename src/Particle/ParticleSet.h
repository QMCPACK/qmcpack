//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PARTICLESET_H
#define QMCPLUSPLUS_PARTICLESET_H

#include <memory>
#include <Configuration.h>
#include "ParticleTags.h"
#include "DynamicCoordinates.h"
#include "Walker.h"
#include "SpeciesSet.h"
#include "Pools/PooledData.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Utilities/TimerManager.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "type_traits/template_types.hpp"
#include "SimulationCell.h"
#include "DTModes.h"

namespace qmcplusplus
{
///forward declaration of DistanceTable
class DistanceTable;
class DistanceTableAA;
class DistanceTableAB;
class ResourceCollection;
class StructFact;

/** Specialized paritlce class for atomistic simulations
 *
 * Derived from QMCTraits, ParticleBase<PtclOnLatticeTraits> and OhmmsElementBase.
 * The ParticleLayout class represents a supercell with/without periodic boundary
 * conditions. The ParticleLayout class also takes care of spatial decompositions
 * for efficient evaluations for the interactions with a finite cutoff.
 */
class ParticleSet : public QMCTraits, public OhmmsElementBase, public PtclOnLatticeTraits
{
public:
  /// walker type
  using Walker_t = Walker<QMCTraits, PtclOnLatticeTraits>;
  /// container type to store the property
  using PropertyContainer_t = Walker_t::PropertyContainer_t;
  /// buffer type for a serialized buffer
  using Buffer_t = PooledData<RealType>;

  enum quantum_domains
  {
    no_quantum_domain = 0,
    classical,
    quantum
  };

  ///quantum_domain of the particles, default = classical
  quantum_domains quantum_domain;

  //@{ public data members
  ///Species ID
  ParticleIndex GroupID;
  ///Position
  ParticlePos R;
  ///internal spin variables for dynamical spin calculations
  ParticleScalar spins;
  ///gradients of the particles
  ParticleGradient G;
  ///laplacians of the particles
  ParticleLaplacian L;
  ///mass of each particle
  ParticleScalar Mass;
  ///charge of each particle
  ParticleScalar Z;

  ///the index of the active bead for particle-by-particle moves
  Index_t activeBead;
  ///the direction reptile traveling
  Index_t direction;

  ///Particle density in G-space for MPC interaction
  std::vector<TinyVector<int, OHMMS_DIM>> DensityReducedGvecs;
  std::vector<ComplexType> Density_G;
  Array<RealType, OHMMS_DIM> Density_r;

  /// DFT potential
  std::vector<TinyVector<int, OHMMS_DIM>> VHXCReducedGvecs;
  std::vector<ComplexType> VHXC_G[2];
  Array<RealType, OHMMS_DIM> VHXC_r[2];

  /** name-value map of Walker Properties
   *
   * PropertyMap is used to keep the name-value mapping of
   * Walker_t::Properties.  PropertyList::Values are not
   * necessarily updated during the simulations.
   */
  PropertySetType PropertyList;

  /** properties of the current walker
   *
   * The internal order is identical to PropertyList, which holds
   * the matching names.
   */
  PropertyContainer_t Properties;

  /** observables in addition to those registered in Properties/PropertyList
   *
   * Such observables as density, gofr, sk are not stored per walker but
   * collected during QMC iterations.
   */
  Buffer_t Collectables;

  ///Property history vector
  std::vector<std::vector<FullPrecRealType>> PropertyHistory;
  std::vector<int> PHindex;
  ///@}

  ///current MC step
  int current_step;

  ///default constructor
  ParticleSet(const SimulationCell& simulation_cell, const DynamicCoordinateKind kind = DynamicCoordinateKind::DC_POS);

  ///copy constructor
  ParticleSet(const ParticleSet& p);

  ///default destructor
  ~ParticleSet() override;

  /** create  particles
   * @param n number of particles
   */
  void create(int n);
  /** create grouped particles
   * @param agroup number of particles per group
   */
  void create(const std::vector<int>& agroup);

  ///write to a std::ostream
  bool get(std::ostream&) const override;

  ///read from std::istream
  bool put(std::istream&) override;

  ///reset member data
  void reset() override;

  ///initialize ParticleSet from xmlNode
  bool put(xmlNodePtr cur) override;

  ///specify quantum_domain of particles
  void setQuantumDomain(quantum_domains qdomain);

  void set_quantum() { quantum_domain = quantum; }

  inline bool is_classical() const { return quantum_domain == classical; }

  inline bool is_quantum() const { return quantum_domain == quantum; }

  ///check whether quantum domain is valid for particles
  inline bool quantumDomainValid(quantum_domains qdomain) const { return qdomain != no_quantum_domain; }

  ///check whether quantum domain is valid for particles
  inline bool quantumDomainValid() const { return quantumDomainValid(quantum_domain); }

  /** add a distance table
   * @param psrc source particle set
   * @param modes bitmask DistanceTable::DTModes
   *
   * if this->myName == psrc.getName(), AA type. Otherwise, AB type.
   */
  int addTable(const ParticleSet& psrc, DTModes modes = DTModes::ALL_OFF);

  ///get a distance table by table_ID
  inline auto& getDistTable(int table_ID) const { return *DistTables[table_ID]; }
  ///get a distance table by table_ID and dyanmic_cast to DistanceTableAA
  const DistanceTableAA& getDistTableAA(int table_ID) const;
  ///get a distance table by table_ID and dyanmic_cast to DistanceTableAB
  const DistanceTableAB& getDistTableAB(int table_ID) const;

  /** reset all the collectable quantities during a MC iteration
   */
  inline void resetCollectables() { std::fill(Collectables.begin(), Collectables.end(), 0.0); }

  /** update the internal data
   *@param skip SK update if skipSK is true
   */
  void update(bool skipSK = false);

  /// batched version of update
  static void mw_update(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK = false);

  /** create Structure Factor with PBCs
   */
  void createSK();

  bool hasSK() const { return bool(structure_factor_); }
  /** return Structure Factor
   */
  const StructFact& getSK() const
  {
    assert(structure_factor_);
    return *structure_factor_;
  };

  /** Turn on per particle storage in Structure Factor
   */
  void turnOnPerParticleSK();

  /** Get state (on/off) of per particle storage in Structure Factor
   */
  bool getPerParticleSKState() const;

  ///retrun the SpeciesSet of this particle set
  inline SpeciesSet& getSpeciesSet() { return my_species_; }
  ///retrun the const SpeciesSet of this particle set
  inline const SpeciesSet& getSpeciesSet() const { return my_species_; }

  ///return parent's name
  inline const std::string& parentName() const { return ParentName; }
  inline void setName(const std::string& aname)
  {
    myName = aname;
    if (ParentName == "0")
    {
      ParentName = aname;
    }
  }

  inline const DynamicCoordinates& getCoordinates() const { return *coordinates_; }
  inline void setCoordinates(const ParticlePos& R) { return coordinates_->setAllParticlePos(R); }

  void resetGroups();


  const auto& getSimulationCell() const { return simulation_cell_; }
  const auto& getLattice() const { return simulation_cell_.getLattice(); }
  auto& getPrimitiveLattice() const { return const_cast<ParticleLayout&>(simulation_cell_.getPrimLattice()); }
  const auto& getLRBox() const { return simulation_cell_.getLRBox(); }

  inline bool isSameMass() const { return same_mass_; }
  inline bool isGrouped() const { return is_grouped_; }
  inline bool isSpinor() const { return is_spinor_; }
  inline void setSpinor(bool is_spinor) { is_spinor_ = is_spinor; }

  /// return active particle id
  inline Index_t getActivePtcl() const { return active_ptcl_; }
  inline const PosType& getActivePos() const { return active_pos_; }
  inline Scalar_t getActiveSpinVal() const { return active_spin_val_; }

  /// return the active position if the particle is active or the return current position if not
  inline const PosType& activeR(int iat) const
  {
    // When active_ptcl_ == iat, a move has been proposed.
    return (active_ptcl_ == iat) ? active_pos_ : R[iat];
  }

  /// return the active spin value if the particle is active or return the current spin value if not
  inline const Scalar_t& activeSpin(int iat) const
  {
    // When active_ptcl_ == iat, a move has been proposed.
    return (active_ptcl_ == iat) ? active_spin_val_ : spins[iat];
  }

  /** move the iat-th particle to active_pos_
   * @param iat the index of the particle to be moved
   * @param displ the displacement of the iat-th particle position
   * @param maybe_accept if false, the caller guarantees that the proposed move will not be accepted.
   *
   * Update active_ptcl_ index and active_pos_ position (R[iat]+displ) for a proposed move.
   * Evaluate the related distance table data DistanceTable::Temp.
   * If maybe_accept = false, certain operations for accepting moves will be skipped for optimal performance.
   */
  void makeMove(Index_t iat, const SingleParticlePos& displ, bool maybe_accept = true);
  /// makeMove, but now includes an update to the spin variable
  void makeMoveWithSpin(Index_t iat, const SingleParticlePos& displ, const Scalar_t& sdispl);

  /// batched version of makeMove
  static void mw_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                          int iat,
                          const std::vector<SingleParticlePos>& displs);

  /// batched version of makeMoveWithSpin
  static void mw_makeMoveWithSpin(const RefVectorWithLeader<ParticleSet>& p_list,
                                  int iat,
                                  const std::vector<SingleParticlePos>& displs,
                                  const std::vector<Scalar_t>& sdispls);

  /** move the iat-th particle to active_pos_
   * @param iat the index of the particle to be moved
   * @param displ random displacement of the iat-th particle
   * @return true, if the move is valid
   *
   * Update active_ptcl_ index and active_pos_ position (R[iat]+displ) for a proposed move.
   * Evaluate the related distance table data DistanceTable::Temp.
   *
   * When a Lattice is defined, passing two checks makes a move valid.
   * outOfBound(displ): invalid move, if displ is larger than half, currently, of the box in any direction
   * isValid(Lattice.toUnit(active_pos_)): invalid move, if active_pos_ goes out of the Lattice in any direction marked with open BC.
   * Note: active_pos_ and distances tables are always evaluated no matter the move is valid or not.
   */
  bool makeMoveAndCheck(Index_t iat, const SingleParticlePos& displ);
  /// makeMoveAndCheck, but now includes an update to the spin variable
  bool makeMoveAndCheckWithSpin(Index_t iat, const SingleParticlePos& displ, const Scalar_t& sdispl);

  /** Handles virtual moves for all the particles to a single newpos.
   *
   * The state active_ptcl_ remains -1 and rejectMove is not needed.
   * acceptMove can not be used.
   * See QMCHamiltonians::MomentumEstimator as an example
   */
  void makeVirtualMoves(const SingleParticlePos& newpos);

  /** move all the particles of a walker
   * @param awalker the walker to operate
   * @param deltaR proposed displacement
   * @param dt  factor of deltaR
   * @return true if all the moves are legal.
   *
   * If big displacements or illegal positions are detected, return false.
   * If all good, R = awalker.R + dt* deltaR
   */
  bool makeMoveAllParticles(const Walker_t& awalker, const ParticlePos& deltaR, RealType dt);

  bool makeMoveAllParticles(const Walker_t& awalker, const ParticlePos& deltaR, const std::vector<RealType>& dt);

  /** move all the particles including the drift
   *
   * Otherwise, everything is the same as makeMove for a walker
   */
  bool makeMoveAllParticlesWithDrift(const Walker_t& awalker,
                                     const ParticlePos& drift,
                                     const ParticlePos& deltaR,
                                     RealType dt);

  bool makeMoveAllParticlesWithDrift(const Walker_t& awalker,
                                     const ParticlePos& drift,
                                     const ParticlePos& deltaR,
                                     const std::vector<RealType>& dt);

  /** accept or reject a proposed move
   *  Two operation modes:
   *  The using and updating distance tables via `ParticleSet` operate in two modes, regular and forward modes.
   *
   *  Regular mode
   *  The regular mode can only be used when the distance tables for particle pairs are fully up-to-date.
   *  This is the case after calling `ParticleSet::update()` in a unit test or after p-by-p moves in a QMC driver.
   *  In this mode, the distance tables remain up-to-date after calling `ParticleSet::acceptMove`
   *  and calling `ParticleSet::rejectMove` is not mandatory.
   *
   *  Forward mode
   *  The forward mode assumes that distance table is not fully up-to-date until every particle is accepted
   *  or rejected to move once in order. This is the mode used in the p-by-p part of drivers.
   *  In this mode, calling `ParticleSet::accept_rejectMove` is required to handle accept/reject rather than
   *  calling individual `ParticleSet::acceptMove` and `ParticleSet::reject`.
   *  `ParticleSet::accept_rejectMove(iel)` ensures the distance tables (jel < iel) part is fully up-to-date
   *  regardless a move is accepted or rejected. For this reason, the rejecting operation inside
   *  `ParticleSet::accept_rejectMove` involves writing the distances with respect to the old particle position.
   */
  void accept_rejectMove(Index_t iat, bool accepted, bool forward_mode = true);

  /** accept the move and update the particle attribute by the proposed move in regular mode
   *@param iat the index of the particle whose position and other attributes to be updated
   */
  void acceptMove(Index_t iat);

  /** reject a proposed move in regular mode
   * @param iat the electron whose proposed move gets rejected.
   */
  void rejectMove(Index_t iat);
  /// batched version of acceptMove and rejectMove fused
  static void mw_accept_rejectMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                   Index_t iat,
                                   const std::vector<bool>& isAccepted,
                                   bool forward_mode = true);

  void initPropertyList();
  inline int addProperty(const std::string& pname) { return PropertyList.add(pname.c_str()); }

  int addPropertyHistory(int leng);
  //        void rejectedMove();
  //        void resetPropertyHistory( );
  //        void addPropertyHistoryPoint(int index, RealType data);

  void convert(const ParticlePos& pin, ParticlePos& pout);
  void convert2Unit(const ParticlePos& pin, ParticlePos& pout);
  void convert2Cart(const ParticlePos& pin, ParticlePos& pout);
  void convert2Unit(ParticlePos& pout);
  void convert2Cart(ParticlePos& pout);
  void convert2UnitInBox(const ParticlePos& pint, ParticlePos& pout);
  void convert2CartInBox(const ParticlePos& pint, ParticlePos& pout);

  void applyBC(const ParticlePos& pin, ParticlePos& pout);
  void applyBC(ParticlePos& pos);
  void applyBC(const ParticlePos& pin, ParticlePos& pout, int first, int last);
  void applyMinimumImage(ParticlePos& pinout);

  /** load a Walker_t to the current ParticleSet
   * @param awalker the reference to the walker to be loaded
   * @param pbyp true if it is used by PbyP update
   *
   * PbyP requires the distance tables and Sk with awalker.R
   */
  void loadWalker(Walker_t& awalker, bool pbyp);
  /** batched version of loadWalker */
  static void mw_loadWalker(const RefVectorWithLeader<ParticleSet>& p_list,
                            const RefVector<Walker_t>& walkers,
                            const std::vector<bool>& recompute,
                            bool pbyp);

  /** save this to awalker
   *
   *  just the R, G, and L
   *  More duplicate data that makes code difficult to reason about should be removed.
   */
  void saveWalker(Walker_t& awalker);

  /** batched version of saveWalker
   *
   *  just the R, G, and L
   */
  static void mw_saveWalker(const RefVectorWithLeader<ParticleSet>& psets, const RefVector<Walker_t>& walkers);

  /** update structure factor and unmark active_ptcl_
   *@param skip SK update if skipSK is true
   *
   * The Coulomb interaction evaluation needs the structure factor.
   * For these reason, call donePbyP after the loop of single
   * electron moves before evaluating the Hamiltonian. Unmark
   * active_ptcl_ is more of a safety measure probably not needed.
   */
  void donePbyP(bool skipSK = false);
  /// batched version of donePbyP
  static void mw_donePbyP(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK = false);

  ///return the address of the values of Hamiltonian terms
  inline FullPrecRealType* restrict getPropertyBase() { return Properties.data(); }

  ///return the address of the values of Hamiltonian terms
  inline const FullPrecRealType* restrict getPropertyBase() const { return Properties.data(); }

  ///return the address of the i-th properties
  inline FullPrecRealType* restrict getPropertyBase(int i) { return Properties[i]; }

  ///return the address of the i-th properties
  inline const FullPrecRealType* restrict getPropertyBase(int i) const { return Properties[i]; }

  inline void setTwist(SingleParticlePos& t) { myTwist = t; }
  inline SingleParticlePos getTwist() const { return myTwist; }

  /** Initialize particles around another ParticleSet
   * Used to initialize an electron ParticleSet by an ion ParticleSet
   */
  void randomizeFromSource(ParticleSet& src);

  /** get species name of particle i
   */
  inline const std::string& species_from_index(int i) { return my_species_.speciesName[GroupID[i]]; }

  inline size_t getTotalNum() const { return TotalNum; }

  inline void resize(size_t numPtcl)
  {
    TotalNum = numPtcl;

    R.resize(numPtcl);
    spins.resize(numPtcl);
    GroupID.resize(numPtcl);
    G.resize(numPtcl);
    L.resize(numPtcl);
    Mass.resize(numPtcl);
    Z.resize(numPtcl);

    coordinates_->resize(numPtcl);
  }

  inline void clear()
  {
    TotalNum = 0;

    R.clear();
    spins.clear();
    GroupID.clear();
    G.clear();
    L.clear();
    Mass.clear();
    Z.clear();

    coordinates_->resize(0);
  }

  inline void assign(const ParticleSet& ptclin)
  {
    resize(ptclin.getTotalNum());
    R.InUnit   = ptclin.R.InUnit;
    R          = ptclin.R;
    spins      = ptclin.spins;
    GroupID    = ptclin.GroupID;
    is_spinor_ = ptclin.is_spinor_;
    if (ptclin.SubPtcl.size())
    {
      SubPtcl.resize(ptclin.SubPtcl.size());
      SubPtcl = ptclin.SubPtcl;
    }
  }

  ///return the number of groups
  inline int groups() const { return SubPtcl.size() - 1; }

  ///return the first index of a group i
  inline int first(int igroup) const { return SubPtcl[igroup]; }

  ///return the last index of a group i
  inline int last(int igroup) const { return SubPtcl[igroup + 1]; }

  ///return the group id of a given particle in the particle set.
  inline int getGroupID(int iat) const
  {
    assert(iat >= 0 && iat < TotalNum);
    return GroupID[iat];
  }

  ///return the size of a group
  inline int groupsize(int igroup) const { return SubPtcl[igroup + 1] - SubPtcl[igroup]; }

  ///add attributes to list for IO
  template<typename ATList>
  inline void createAttributeList(ATList& AttribList)
  {
    R.setTypeName(ParticleTags::postype_tag);
    R.setObjName(ParticleTags::position_tag);
    spins.setTypeName(ParticleTags::scalartype_tag);
    spins.setObjName(ParticleTags::spins_tag);
    GroupID.setTypeName(ParticleTags::indextype_tag);
    GroupID.setObjName(ParticleTags::ionid_tag);
    //add basic attributes
    AttribList.add(R);
    AttribList.add(spins);
    AttribList.add(GroupID);

    G.setTypeName(ParticleTags::gradtype_tag);
    L.setTypeName(ParticleTags::laptype_tag);

    G.setObjName("grad");
    L.setObjName("lap");

    AttribList.add(G);
    AttribList.add(L);

    //more particle attributes
    Mass.setTypeName(ParticleTags::scalartype_tag);
    Mass.setObjName("mass");
    AttribList.add(Mass);

    Z.setTypeName(ParticleTags::scalartype_tag);
    Z.setObjName("charge");
    AttribList.add(Z);
  }

  inline int getNumDistTables() const { return DistTables.size(); }

  /// initialize a shared resource and hand it to a collection
  void createResource(ResourceCollection& collection) const;
  /** acquire external resource and assocaite it with the list of ParticleSet
   * Note: use RAII ResourceCollectionTeamLock whenever possible
   */
  static void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<ParticleSet>& p_list);
  /** release external resource
   * Note: use RAII ResourceCollectionTeamLock whenever possible
   */
  static void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<ParticleSet>& p_list);

  static RefVectorWithLeader<DistanceTable> extractDTRefList(const RefVectorWithLeader<ParticleSet>& p_list, int id);
  static RefVectorWithLeader<DynamicCoordinates> extractCoordsRefList(const RefVectorWithLeader<ParticleSet>& p_list);
  static RefVectorWithLeader<StructFact> extractSKRefList(const RefVectorWithLeader<ParticleSet>& p_list);

protected:
  /// reference to global simulation cell
  const SimulationCell& simulation_cell_;

  ///true if the particles are grouped
  bool is_grouped_;
  ///true if the particles have the same mass
  bool same_mass_;
  ///true is a dynamic spin calculation
  bool is_spinor_;
  /** the index of the active particle during particle-by-particle moves
   *
   * when a single particle move is proposed, the particle id is assigned to active_ptcl_
   * No matter the move is accepted or rejected, active_ptcl_ is marked back to -1.
   * This state flag is used for picking coordinates and distances for SPO evaluation.
   */
  Index_t active_ptcl_;
  ///the proposed position of active_ptcl_ during particle-by-particle moves
  SingleParticlePos active_pos_;
  ///the proposed spin of active_ptcl_ during particle-by-particle moves
  Scalar_t active_spin_val_;

  ///SpeciesSet of particles
  SpeciesSet my_species_;

  ///Structure factor
  std::unique_ptr<StructFact> structure_factor_;

  /** map to handle distance tables
   *
   * myDistTableMap[source-particle-tag]= locator in the distance table
   * myDistTableMap[ObjectTag] === 0
   */
  std::map<std::string, int> myDistTableMap;

  /// distance tables that need to be updated by moving this ParticleSet
  std::vector<std::unique_ptr<DistanceTable>> DistTables;

  /// Descriptions from distance table creation.  Same order as DistTables.
  std::vector<std::string> distTableDescriptions;

  TimerList_t myTimers;

  SingleParticlePos myTwist;

  std::string ParentName;

  ///total number of particles
  size_t TotalNum;

  ///array to handle a group of distinct particles per species
  ParticleIndex SubPtcl;
  ///internal representation of R. It can be an SoA copy of R
  std::unique_ptr<DynamicCoordinates> coordinates_;

  /** compute temporal DistTables and SK for a new particle position
   *
   * @param iat the particle that is moved on a sphere
   * @param newpos a new particle position
   * @param maybe_accept if false, the caller guarantees that the proposed move will not be accepted.
   */
  void computeNewPosDistTables(Index_t iat, const SingleParticlePos& newpos, bool maybe_accept = true);


  /** compute temporal DistTables and SK for a new particle position for each walker in a batch
   *
   * @param p_list the list of wrapped ParticleSet references in a walker batch
   * @param iat the particle that is moved on a sphere
   * @param new_positions new particle positions
   * @param maybe_accept if false, the caller guarantees that the proposed move will not be accepted.
   */
  static void mw_computeNewPosDistTables(const RefVectorWithLeader<ParticleSet>& p_list,
                                         Index_t iat,
                                         const std::vector<SingleParticlePos>& new_positions,
                                         bool maybe_accept = true);

  /** actual implemenation for accepting a proposed move in forward mode
   *
   * @param iat the index of the particle whose position and other attributes to be updated
   */
  void acceptMoveForwardMode(Index_t iat);

  /** reject a proposed move in forward mode
   * @param iat the electron whose proposed move gets rejected.
   */
  void rejectMoveForwardMode(Index_t iat);
};

} // namespace qmcplusplus
#endif
