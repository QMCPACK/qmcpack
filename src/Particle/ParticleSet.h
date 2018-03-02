//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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

#include <Configuration.h>
#include <ParticleTags.h>
#include <Particle/Walker.h>
#include <Utilities/SpeciesSet.h>
#include <Utilities/PooledData.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Utilities/NewTimer.h>
#include <OhmmsSoA/Container.h>

namespace qmcplusplus
{

///forward declaration of DistanceTableData
class DistanceTableData;

class StructFact;

/** Monte Carlo Data of an ensemble
 *
 * The quantities are shared by all the nodes in a group
 * - NumSamples number of samples
 * - Weight     total weight of a sample
 * - Energy     average energy of a sample
 * - Variance   variance
 * - LivingFraction fraction of walkers alive each step.
 */
template<typename T>
struct MCDataType
{
  T NumSamples;
  T RNSamples;
  T Weight;
  T Energy;
  T AlternateEnergy;
  T Variance;
  T R2Accepted;
  T R2Proposed;
  T LivingFraction;
};


/** Specialized paritlce class for atomistic simulations
 *
 * Derived from QMCTraits, ParticleBase<PtclOnLatticeTraits> and OhmmsElementBase.
 * The ParticleLayout class represents a supercell with/without periodic boundary
 * conditions. The ParticleLayout class also takes care of spatial decompositions
 * for efficient evaluations for the interactions with a finite cutoff.
 */
class ParticleSet
  :  public QMCTraits
  , public OhmmsElementBase
  , public PtclOnLatticeTraits
{
public:
  ///@typedef walker type
  typedef Walker<QMCTraits,PtclOnLatticeTraits>  Walker_t;
  ///@typedef container type to store the property
  typedef Walker_t::PropertyContainer_t          PropertyContainer_t;
  ///@typedef buffer type for a serialized buffer
  typedef PooledData<RealType>                   Buffer_t;

  enum quantum_domains {no_quantum_domain=0,classical,quantum};

  ///quantum_domain of the particles, default = classical
  quantum_domains quantum_domain;

  //@{ public data members
  ///property of an ensemble represented by this ParticleSet
  MCDataType<EstimatorRealType> EnsembleProperty;

  ///ParticleLayout
  ParticleLayout_t Lattice, PrimitiveLattice;
  ///Long-range box
  ParticleLayout_t LRBox;

  ///unique, persistent ID for each particle
  ParticleIndex_t ID;
  ///index to the primitice cell with tiling
  ParticleIndex_t PCID;
  /** ID map that reflects species group
   *
   * IsGrouped=true, if ID==IndirectID
   */
  ParticleIndex_t IndirectID;
  ///Species ID
  ParticleIndex_t GroupID;
  ///Position
  ParticlePos_t R;
  ///SoA copy of R
  VectorSoaContainer<RealType,DIM> RSoA;
  ///gradients of the particles
  ParticleGradient_t G;
  ///laplacians of the particles
  ParticleLaplacian_t L;
  ///differential gradients of the particles
  ParticleGradient_t dG;
  ///differential laplacians of the particles
  ParticleLaplacian_t dL;
  ///mass of each particle
  ParticleScalar_t Mass;
  ///charge of each particle
  ParticleScalar_t Z;

  ///true, if a physical or local bounding box is used
  bool UseBoundBox;
  ///true if fast update for sphere moves
  bool UseSphereUpdate;
  ///true if the particles are grouped
  bool IsGrouped;
  ///true if the particles have the same mass
  bool SameMass;
  /// true if all the internal state is ready for estimators
  bool Ready4Measure;
  ///threa id
  Index_t ThreadID;
  ///the index of the active particle for particle-by-particle moves
  Index_t activePtcl;
  ///the group of the active particle for particle-by-particle moves
  Index_t activeGroup;
  ///the index of the active bead for particle-by-particle moves
  Index_t activeBead;
  ///the direction reptile traveling
  Index_t direction;

  /** the position of the active particle for particle-by-particle moves
   *
   * Saves the position before making a move to handle rejectMove
   */
  SingleParticlePos_t activePos;

  /** the proposed position in the Lattice unit
   */
  SingleParticlePos_t newRedPos;

  ///SpeciesSet of particles
  SpeciesSet mySpecies;

  ///Structure factor
  StructFact *SK;

  ///distance tables that need to be updated by moving this ParticleSet
  std::vector<DistanceTableData*> DistTables;

  ///spherical-grids for non-local PP
  std::vector<ParticlePos_t*> Sphere;

  ///Particle density in G-space for MPC interaction
  std::vector<TinyVector<int,OHMMS_DIM> > DensityReducedGvecs;
  std::vector<ComplexType>   Density_G;
  Array<RealType,OHMMS_DIM> Density_r;

  /// DFT potential
  std::vector<TinyVector<int,OHMMS_DIM> > VHXCReducedGvecs;
  std::vector<ComplexType>   VHXC_G[2];
  Array<RealType,OHMMS_DIM> VHXC_r[2];

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
  PropertyContainer_t  Properties;

  /** observables in addition to those registered in Properties/PropertyList
   *
   * Such observables as density, gofr, sk are not stored per walker but
   * collected during QMC iterations.
   */
  Buffer_t Collectables;

  ///clones of this object: used by the thread pool
  std::vector<ParticleSet*> myClones;

  ///Property history vector
  std::vector<std::vector<EstimatorRealType> >  PropertyHistory;
  std::vector<int> PHindex;
  ///@}

  ///current MC step
  int current_step;

  ///default constructor
  ParticleSet();

  ///copy constructor
  ParticleSet(const ParticleSet& p);

  ///default destructor
  virtual ~ParticleSet();

  /** create  particles
   * @param n number of particles
   */
  void create(int n);
  /** create grouped particles
   * @param agroup number of particles per group
   */
  void create(const std::vector<int>& agroup);

  ///write to a std::ostream
  bool get(std::ostream& ) const;

  ///read from std::istream
  bool put( std::istream& );

  ///reset member data
  void reset();

  ///initialize ParticleSet from xmlNode
  bool put(xmlNodePtr cur);

  ///specify quantum_domain of particles
  void set_quantum_domain(quantum_domains qdomain);

  void set_quantum()
  {
    quantum_domain=quantum;
  }

  inline bool is_classical() const
  {
    return quantum_domain==classical;
  }

  inline bool is_quantum() const
  {
    return quantum_domain==quantum;
  }

  ///check whether quantum domain is valid for particles
  inline bool quantum_domain_valid(quantum_domains qdomain) const
  {
    return qdomain!=no_quantum_domain;
  }

  ///check whether quantum domain is valid for particles
  inline bool quantum_domain_valid() const
  {
    return quantum_domain_valid(quantum_domain);
  }

  ///set UseBoundBox
  void setBoundBox(bool yes);

  /** check bounding box
   * @param rb cutoff radius to check the condition
   */
  void checkBoundBox(RealType rb);

  /**  add a distance table
   * @param psrc source particle set
   *
   * Ensure that the distance for this-this is always created first.
   */
  int  addTable(const ParticleSet& psrc, int dt_type);

  /** returns index of a distance table, -1 if not present
   * @param psrc source particle set
   */
  int getTable(const ParticleSet& psrc);

  /** reset all the collectable quantities during a MC iteration
   */
  inline void resetCollectables()
  {
    std::fill(Collectables.begin(),Collectables.end(),0.0);
  }

  /** update the internal data
   *@param skip SK update if skipSK is true
   */
  void update(bool skipSK=false);

  /**update the internal data with new position
   *@param pos position vector assigned to R
   */
  void update(const ParticlePos_t& pos);

  /** create Structure Factor with PBCs
   */
  void createSK();

  /** Turn on per particle storage in Structure Factor
   */
  void turnOnPerParticleSK();

  ///retrun the SpeciesSet of this particle set
  inline SpeciesSet& getSpeciesSet()
  {
    return mySpecies;
  }
  ///retrun the const SpeciesSet of this particle set
  inline const SpeciesSet& getSpeciesSet() const
  {
    return mySpecies;
  }

  ///return parent's name
  inline const std::string& parentName() const
  {
    return ParentName;
  }
  inline void setName(const std::string& aname)
  {
    myName     = aname;
    if(ParentName=="0")
    {
      ParentName = aname;
    }
  }

  inline RealType getTotalWeight() const
  {
    return EnsembleProperty.Weight;
  }

  void resetGroups();

  /** set active particle
   * @param iat particle index
   *
   * Compute internal data based on current R[iat]
   * Introduced to work with update-only methods.
   */
  void setActive(int iat);

  /** return the position of the active partice
   *
   * activePtcl=-1 is used to flag non-physical moves
   */
  inline const PosType& activeR(int iat) const
  {
    return (activePtcl == iat)? activePos:R[iat];
  }
  
  /**move a particle
   *@param iat the index of the particle to be moved
   *@param displ random displacement of the iat-th particle
   *
   * Update activePos  by  R[iat]+displ
   */
  SingleParticlePos_t makeMove(Index_t iat, const SingleParticlePos_t& displ);

  /** move a particle
   * @param iat the index of the particle to be moved
   * @param displ random displacement of the iat-th particle
   * @return true, if the move is valid
   */
  bool makeMoveAndCheck(Index_t iat, const SingleParticlePos_t& displ);

  /** move all the particles of a walker
   * @param awalker the walker to operate
   * @param deltaR proposed displacement
   * @param dt  factor of deltaR
   * @return true if all the moves are legal.
   *
   * If big displacements or illegal positions are detected, return false.
   * If all good, R = awalker.R + dt* deltaR
   */
  bool makeMove(const Walker_t& awalker, const ParticlePos_t& deltaR, RealType dt);

  bool makeMove(const Walker_t& awalker, const ParticlePos_t& deltaR, const std::vector<RealType>& dt);
  /** move all the particles including the drift
   *
   * Otherwise, everything is the same as makeMove for a walker
   */
  bool makeMoveWithDrift(const Walker_t& awalker
                         , const ParticlePos_t& drift, const ParticlePos_t& deltaR, RealType dt);

  bool makeMoveWithDrift(const Walker_t& awalker
                         , const ParticlePos_t& drift, const ParticlePos_t& deltaR, const std::vector<RealType>& dt);

  void makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ);

  /** Handles a virtual move for all the particles to ru.
   * @param ru position in the reduced cordinate
   *
   * The data of the 0-th particle is overwritten by the new position
   * and the rejectMove should be called for correct use.
   * See QMCHamiltonians::MomentumEstimator
   */
  void makeVirtualMoves(const SingleParticlePos_t& newpos);

  /** accept the move
   *@param iat the index of the particle whose position and other attributes to be updated
   */
  void acceptMove(Index_t iat);

  /** reject the move
   */
  void rejectMove(Index_t iat);

  void initPropertyList();
  inline int addProperty(const std::string& pname)
  {
    return PropertyList.add(pname.c_str());
  }

  int addPropertyHistory(int leng);
  //        void rejectedMove();
  //        void resetPropertyHistory( );
  //        void addPropertyHistoryPoint(int index, RealType data);

  void clearDistanceTables();
  void resizeSphere(int nc);

  void convert(const ParticlePos_t& pin, ParticlePos_t& pout);
  void convert2Unit(const ParticlePos_t& pin, ParticlePos_t& pout);
  void convert2Cart(const ParticlePos_t& pin, ParticlePos_t& pout);
  void convert2Unit(ParticlePos_t& pout);
  void convert2Cart(ParticlePos_t& pout);
  void convert2UnitInBox(const ParticlePos_t& pint, ParticlePos_t& pout);
  void convert2CartInBox(const ParticlePos_t& pint, ParticlePos_t& pout);

  void applyBC(const ParticlePos_t& pin, ParticlePos_t& pout);
  void applyBC(ParticlePos_t& pos);
  void applyBC(const ParticlePos_t& pin, ParticlePos_t& pout, int first, int last);
  void applyMinimumImage(ParticlePos_t& pinout);

  /** load a Walker_t to the current ParticleSet
   * @param awalker the reference to the walker to be loaded
   * @param pbyp true if it is used by PbyP update
   *
   * PbyP requires the distance tables and Sk with awalker.R
   */
  void loadWalker(Walker_t& awalker, bool pbyp);
  /** save this to awalker
   */
  void saveWalker(Walker_t& awalker);

  /** update the buffer
   *@param skip SK update if skipSK is true
   */
  void donePbyP(bool skipSK=false);

  //return the address of the values of Hamiltonian terms
  inline EstimatorRealType* restrict getPropertyBase()
  {
    return Properties.data();
  }

  //return the address of the values of Hamiltonian terms
  inline const EstimatorRealType* restrict getPropertyBase() const
  {
    return Properties.data();
  }

  ///return the address of the i-th properties
  inline EstimatorRealType* restrict getPropertyBase(int i)
  {
    return Properties[i];
  }

  ///return the address of the i-th properties
  inline const EstimatorRealType* restrict getPropertyBase(int i) const
  {
    return Properties[i];
  }

  inline void setTwist(SingleParticlePos_t& t)
  {
    myTwist=t;
  }
  inline SingleParticlePos_t getTwist() const
  {
    return myTwist;
  }

  /** Initialize particles around another ParticleSet
   * Used to initialize an electron ParticleSet by an ion ParticleSet
   */
  void randomizeFromSource (ParticleSet &src);

  /** make clones
   * @param n number of clones including itself
   */
  virtual void make_clones(int n);

  /** return the ip-th clone
   * @param ip thread number
   *
   * Return itself if ip==0
   */
  inline ParticleSet* get_clone(int ip)
  {
    if(ip >= myClones.size())
      return 0;
    return (ip)? myClones[ip]:this;
  }

  inline const ParticleSet* get_clone(int ip) const
  {
    if(ip >= myClones.size())
      return 0;
    return (ip)? myClones[ip]:this;
  }

  inline int clones_size() const
  {
    return myClones.size();
  }

  /** update R of its own and its clones
   * @param rnew new position array of N
   */
  template<typename PAT>
  inline void update_clones(const PAT& rnew)
  {
    if(R.size() != rnew.size())
      APP_ABORT("ParticleSet::updateR failed due to different sizes");
    R=rnew;
    for(int ip=1; ip<myClones.size(); ++ip)
      myClones[ip]->R=rnew;
  }

  /** reset internal data of clones including itself
   */
  void reset_clones();
      
  /** get species name of particle i
   */
  inline const std::string& species_from_index(int i)
  {
    return mySpecies.speciesName[GroupID[i]];
  }

  inline int getTotalNum() const
  {
    return TotalNum;
  }

  inline void resize(int numPtcl)
  {
    TotalNum = numPtcl;

    R.resize(numPtcl);
    ID.resize(numPtcl);
    PCID.resize(numPtcl);
    GroupID.resize(numPtcl);
    G.resize(numPtcl);
    dG.resize(numPtcl);
    L.resize(numPtcl);
    dL.resize(numPtcl);
    Mass.resize(numPtcl);
    Z.resize(numPtcl);
    IndirectID.resize(numPtcl);

    RSoA.resize(numPtcl);
  }

  inline void clear()
  {
    TotalNum = 0;

    R.clear();
    ID.clear();
    PCID.clear();
    GroupID.clear();
    G.clear();
    dG.clear();
    L.clear();
    dL.clear();
    Mass.clear();
    Z.clear();
    IndirectID.clear();

    RSoA.resize(0);
  }

  inline void assign(const ParticleSet& ptclin)
  {
    resize(ptclin.getTotalNum());
    Lattice = ptclin.Lattice;
    PrimitiveLattice = ptclin.PrimitiveLattice;
    R.InUnit = ptclin.R.InUnit;
    R = ptclin.R;
    ID = ptclin.ID;
    GroupID = ptclin.GroupID;
    if(ptclin.SubPtcl.size())
    {
      SubPtcl.resize(ptclin.SubPtcl.size());
      SubPtcl =ptclin.SubPtcl;
    }
  }

  ///return the number of groups
  inline int groups() const
  {
    return SubPtcl.size()-1;
  }

  ///return the first index of a group i
  inline int first(int igroup) const
  {
    return SubPtcl[igroup];
  }

  ///return the last index of a group i
  inline int last(int igroup) const
  {
    return SubPtcl[igroup+1];
  }

  ///add attributes to list for IO
  template<typename ATList>
  inline void createAttributeList(ATList& AttribList)
  {
    R.setTypeName(ParticleTags::postype_tag);
    R.setObjName(ParticleTags::position_tag);
    ID.setTypeName(ParticleTags::indextype_tag);
    ID.setObjName(ParticleTags::id_tag);
    GroupID.setTypeName(ParticleTags::indextype_tag);
    GroupID.setObjName(ParticleTags::ionid_tag);
    //add basic attributes
    AttribList.add(R);
    AttribList.add(ID);
    AttribList.add(GroupID);

    G.setTypeName(ParticleTags::gradtype_tag);
    L.setTypeName(ParticleTags::laptype_tag);
    dG.setTypeName(ParticleTags::gradtype_tag);
    dL.setTypeName(ParticleTags::laptype_tag);

    G.setObjName("grad");
    L.setObjName("lap");
    dG.setObjName("dgrad");
    dL.setObjName("dlap");

    AttribList.add(G);
    AttribList.add(L);
    AttribList.add(dG);
    AttribList.add(dL);

    //more particle attributes
    Mass.setTypeName(ParticleTags::scalartype_tag);
    Mass.setObjName("mass");
    AttribList.add(Mass);

    Z.setTypeName(ParticleTags::scalartype_tag);
    Z.setObjName("charge");
    AttribList.add(Z);

    PCID.setTypeName(ParticleTags::indextype_tag); //add PCID tags
    PCID.setObjName("pcid");
    AttribList.add(PCID);

    IndirectID.setTypeName(ParticleTags::indextype_tag); //add IndirectID tags
    IndirectID.setObjName("id1");
    AttribList.add(IndirectID);
  }

protected:
  /** map to handle distance tables
   *
   * myDistTableMap[source-particle-tag]= locator in the distance table
   * myDistTableMap[ObjectTag] === 0
   */
  std::map<std::string,int> myDistTableMap;

  std::vector<NewTimer*> myTimers;
  SingleParticlePos_t myTwist;

  std::string ParentName;

  ///total number of particles
  int TotalNum;

  ///array to handle a group of distinct particles per species
  ParticleIndex_t SubPtcl;
};
}
#endif
