//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PARTICLESET_H
#define QMCPLUSPLUS_PARTICLESET_H

#include <Configuration.h>
#include <Particle/Walker.h>
#include <Utilities/SpeciesSet.h>
#include <Utilities/PooledData.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Utilities/NewTimer.h>

namespace qmcplusplus {

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
       , public ParticleBase<PtclOnLatticeTraits>
  {
  public:
    ///@typedef walker type
    typedef Walker<QMCTraits,PtclOnLatticeTraits> Walker_t;
    ///@typedef container type to store the property
    typedef Walker_t::PropertyContainer_t  PropertyContainer_t;
    ///@typedef buufer type for a serialized buffer
    typedef Walker_t::Buffer_t             Buffer_t;

    //@{ public data members
    ///property of an ensemble represented by this ParticleSet
    MCDataType<RealType> EnsembleProperty;

    ///gradients of the particles
    ParticleGradient_t G;

    ///laplacians of the particles
    ParticleLaplacian_t L;

    ///differential gradients of the particles
    ParticleGradient_t dG;

    ///differential laplacians of the particles
    ParticleLaplacian_t dL;

    ///current position after applying PBC in the Lattice Unit
    ParticlePos_t redR;

    /** ID map that reflects species group
     *
     * IsGrouped=true, if ID==IndirectID
     */
    ParticleIndex_t IndirectID;
    ///mass of each particle
    ParticleScalar_t Mass;

    ///true, if a physical or local bounding box is used
    bool UseBoundBox;
    ///true if fast update for sphere moves
    bool UseSphereUpdate;
    ///true if the particles are grouped
    bool IsGrouped;
    ///threa id
    Index_t ThreadID;
    ///the indexp of the active particle for particle-by-particle moves
    Index_t activePtcl;

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
    vector<DistanceTableData*> DistTables;

    ///spherical-grids for non-local PP
    vector<ParticlePos_t*> Sphere;

    ///Particle density in G-space for MPC interaction
    vector<TinyVector<int,OHMMS_DIM> > DensityReducedGvecs;
    vector<ComplexType>   Density_G;
    Array<RealType,OHMMS_DIM> Density_r;

    /// DFT potential
    vector<TinyVector<int,OHMMS_DIM> > VHXCReducedGvecs;
    vector<ComplexType>   VHXC_G[2];
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
    vector<ParticleSet*> myClones;

    ///Property history vector
    vector<vector<RealType> >  PropertyHistory;
    vector<int> PHindex;
    ///@} 

    ///default constructor
    ParticleSet();

    ///copy constructor
    ParticleSet(const ParticleSet& p);

    ///default destructor
    virtual ~ParticleSet();

    ///write to a ostream
    bool get(ostream& ) const;

    ///read from istream
    bool put(istream& );

    ///reset member data
    void reset();

    ///initialize ParticleSet from xmlNode
    bool put(xmlNodePtr cur);

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
    int  addTable(const ParticleSet& psrc);

    /** reset all the collectable quantities during a MC iteration 
     */
    inline void resetCollectables()
    {
      std::fill(Collectables.begin(),Collectables.end(),0.0);
    }

    /** update the internal data
     *@param iflag index for the update mode
     */
    void update(int iflag=0);

    /**update the internal data with new position
     *@param pos position vector assigned to R
     */
    void update(const ParticlePos_t& pos);

    /** create Structure Factor with PBCs
     */
    void createSK();

    ///retrun the SpeciesSet of this particle set
    inline SpeciesSet& getSpeciesSet() { return mySpecies;}
    ///retrun the const SpeciesSet of this particle set
    inline const SpeciesSet& getSpeciesSet() const { return mySpecies;}

    ///return this id
    inline int tag() const { return ObjectTag;}

    ///return parent's id 
    inline int parent() const { return ParentTag;}

    inline RealType getTotalWeight() const 
    {
      return EnsembleProperty.Weight;
    }

    void resetGroups();

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
    /** move all the particles including the drift
     *
     * Otherwise, everything is the same as makeMove for a walker
     */
    bool makeMoveWithDrift(const Walker_t& awalker
        , const ParticlePos_t& drift, const ParticlePos_t& deltaR, RealType dt);

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

    inline SingleParticlePos_t getOldPos() const
    {
      return activePos;
    }

    void initPropertyList();
    inline int addProperty(const string& pname) {
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

    //void registerData(Buffer_t& buf);
    //void registerData(Walker_t& awalker, Buffer_t& buf);
    //void updateBuffer(Walker_t& awalker, Buffer_t& buf);
    //void updateBuffer(Buffer_t& buf);
    //void copyToBuffer(Buffer_t& buf);
    //void copyFromBuffer(Buffer_t& buf);

    //return the address of the values of Hamiltonian terms
    inline RealType* restrict getPropertyBase() 
    {
      return Properties.data();
    }

    //return the address of the values of Hamiltonian terms
    inline const RealType* restrict getPropertyBase() const 
    {
      return Properties.data();
    }

    ///return the address of the i-th properties
    inline RealType* restrict getPropertyBase(int i) 
    {
      return Properties[i];
    }

    ///return the address of the i-th properties
    inline const RealType* restrict getPropertyBase(int i) const 
    {
      return Properties[i];
    }
    
    inline void setTwist(SingleParticlePos_t t){myTwist=t;}
    inline SingleParticlePos_t getTwist(){return myTwist;}

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

  protected:
    ///the number of particle objects
    static Index_t PtclObjectCounter;

    ///id of this object    
    Index_t ObjectTag;

    ///id of the parent
    Index_t ParentTag;
    
    /** map to handle distance tables
     *
     * myDistTableMap[source-particle-tag]= locator in the distance table
     * myDistTableMap[ObjectTag] === 0
     */
    map<int,int> myDistTableMap;
    void initParticleSet();

    vector<NewTimer*> myTimers;
    SingleParticlePos_t myTwist;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
