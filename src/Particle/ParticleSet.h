//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PARTICLESET_H
#define QMCPLUSPLUS_PARTICLESET_H

#include "Configuration.h"
#include "Utilities/SpeciesSet.h"
#include "Particle/Walker.h"
#include "OhmmsData/RecordProperty.h"

namespace qmcplusplus {

  ///forward declaration of DistanceTableData
  class DistanceTableData;  

  class StructFact;

  /** Specialized paritlce class for atomistic simulations
   *
   *Derived from QMCTraits, ParticleBase<PtclOnLatticeTraits> and OhmmsElementBase.
   *The ParticleLayout class represents a supercell with/without periodic boundary
   *conditions. The ParticleLayout class also takes care of spatial decompositions
   *for efficient evaluations for the interactions with a finite cutoff.
   */
  class ParticleSet:  
    public QMCTraits,
    public OhmmsElementBase,
    public ParticleBase<PtclOnLatticeTraits>
  {
   
  public:
    
    typedef ParticleAttrib<GradType>  ParticleGradient_t;
    typedef ParticleAttrib<ValueType> ParticleLaplacian_t;
    typedef Walker<RealType,ParticlePos_t> Walker_t;
    
    ///gradients of the particles
    ParticleGradient_t G;
    
    ///laplacians of the particles
    ParticleLaplacian_t L;
   
    ///SpeciesSet of particles
    SpeciesSet mySpecies;

    ///Structure factor
    StructFact *SK, *SKOld;

    ///spherical-grids for non-local PP
    vector<ParticlePos_t*> Sphere;

    ///default constructor
    ParticleSet();

    ///default constructor
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
    
    /** set the update mode
     * @param updatemode
     */
    //void setUpdateMode(int updatenode);

    /** add a distance table */
    void addTable(DistanceTableData* d_table);

    /**update the internal data
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

    inline SpeciesSet& getSpeciesSet() { return mySpecies;}
    inline const SpeciesSet& getSpeciesSet() const { return mySpecies;}

    ///return the id
    inline int tag() const { return ObjectTag;}
    
    /**move a particle
     *@param iat the index of the particle to be moved
     *@param displ random displacement of the iat-th particle
     *
     * Update activePos  by  R[iat]+displ
     */
    SingleParticlePos_t makeMove(Index_t iat, const SingleParticlePos_t& displ);

    void makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ);

    /** accept the move
     *@param iat the index of the particle whose position and other attributes to be updated
     */
    void acceptMove(Index_t iat);

    /** reject the move
     */
    void rejectMove(Index_t iat);

    void initPropertyList();
    inline int addProperty(const string& pname) {
      return PropertyList.add(pname.c_str());
    }

   void resizeSphere(int nc);

   void convert(const ParticlePos_t& pin, ParticlePos_t& pout);
   void convert2Unit(const ParticlePos_t& pin, ParticlePos_t& pout);
   void convert2Cart(const ParticlePos_t& pin, ParticlePos_t& pout);
   void convert2Unit(ParticlePos_t& pout);
   void convert2Cart(ParticlePos_t& pout);

   void applyBC(const ParticlePos_t& pin, ParticlePos_t& pout);
   void applyBC(ParticlePos_t& pos);
   void applyBC(const ParticlePos_t& pin, ParticlePos_t& pout, int first, int last);

   void registerData(PooledData<RealType>& buf);
   void registerData(Walker_t& awalker, PooledData<RealType>& buf);
   void updateBuffer(Walker_t& awalker, PooledData<RealType>& buf);
   void copyToBuffer(PooledData<RealType>& buf);
   void copyFromBuffer(PooledData<RealType>& buf);

  protected:
    ///the number of particle objects
    static Index_t PtclObjectCounter;

    ///id of this object    
    Index_t ObjectTag;

    ///the indexp of the active particle for particle-by-particle moves
    Index_t activePtcl;

    ///the position of the active particle for particle-by-particle moves
    SingleParticlePos_t activePos;

    ///distance tables that need to be updated by moving this ParticleSet
    vector<DistanceTableData*> DistTables;

    /** name-value map of Walker Properties
     *
     * PropertyMap is used to keep the name-value mapping of
     * Walker_t::Properties.
     */ 
    RecordNamedProperty<RealType> PropertyList;

    void initParticleSet();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
