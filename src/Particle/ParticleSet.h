/////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_QMC_PARTICLESET_H
#define OHMMS_QMC_PARTICLESET_H

#include "Configuration.h"
#include "Utilities/OhmmsObject.h"
#include "Utilities/SpeciesSet.h"
#include "Particle/Walker.h"

namespace ohmmsqmc {

  ///forward declaration of DistanceTableData
  class DistanceTableData;  

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
    
    Walker<RealType,ParticlePos_t>::PropertyContainer_t  Properties;

    ///gradients of the particles
    ParticleGradient_t G;
    
    ///laplacians of the particles
    ParticleLaplacian_t L;
   
    ///SpeciesSet of particles
    SpeciesSet Species;
    
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
    
    /**update the internal data
     *@param iflag index for the update mode
     */
    void update(int iflag=0);

    ///return the id
    inline int tag() const { return ObjectTag;}

    /**move a particle
     *@param iat the index of the particle to be moved
     *@param newpos new position of the iat-th particle
     */
    SingleParticlePos_t makeMove(Index_t iat, const SingleParticlePos_t& displ);

    /**accept the move
     *@param iat the index of the particle whose position and other attributes to be updated
     */
    void acceptMove(Index_t iat);

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

    void initParticleSet();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
