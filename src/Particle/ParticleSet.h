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
#ifndef OHMMS_QMC_PARTICLESET_H
#define OHMMS_QMC_PARTICLESET_H

#include "Configuration.h"
#include "Utilities/OhmmsObject.h"
#include "Utilities/SpeciesSet.h"

namespace ohmmsqmc {
  
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

   ///gradients of the particles
   ParticleGradient_t G;
   
   ///laplacians of the particles
   ParticleLaplacian_t L;
   
   ///SpeciesSet of particles
   SpeciesSet Species;
   
   ///default constructor
   ParticleSet();
   
   ///default destructor
   virtual ~ParticleSet();

//    /**function to create N-particle system 
//      *@param n an integer array containing the number of particles belonging
//      *to each subgroup.
//      *@brief Allocate the ParticleAttributes, such as R, G, L. 
//      *The size of n is the number of distinct subgroups. SubPtcl class
//      *is used to efficient evaluate the subgroup properties.
//      */
//    template<class IV>
//     void create(IV n) {
//       SubPtcl.resize(n.size()+1);
//       SubPtcl[0] = 0;
//       for(int is=0; is<n.size(); is++) SubPtcl[is+1] = SubPtcl[is]+n[is];
//       makeGroups();
//     }
   
   ///write to a ostream
   bool get(ostream& ) const;
   
   ///read from istream
   bool put(istream& );
   
   ///reset member data
   void reset();
   
   ///initialize ParticleSet from xmlNode
   bool put(xmlNodePtr cur);

   void update(int i=0);

  ///return the id
  inline int tag() const { return ObjectTag;}

 protected:

   int ObjectTag;
   static int PtclObjectCounter;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
