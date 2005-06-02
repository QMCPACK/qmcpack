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
#include <numeric>
#include <iomanip>
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
namespace ohmmsqmc {

  int  ParticleSet::PtclObjectCounter = 0;

  ParticleSet::ParticleSet() {
    initParticleSet();
  }

  ParticleSet::ParticleSet(const ParticleSet& p) {
    initBase();
    initParticleSet();
    assign(p);
  }

  ParticleSet::~ParticleSet() {}

  void ParticleSet::initParticleSet() {
    ObjectTag = PtclObjectCounter;
    PtclObjectCounter++;
    G.setTypeName(ParticleTags::postype_tag);
    G.setObjName("grad");
    L.setTypeName(ParticleTags::scalartype_tag);
    L.setObjName("laplacian");
    addAttribute(G);
    addAttribute(L);
  }

  ///write to a ostream
  bool ParticleSet::get(ostream& os) const {
    for(int i=0; i<LocalNum; i++) {
      os << mySpecies.speciesName[GroupID[i]]  << R[i] << endl;
    }
    os << "sub-particle " << SubPtcl.size() << endl;
    for(int i=0; i<SubPtcl.size(); i++) os << SubPtcl[i] << " ";
    os << endl;
    return true;
  }
    
  ///read from istream
  bool ParticleSet::put(istream& is) { 
    return true;
  }

  ///reset member data
  void ParticleSet::reset() { 
  }

  ///read the particleset
  bool ParticleSet::put(xmlNodePtr cur){
    return true;
  }

  void ParticleSet::update(int iflag) { 
    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
  }  

  /** move a particle iat
   *@param iat the index of the particle to be moved
   *@param displ the displacement of the iath-particle position
   *@return the proposed position
   *Update activePtcl index and activePos position for the proposed move.
   *Evaluate the related distance table data DistanceTableData::Temp.
   */
  ParticleSet::SingleParticlePos_t 
  ParticleSet::makeMove(Index_t iat, const SingleParticlePos_t& displ) {
    activePtcl=iat;
    activePos=R[iat]+displ;
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->move(*this,activePos,iat);
    }
    return activePos;
  }
  
  /** update the particle attribute by the proposed move
   *@param iat the particle index
   *
   *When the activePtcl is equal to iat, overwrite the position and update the
   *content of the distance tables.
   */
  void ParticleSet::acceptMove(Index_t iat) {
    if(iat == activePtcl) {
      R[iat]=activePos; 
      for(int i=0; i< DistTables.size(); i++) {
        DistTables[i]->update(iat);
      }
    }
  }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
