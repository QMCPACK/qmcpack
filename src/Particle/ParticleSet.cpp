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
    ObjectTag = PtclObjectCounter;
    PtclObjectCounter++;
    G.setTypeName(ParticleTags::postype_tag);
    G.setObjName("grad");
    L.setTypeName(ParticleTags::scalartype_tag);
    L.setObjName("laplacian");
    addAttribute(G);
    addAttribute(L);
  }

  ParticleSet::~ParticleSet() {}

  ///write to a ostream
  bool ParticleSet::get(ostream& os) const {
    for(int i=0; i<LocalNum; i++) {
      os << Species.speciesName[GroupID[i]]  << R[i] << endl;
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


}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
