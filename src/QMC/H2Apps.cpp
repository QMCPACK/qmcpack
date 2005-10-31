//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMC/H2Apps.h"
#include "QMCBase/RandomSeqGenerator.h"
#include "QMCHamiltonians/HartreePotential.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/HarmonicPotential.h"
#include "QMCHamiltonians/BareKineticEnergy.h"

namespace qmcplusplus {

  using namespace xmlpp;

  bool H2Apps::init() {

    setParticleSets(d_root);

    setWavefunctions(d_root);
    
    setHamiltonian(d_root);
  }   

  bool  H2Apps::setParticleSets(xmlpp::Node* root) {
    
    el.m_name = "e";
    int iu = el.Species.addSpecies("u");
    int id = el.Species.addSpecies("d");
    int icharge = el.Species.addAttribute("charge");
    el.Species(icharge,iu) = -1.0;
    el.Species(icharge,id) = -1.0;
    vector<int> num(2,1);
    el.create(num);  

    ion.m_name = "i";
    icharge = ion.Species.addAttribute("charge");
    ion.Species(icharge,0) = 1.0;

    ion.create(2);
    ion.R[0] = PosType(0.0,0.0,-1.0);
    ion.R[1] = PosType(0.0,0.0,1.0);
  }

  /*@fn  bool H2Apps::setHamiltonian(xmlpp::Node* root)
   *@brief initialize the Hamiltonian
   */
  bool H2Apps::setHamiltonian(xmlpp::Node* root){

    H.add(new HartreePotential(el));
    H.add(new CoulombPotential(ion,el));
    H.add(new BareKineticEnergy);
    if(ion.getTotalNum()>1) H.add(new IonIonPotential(ion));
    return true;
  }
  
  bool H2Apps::setWavefunctions(xmlpp::Node* root) {

    ///to be added.
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
