//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
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
#include "QMCHamiltonians/CorePolPotential.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"

namespace ohmmsqmc {
 
  GeCorePolPotential::GeCorePolPotential(ParticleSet& ions, ParticleSet& els): 
    d_ie(NULL), alpha(0.3558), r_b(0.7048), eCoreCore(0.0) { 
  
    //set the distance tables
    d_ie = DistanceTable::getTable(DistanceTable::add(ions,els));
    nCenters = ions.getTotalNum();
    nParticles = els.getTotalNum();
    C = -0.5*alpha;
    r_binv = 1.0/r_b;
  
    CoreCoef.resize(nCenters);
    CoreCoreDipole.resize(nCenters,nCenters);
    ElCoreDipole.resize(nCenters,d_ie->getTotNadj());
  
    CoreCoreDipole = 0.0;
    ElCoreDipole = 0.0;
  
    //only calculate the cpp for Ge atoms
    int GeCounter = 0;
    for(int iat=0; iat<nCenters; iat++){
      string sname = ions.Species.speciesName[ions.GroupID[iat]];
      if(sname == "Ge"){
	LOGMSG("Adding a core-electron potential for " << sname << " #" << GeCounter++)
	  CoreCoef[iat] = true;
      }
      else CoreCoef[iat] = false;
    }
  
    //index for attribute charge
    int iz = ions.Species.addAttribute("charge");

    //calculate the Core-Core Dipole matrix
    int nn=0;
    for(int iat=0; iat<nCenters; iat++) {
      for(int jat=iat+1;jat<nCenters; jat++, nn++) {
	//check to see if both ions are Ge
	if(CoreCoef[iat]*CoreCoef[jat]){
	  RealType rinv3 = pow(d_ii->rinv(nn),3);//(1/R^3)
	  PosType dipole = rinv3*d_ii->dr(nn);//(\vec{R}/R^3)
	  CoreCoreDipole(iat,jat) = dipole*ions.Species(iz,ions.GroupID[jat]);//charge of jat
	  CoreCoreDipole(jat,iat) = dipole*ions.Species(iz,ions.GroupID[iat]);//charge of iat
	}
      }
    }

    //calculate the core-core term (constant)
    nn = 0;
    for(int iat=0; iat<nCenters; iat++) {
      for(int jat=iat+1;jat<nCenters; jat++, nn++) {
	eCoreCore += dot(CoreCoreDipole(iat,jat),CoreCoreDipole(iat,jat));
	eCoreCore += dot(CoreCoreDipole(jat,iat),CoreCoreDipole(jat,iat));
      }
    }
    eCoreCore*=C;

  }
    
  GeCorePolPotential::~GeCorePolPotential() { }

  GeCorePolPotential::ValueType 
  GeCorePolPotential::evaluate(ParticleSet& P) {

    RealType esum=0.0;
    //calculate the Electron-Core Dipole matrix
    int nn=0;
    for(int iat=0; iat<nCenters; iat++) {
      if(CoreCoef[iat]){
	for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++){
	  RealType rinv3 = pow(d_ie->rinv(nn),3);//(1/r^3)
	  PosType dipole = rinv3*d_ie->dr(nn);//(\vec{r}/r^3)
	  ElCoreDipole(iat,nn) = dipole*fcpp(d_ie->r(nn)*r_binv);
	}
      }
    }
  
    //now loop over the ions
    for(int iat=0; iat<nCenters; iat++) {
      //loop over the electrons
      for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++)
	esum += dot(ElCoreDipole(iat,nn),ElCoreDipole(iat,nn));
    
      //loop over distinct pairs of electrons
      for(int nnj=d_ie->M[iat]; nnj<d_ie->M[iat+1]; nnj++){
	for(int nnk=nnj+1; nnk<d_ie->M[iat+1]; nnk++)
	  esum += 2.0*dot(ElCoreDipole(iat,nnj),ElCoreDipole(iat,nnk));
      }
    
      //loop over ions and electrons 
      for(int jat=iat+1; jat<nCenters; jat++) {
	int nni = d_ie->M[iat];
	int nnj = d_ie->M[jat];
	for(int k=0; k<nParticles; k++, nni++, nnj++){
	  esum -= 2.0*dot(CoreCoreDipole(iat,jat),ElCoreDipole(iat,nni));
	  esum -= 2.0*dot(CoreCoreDipole(jat,iat),ElCoreDipole(jat,nnj));
	}
      }
    }//iat
    return C*esum + eCoreCore;
  }
}

  /***************************************************************************
   * $RCSfile$   $Author$
   * $Revision$   $Date$
   * $Id$ 
   ***************************************************************************/

  
