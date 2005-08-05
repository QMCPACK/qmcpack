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
#include "QMCHamiltonians/GeCorePolPotential.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"
#include <fstream>
#


namespace ohmmsqmc {
 
  GeCorePolPotential::GeCorePolPotential(ParticleSet& ions, ParticleSet& els): 
    d_ie(NULL), d_ii(NULL), alpha(0.3558), r_b(0.7048), eCoreCore(0.0) { 
    
    //set the distance tables
    d_ie = DistanceTable::getTable(DistanceTable::add(ions,els));
    d_ie->create(1);
    d_ii = DistanceTable::getTable(DistanceTable::add(ions));
    d_ii->create(1);
    d_ii->evaluate(ions);

    nCenters = ions.getTotalNum();
    nParticles = els.getTotalNum();
    C = -0.5*alpha;
    r_binv = 1.0/r_b;
  
    CoreCoef.resize(nCenters);
    CoreCoreDipole.resize(nCenters,0.0);
    CoreElDipole.resize(nCenters,nParticles);

    //   CoreCoef = false;
    CoreElDipole = 0.0;
  
    int GeCounter = 0;
    SpeciesSet& Species(ions.getSpeciesSet());
    for(int iat=0; iat<nCenters; iat++){
      CoreCoef[iat] = false;
      string sname = Species.speciesName[ions.GroupID[iat]];
      if(sname == "Ge"){
	LOGMSG("Adding a core-electron potential for " << sname << " #" << GeCounter++)
	  CoreCoef[iat] = true;
      }
    }
    
  
    //index for attribute charge
    int iz = Species.addAttribute("charge");
    //calculate the Core-Core Dipole matrix
    for(int iat=0; iat<nCenters; iat++) {
      for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++) { 
	int jat(d_ii->J[nn]);
	//check to see if both ions are Ge
	RealType rinv3 = pow(d_ii->rinv(nn),3);//(1/R_{JI}^3) R_{JI} = R_J-R_I
	PosType dipole(rinv3*d_ii->dr(nn));//(\vec{R_{JI}}/R_{JI}^3)

        //Sign is here.
	CoreCoreDipole[iat] -= dipole*Species(iz,ions.GroupID[jat]);//charge of jat
	CoreCoreDipole[jat] += dipole*Species(iz,ions.GroupID[iat]);//charge of iat
      }
    }

    RealType corecore(0.0);
    for(int iat=0; iat<nCenters; iat++) {
      if(CoreCoef[iat]) corecore+=dot(CoreCoreDipole[iat],CoreCoreDipole[iat]);
    }
    LOGMSG("Core-Core Dipole = " << C*corecore);
  }
    
  GeCorePolPotential::~GeCorePolPotential() { }

  GeCorePolPotential::ValueType 
  GeCorePolPotential::evaluate(ParticleSet& P) {
    if(Primary) {
      //calculate the Electron-Core Dipole matrix
      //CoreElDipole=0.0;
      RealType e = 0.0;
      for(int iat=0; iat<nCenters; iat++) {
        if(CoreCoef[iat]){
          PosType cc(CoreCoreDipole[iat]);
          for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++){
            int eid(d_ie->J[nn]);
            RealType rinv3 = pow(d_ie->rinv(nn),3);//(1/r^3)
            PosType dipole = rinv3*d_ie->dr(nn);//(\vec{r}/r^3)
            cc +=  dipole*fcpp(d_ie->r(nn)*r_binv); //CoreElDipole(iat,eid) = dipole*fcpp(d_ie->r(nn)*r_binv);
          }
          //PosType cc(CoreCoreDipole[iat]+core_el);
          //for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++) { 
          //  cc += CoreElDipole(iat,d_ie->J[nn]);
          //}
          e += dot(cc,cc);
        }
      }
      Value=C*e;
    } 
    return Value;
  }

}

  /***************************************************************************
   * $RCSfile$   $Author$
   * $Revision$   $Date$
   * $Id$ 
   ***************************************************************************/
