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
#include "QMCHamiltonians/LocalCorePolPotential.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/IteratorUtility.h"

namespace ohmmsqmc {
 
  LocalCorePolPotential::LocalCorePolPotential(ParticleSet& ions, 
					       ParticleSet& els): 
    FirstTime(true), eCoreCore(0.0), IonSys(ions), d_ie(0), d_ii(0) { 
    
    //set the distance tables
    d_ie = DistanceTable::getTable(DistanceTable::add(ions,els));
    d_ii = DistanceTable::getTable(DistanceTable::add(ions));

    nCenters = ions.getTotalNum();
    nParticles = els.getTotalNum();
  
    CoreCoreDipole.resize(nCenters,0.0);
    CoreElDipole.resize(nCenters,nParticles);
    CoreElDipole = 0.0;
  }
  
  /** destructor
   *
   * Delete InpCPP.
   */
  LocalCorePolPotential::~LocalCorePolPotential() { 
    for(int i=0; i<InpCPP.size(); i++) if(InpCPP[i]) delete InpCPP[i];
  }
  

  void LocalCorePolPotential::resetTargetParticleSet(ParticleSet& P) {
    d_ie = DistanceTable::getTable(DistanceTable::add(d_ie->origin(),P));
  }

  /** process xml node for each element
   * @param cur xmlnode <element name="string" alpha="double" rb="double"/> 
   */
  bool LocalCorePolPotential::CPP_Param::put(xmlNodePtr cur) {
    const xmlChar* a_ptr = xmlGetProp(cur,(const xmlChar *)"alpha");
    const xmlChar* b_ptr = xmlGetProp(cur,(const xmlChar *)"rb");
    if(a_ptr) alpha = atof((const char*)a_ptr);
    if(b_ptr) r_b = atof((const char*)b_ptr);
    C = -0.5*alpha;
    one_over_rr = 1.0/r_b/r_b;
    LOGMSG("\talpha = " << alpha << " rb = " << r_b)
    return true;
  }
  
  /** process xml node for CPP
   * @param cur xmlnode containing element+
   * 
   * element/@name is used to find the index of the element of the 
   * IonSys::SpeciesSet. The size of InpCPP is the number of species.
   * The size of Centers is the number of ions.
   */
  bool LocalCorePolPotential::put(xmlNodePtr cur){
    string ename;

    InpCPP.resize(IonSys.getSpeciesSet().getTotalNum(),0);
    cur= cur->children;
    bool success(false);
    while(cur != NULL){
      string cname((const char*)cur->name);
      if(cname == "element"){
      	const xmlChar* e_ptr = xmlGetProp(cur,(const xmlChar*)"name");
      	if(e_ptr){
          int itype = IonSys.getSpeciesSet().addSpecies((const char*)e_ptr);
          if(InpCPP[itype]==0) InpCPP[itype] = new CPP_Param;
          LOGMSG("CPP parameters for " << IonSys.getSpeciesSet().speciesName[itype])
          success &= InpCPP[itype]->put(cur);
      	}
      }
      cur=cur->next;
    }

    //resize Centers by the number of centers
    Centers.resize(nCenters,0);
    for(int iat=0; iat<nCenters; iat++) {
      Centers[iat]=InpCPP[IonSys.GroupID[iat]];
    }
    return success;
  }

  LocalCorePolPotential::Return_t
  LocalCorePolPotential::evaluate(ParticleSet& P) {

    if(FirstTime) {
      //index for attribute charge
      SpeciesSet& Species(IonSys.getSpeciesSet());
      int iz = Species.addAttribute("charge");
      //calculate the Core-Core Dipole matrix
      for(int iat=0; iat<nCenters; iat++) {
        for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++) { 
          int jat(d_ii->J[nn]);
          RealType rinv3 = pow(d_ii->rinv(nn),3);//(1/R_{JI}^3) R_{JI} = R_J-R_I
          PosType dipole(rinv3*d_ii->dr(nn));//(\vec{R_{JI}}/R_{JI}^3)
          //Sign and the charge of the paired ion are taken into account here 
          CoreCoreDipole[iat] -= dipole*Species(iz,IonSys.GroupID[jat]);
          CoreCoreDipole[jat] += dipole*Species(iz,IonSys.GroupID[iat]);
        }
      }

      RealType corecore(0.0);
      for(int iat=0; iat<nCenters; iat++) {
        if(Centers[iat]) 
          corecore+= Centers[iat]->C*dot(CoreCoreDipole[iat],CoreCoreDipole[iat]);
      }
      LOGMSG("Core-Core Dipole = " << corecore);
      FirstTime=false;
    }
    //calculate the Electron-Core Dipole matrix
    //CoreElDipole=0.0;

    RealType e = 0.0;
    for(int iat=0; iat<nCenters; iat++) {
      if(Centers[iat]) {
        PosType cc(CoreCoreDipole[iat]);
        for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++){
          int eid(d_ie->J[nn]);
          RealType rinv3 = pow(d_ie->rinv(nn),3);//(1/r^3)
          PosType dipole = rinv3*d_ie->dr(nn);//(\vec{r}/r^3)
          //cc +=  dipole*fcpp(d_ie->r(nn)*r_binv); 
          cc += dipole*((*Centers[iat])(d_ie->r(nn)));
        }
        e += Centers[iat]->C*dot(cc,cc);
      }
    }
    return Value=e;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
