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
#include "Configuration.h"
#include "QMC/QDApps.h"
#include "QMCBase/RandomSeqGenerator.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
//#include "QMCHamiltonians/HartreePotential.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/HarmonicPotential.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/TriCubicSplinePotential.h"
//#include "QMCWaveFunctions/TriCubicSplineBuilder.h"
#include "QMCWaveFunctions/QDwfBuilder.h"
#include "QMCWaveFunctions/JastrowBuilder.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMC/QMCUtilities.h"


#include <iostream>

namespace qmcplusplus{

  using namespace xmlpp;
  
  /// Constructor
  QDApps::QDApps(int argc, char** argv): QMCApps(argc,argv), 
					 DeviceGrid(NULL) { 
    el.m_name = "e";
    int iu = el.Species.addSpecies("u");
    int id = el.Species.addSpecies("d");
    int icharge = el.Species.addAttribute("charge");
    el.Species(icharge,iu) = -1;
    el.Species(icharge,id) = -1;
  }

  /// Destructor
  QDApps::~QDApps() {
    DEBUGMSG("QDApps::~QDApps")
  }

  bool QDApps::init(){

    setParticleSets(d_root);

    setQDot(d_root);

    setRefParams();

    readGrid(d_root);

    setWavefunctions(d_root);

    setHamiltonian(d_root);

    if(!setMCWalkers(d_root)) {
      
      int nup = el.last(0);
      int nions = ion.getTotalNum();
      cout << "ulength is " << ulength << endl;
      double r=0.0;//1.2*ulength;
      for (int ipart=0; ipart<el.getTotalNum(); ++ipart) {
	double costheta=2*Random()-1;
	double sintheta=sqrt(1-costheta*costheta);
	double phi=2*M_PI*Random();
	PosType rpos = DeviceGrid->ptr(int(el.R[ipart][0]),
				       int(el.R[ipart][1]),
				       int(el.R[ipart][2]));
	el.R[ipart] =  rpos + PosType
	  (r*cos(phi)*sintheta, r*sin(phi)*sintheta,r*costheta);
	cout << "particle: " << el.R[ipart] << endl;
      }
      
    }

    /*
    //// -----------------------------
    DistanceTable::create(1);
    gridvec_t ir(79,19,48);

    for(int ix = 43; ix <= 43; ix++){
      ir[0] = ix; 
      elocal(ir);
    }
    exit(-1);
    //// -----------------------------
    */

    return true;
  }

  void QDApps::setRefParams(){
    double zref = el.R[0][2]; 
    for(int i = 0; i < inv_meff.size(); i++){
      if( zref > m_M[i] && zref < m_M[i+1] ){
	minv_ref = inv_meff[i]; eps_ref = eps_m[i];
      }
    }
    epsbym = eps_ref * minv_ref;
    mbyepsq = 1.0/(eps_ref*epsbym);

    return;
  }

  void QDApps::elocal(const gridvec_t& ir){

    int n = 0; double elocal;
    el.R[0] = DeviceGrid->ptr(ir);
    //    el.R[0][0] = 79.626;
    //    el.R[0][1] = 26.638;
    //    el.R[0][2] = 14.463;
    //el.R[0][0] = 8.0;
    for(int ii = 0; ii < 5000; ii++){
      double psi = Psi.evaluate(el);
      elocal = H.evaluate(el);
      cout << el.R[0][n] << '\t' << H[1] << '\t' << H[2] << '\t' << psi << '\t' << elocal << endl;
      el.R[0][n] += 0.012;
    }

    return;
  }

  bool QDApps::setQDot(xmlpp::Node* root){

    /// find "Structure" and point to the child node
    NodeSet qdset = root->find("//Structure");    

    xmlNodePtr cur = qdset[0]->cobj()->xmlChildrenNode;

    int zf;
    while( cur != NULL ){
      string cname((const char*)(cur->name));
      if( cname == "layer"){
	zf = atoi((char*)xmlGetProp(cur,(xmlChar*)"zf"));
	m_M.push_back(atoi((char*)xmlGetProp(cur,(xmlChar*)"zi")));
	mat_m.push_back(atoi((char*)xmlGetProp(cur,(xmlChar*)"mat")));
      }
      cur = cur->next;
    }
    m_M.push_back(zf);

    NodeSet qdset1 = root->find("//Materials");
    xmlNodePtr cur1 = qdset1[0]->cobj()->xmlChildrenNode;

    while( cur1 != NULL){
      string cname((const char*)(cur1->name));
      if( cname == "mat" ){
	double meff = atof((char*)xmlGetProp(cur1,(xmlChar*)"meff"));
	inv_meff.push_back(1.0/meff);
	e_gap.push_back(atof((char*)xmlGetProp(cur1,(xmlChar*)"Eg")));
	eps_m.push_back(atof((char*)xmlGetProp(cur1,(xmlChar*)"eps")));
	phi_s.push_back(atof((char*)xmlGetProp(cur1,(xmlChar*)"phis")));
	prior_m.push_back(atoi((char*)xmlGetProp(cur1,(xmlChar*)"prior")));
      }
      cur1 = cur1->next;
    }

    /// rearrange according to mat
    rearrange(prior_m);    rearrange(inv_meff);
    rearrange(e_gap);    rearrange(eps_m);    rearrange(phi_s);

    return true;

  }

  bool QDApps::setParticleSets(xmlpp::Node* root){

    bool init_els = determineNumOfElectrons(el,root);
    
    NodeSet pset = root->find("//ParticleSet");
    for(int i=0; i<pset.size(); i++) {
      ///first the name of the current ParticleSet
      Attribute* aname = dynamic_cast<Element*>(pset[i])->get_attribute("name");  
      char fc = aname->get_value()[0];
      if(fc == 'e') {
	if(init_els) {
	  el.m_name = "e";
	  XMLReport("The configuration for electrons is already determined by the wave function")
          XMLParticleIO pread(el,true);
          pread.put(pset[i]);
 	} else {
          XMLParticleIO pread(el);
          pread.put(pset[i]);
        }
      }
    }
    
    cout << "Electronic configuration : " << el.m_name << endl;
    for(int iat=0; iat<el.getTotalNum(); iat++) 
      cout << el.GroupID[iat] << " " << el.R[iat] << endl;


    return true; 


  }

  bool QDApps::setIons(xmlpp::Node* pnode) {
    XMLParticleIO pread(ion);
    pread.put(pnode);
    return true;
  }


  bool QDApps::rearrange(std::vector<double>& array){

    std::vector<double> temp;

    for(int i = 0; i < array.size(); i++) temp.push_back(array[i]);
    for(int i = 0; i < mat_m.size(); i++) array[i] = temp[mat_m[i]];

    return true;
  }

  bool QDApps::rearrange(std::vector<int>& array){

    std::vector<int> temp;

    for(int i = 0; i < array.size(); i++) temp.push_back(array[i]);
    for(int i = 0; i < mat_m.size(); i++) array[i] = temp[mat_m[i]];

    return true;

  }

  bool QDApps::setHamiltonian(xmlpp::Node* root){

    cout << "The name of potential file " << v_file << endl;
    ///H.add(new CoulombPotentialAA(el), "hartree");
    H.add(new BareKineticEnergy(),"KE");
    H.add(new TriCubicSplinePotential(mbyepsq,DeviceGrid,v_file),"Pot");
    return true;
  }

  bool QDApps::readGrid(xmlpp::Node* root){

    /// If the grid has not been created yet, create it
    if(!DeviceGrid)  {
       WARNMSG("Creating device grid")
       DeviceGrid = new Grid3D();
    }

    NodeSet pset = root->find("//HeteroStructure");
    if(pset.empty()){
      ERRORMSG("HeteroStructure is missing. Exiting ..." << endl)
      return false;
    } else {
      xmlNodePtr cur = pset[0]->cobj()->xmlChildrenNode;
      while ( cur != NULL ) {
	if(!xmlStrcmp( cur->name,(const xmlChar *) "Grid3D"))
	  ulength = DeviceGrid->put(cur,epsbym); 
	cur = cur->next;
      }
      LOGMSG("Completed the initialization of the Grid.")
    }
    return true;
  }


  bool QDApps::setWavefunctions(xmlpp::Node* root){

    NodeSet pset = root->find("//Wavefunction");
    vector<ParticleSet*> PtclSets;
    PtclSets.push_back(&ion);
    PtclSets.push_back(&el);
    if(pset.empty()) {
      ERRORMSG("Wavefunction is missing. Exit." << endl)
	return false;
    } else {
      if(pset.size() > 1) {
	WARNMSG("Mulptiple wavefunction is provided. Only the first wavefunction counts")
      }
      xmlNodePtr cur = pset[0]->cobj()->xmlChildrenNode;
      while(cur != NULL) {
	if ((!xmlStrcmp(cur->name, (const xmlChar *)"DeterminantSet"))) {
	  //TriCubicSplineBuilder s(Psi,DeviceGrid);
	  //s.put(cur);
	  QDwfBuilder QD(Psi);
	  QD.put(cur);
	  XMLReport("Done with the initialization of SlaterDeterminat.")
	 } // if DeterminantSet      
	if((!xmlStrcmp(cur->name, (const xmlChar *)"Jastrow"))) {
	  JastrowBuilder a(Psi,PtclSets);
	  a.put(cur);
	}
	cur = cur->next;
      }
      LOGMSG("Completed the initialization of a many-body wave function")
    }

    /// Now read in the file for the potential data
    NodeSet fset = root->find("//PotFile");
    v_file = string((char*)xmlGetProp(fset[0]->cobj(),(xmlChar*)"vfile"));



    return true;
  }





}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
