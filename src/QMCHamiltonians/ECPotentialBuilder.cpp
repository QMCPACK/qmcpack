//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {
  /** constructor
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\param psi trial wavefunction
   */
  ECPotentialBuilder::ECPotentialBuilder(QMCHamiltonian& h,
      ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi): 
    hasLocalPot(false),hasNonLocalPot(false),
    targetH(h), IonConfig(ions), targetPtcl(els), targetPsi(psi)
  { }

  bool ECPotentialBuilder::put(xmlNodePtr cur) {

    if(localPot.empty()) {
      int ng(IonConfig.getSpeciesSet().getTotalNum());
      localZeff.resize(ng,1);
      localPot.resize(ng,0);
      nonLocalPot.resize(ng,0);
    }

    string ecpFormat("table");
    const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"format");
    if(t != NULL) {
      ecpFormat= (const char*)t;
    } 
    if(ecpFormat == "xml")  {
      useXmlFormat(cur);
    } else {
      useSimpleTableFormat();
    }

    ///create LocalECPotential
    if(hasLocalPot) {
      LocalECPotential* apot = new LocalECPotential(IonConfig,targetPtcl);
      for(int i=0; i<localPot.size(); i++) {
        if(localPot[i]) apot->add(i,localPot[i],localZeff[i]);
      }
      targetH.addOperator(apot,"LocalECP");
    }

    if(hasNonLocalPot) {
      //resize the sphere
      targetPtcl.resizeSphere(IonConfig.getTotalNum());

      NonLocalECPotential* apot = new NonLocalECPotential(IonConfig,targetPtcl,targetPsi);
      for(int i=0; i<nonLocalPot.size(); i++) {
        if(nonLocalPot[i]) apot->add(i,nonLocalPot[i]);
      }
      targetH.addOperator(apot,"NonLocalECP");

      for(int ic=0; ic<IonConfig.getTotalNum(); ic++) {
        int ig=IonConfig.GroupID[ic];
        if(nonLocalPot[ig]) { 
          if(nonLocalPot[ig]->nknot) targetPtcl.Sphere[ic]->resize(nonLocalPot[ig]->nknot);
        }
      }
    }
    return true;
  }

  void ECPotentialBuilder::useXmlFormat(xmlNodePtr cur) {

    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "pseudo") {
        string href("none");
        string ionName("none");
        OhmmsAttributeSet hAttrib;
        hAttrib.add(href,"href");
        hAttrib.add(ionName,"elementType"); hAttrib.add(ionName,"symbol");
        hAttrib.put(cur);

        int speciesIndex=IonConfig.getSpeciesSet().findSpecies(ionName);
        if(speciesIndex < IonConfig.getSpeciesSet().getTotalNum()) {
          app_log() << "  Adding pseudopotential for " << ionName << endl;
          ECPComponentBuilder ecp(ionName);
          bool success=false;
          if(href == "none") {
            success=ecp.put(cur);
          } else {
            success=ecp.parse(href);
          }
          if(success) {
            if(ecp.pp_loc) {
              localPot[speciesIndex]=ecp.pp_loc;
              localZeff[speciesIndex]=ecp.Zeff;
              hasLocalPot=true;
            }
            if(ecp.pp_nonloc) {
              nonLocalPot[speciesIndex]=ecp.pp_nonloc;
              hasNonLocalPot=true;
            }
          }
        } else {
          app_error() << "  Ion species " << ionName << " is not found." << endl;
        }
      } 
      cur=cur->next;
    }
  }

  /** reimplement simple table format used by NonLocalPPotential
   */
  void ECPotentialBuilder::useSimpleTableFormat() {

    const SpeciesSet& Species(IonConfig.getSpeciesSet());
    int ng(Species.getTotalNum());
    for(int ig=0; ig<ng;ig++) {
      vector<RealType> grid_temp, pp_temp;
      string species(Species.speciesName[ig]);
      string fname = species+".psf";
      ifstream fin(fname.c_str(),ios_base::in);
      if(!fin){
	ERRORMSG("Could not open file " << fname)
        exit(-1);
      }      

      // Read Number of potentials (local and non) for this atom
      int npotentials;
      fin >> npotentials;
      RealType r, f1;

      int lmax=-1;
      int numnonloc=0;
      RealType rmax(0.0);

      app_log() << "  ECPotential for " << species << endl;
      NonLocalECPComponent* mynnloc=0;

      for (int ij=0; ij<npotentials; ij++){
	int angmom,npoints;
	fin >> angmom >> npoints;
        if(grid_temp.size()<npoints) grid_temp.resize(npoints);
        if(pp_temp.size()<npoints) pp_temp.resize(npoints);
	for (int j=0; j<npoints; j++){
          fin >> grid_temp[j] >> pp_temp[j];
	}

        GridType *agrid = new NumericalGrid<ValueType>(grid_temp);
	RadialPotentialType *app = new RadialPotentialType(agrid,pp_temp);

	int imin = 0;
	RealType yprime_i = ((*app)(imin+1)-(*app)(imin))/app->dr(imin);
	if(angmom < 0) {
          hasLocalPot=true; //will create LocalECPotential
	  app->spline(imin,yprime_i,app->size()-1,0.0);
          localPot[ig]=app;
          app_log() << "    LocalECP l=" << angmom << " deriv= " << yprime_i << endl;
        } else {
          app_log() << "    NonLocalECP l=" << angmom << " rmax = " << agrid->rmax() << endl;
          hasNonLocalPot=true; //will create NonLocalECPotential
          if(mynnloc == 0) mynnloc = new NonLocalECPComponent;
	  app->spline(imin,yprime_i,app->size()-1,0.0);
	  mynnloc->add(angmom,app);
	  lmax=std::max(lmax,angmom);
	  rmax=std::max(rmax,agrid->rmax());
          numnonloc++;
	}

        if(mynnloc) {
          mynnloc->lmax=lmax; 
          mynnloc->Rmax=rmax;
          app_log() << "    Maximum cutoff of NonLocalECP " << rmax << endl;
        }
      } 
      fin.close();

      if(mynnloc) {
        nonLocalPot[ig]=mynnloc;
        int numsgridpts=0;

        string fname = species+".sgr";
        ifstream fin(fname.c_str(),ios_base::in);
        if(!fin){
          app_error() << "Could not open file " << fname << endl;
          exit(-1);
        }
        PosType xyz;
        ValueType weight;
        while(fin >> xyz >> weight){
          mynnloc->addknot(xyz,weight);
          numsgridpts++;
        }
        //cout << "Spherical grid : " << numsgridpts << " points" <<endl;
        mynnloc->resize_warrays(numsgridpts,numnonloc,lmax);
      }
    }//species
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
