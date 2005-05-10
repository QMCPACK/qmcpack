//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "QMCWaveFunctions/NumericalOrbitalSetBuilder.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "QMC/MO2Grid3D.h"

namespace ohmmsqmc {

  NumericalOrbitalSetBuilder::NumericalOrbitalSetBuilder(TrialWaveFunction& wfs): 
    OrbitalBuilderBase(wfs),GridXYZ(0) {
  }   

  bool NumericalOrbitalSetBuilder::put(xmlNodePtr cur){

    MO2Grid3D converter;

    //save current node
    xmlNodePtr curSave=cur;

    //MO2Grid3D returns a modified xml node
    xmlNodePtr curMod = converter.generateNumericalOrbitals(cur);

    if(curMod) {
      //copy the orbitals that have been already created
      converter.copyOrbitalSet(SPOSet);
      curSave=curMod;
    }

    if(SPOSet.empty()) initOrbitalSet(curSave);
    return addSlaterDeterminantSet(curSave);
  }

  void NumericalOrbitalSetBuilder::initOrbitalSet(xmlNodePtr cur){

    std::vector<RealType> ri(3,-5.0);
    std::vector<RealType> rf(3,5.0);
    std::vector<int> npts(3,101);
    xmlNodePtr curSave=cur;
    cur = cur->xmlChildrenNode;
    int idir=0;

    //first pass to check if basisset is used
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "spline") {
        splinePtr=cur; // save spline data
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
          string tname((const char*)(tcur->name));
          if(tname == "grid") {
            a=xmlGetProp(tcur,(const xmlChar*)"dir");
            if(a) { idir=atoi((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"ri");
            if(a) { ri[idir]=atof((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"rf");
            if(a) { rf[idir]=atof((const char*)a);}
            a=xmlGetProp(tcur,(const xmlChar*)"npts");
            if(a) { npts[idir]=atoi((const char*)a);}
          }
          tcur = tcur->next;
        }
      } 
      cur = cur->next;
    }

    typedef LinearGrid<RealType> GridType;
    GridType *gridX=new GridType;
    GridType *gridY=new GridType;
    GridType *gridZ=new GridType;

    gridX->set(ri[0],rf[0],npts[0]);
    gridY->set(ri[1],rf[1],npts[1]);
    gridZ->set(ri[2],rf[2],npts[2]);
    GridXYZ = new XYZCubicGrid<RealType>(gridX,gridY,gridZ);

  }

  bool NumericalOrbitalSetBuilder::addSlaterDeterminantSet(xmlNodePtr cur){

    typedef DiracDeterminant<SPOSetType>  Det_t;
    typedef SlaterDeterminant<SPOSetType> SlaterDeterminant_t;
    std::vector<SlaterDeterminant_t*> slaterdets;
    int is=0;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == OrbitalBuilderBase::sd_tag) {
	slaterdets.push_back(new SlaterDeterminant_t);
	sdet_coeff.push_back(1.0);
	first = 0;
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
	  string tname((const char*)(tcur->name));
	  if(tname == OrbitalBuilderBase::det_tag) {
            SPOSetType* psi=createSPOSet(tcur);
	    Det_t *adet = new Det_t(*psi,first);
	    adet->set(first,psi->size());
	    slaterdets[is]->add(adet);
            first+=psi->size();
	  }
	  tcur = tcur->next;
	}
      }
      cur = cur->next;
    }

    if(slaterdets.size() > 1) {
      XMLReport("Creating a multi-determinant wavefunction")
      MultiSlaterDeterminant<SPOSetType>
        *multidet= new MultiSlaterDeterminant<SPOSetType>;
      for(int i=0; i<slaterdets.size(); i++) {
        XMLReport("Coefficient for a SlaterDeterminant " << sdet_coeff[i])
        slaterdets[i]->setOptimizable(optimizeit);
        multidet->add(slaterdets[i],sdet_coeff[i]);
      }
      multidet->setOptimizable(optimizeit);
      //add a MultiDeterminant to the trial wavefuntion
      wfs_ref.add(multidet);
    } else {
      XMLReport("Creating a SlaterDeterminant wavefunction")
      wfs_ref.add(slaterdets[0]);
    }
  }

  NumericalOrbitalSetBuilder::SPOSetType*
  NumericalOrbitalSetBuilder::createSPOSet(xmlNodePtr cur) {

    //get id attribute
    string detname("det");
    a=xmlGetProp(cur,(const xmlChar*)"id");
    if(a) { detname = (const char*)a; }

    //get ref attribute
    string detref(detname);
    a=xmlGetProp(cur,(const xmlChar*)"ref");
    if(a) { detref=(const char*)a; }

    //add src attribute
    string detsrc;
    a=xmlGetProp(cur,(const xmlChar*)"src");
    if(a) {detsrc=(const char*)a;}

    //create a new Single-Particle-Orbital Set that will be connected to a DiracDeterminant
    SPOSetType* psi=new SPOSetType();

    //count the number of orbitals of this DiracDeterminant
    a=xmlGetProp(cur,(const xmlChar*)"orbitals");
    int norbs(atoi((const char*)a));

    //first, add the psi's explicitly added
    cur=cur->children;
    while(cur != NULL && norbs) {
      string cname((const char*)cur->name);
      if(cname == "psi") {
        a=xmlGetProp(cur,(const xmlChar*)"index");
        if(a) { 
          SPOType* orb = createSPO(detsrc,atoi(const char*)a);
          if(orb) { psi->add(orb); --norbs;}
        }
      }
      cur=cur->next;
    }

    //second, add the rest from the bottom of the eigen states
    int iorb(0);
    while(norbs) {
      SPOType* orb = createSPO(detsrc,iorb);
      if(orb) { psi->add(orb); --norbs; ++iorb;}
    }

    return psi;
  }

  NumericalOrbitalSetBuilder::SPOType*
  NumericalOrbitalSetBuilder::createSPO(const string& srcfile, int iorb) {

    char oname[128];
    sprintf(oname,"%s.wf%04d",srcfile.c_str(),iorb);
    SPOType* orb=0;
    map<string,TriCubicSplineT<ValueType>* >::iterator sit(SPOSet.find());
    if(sit != SPOSet.end()) {
      orb=(*sit).second;
    } else {
      //create one
    }
    return orb;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
