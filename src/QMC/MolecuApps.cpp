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
#include "QMC/MolecuApps.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/LocalPPotential.h"
#include "QMCHamiltonians/HarmonicPotential.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.h"
#include "QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.h"
#include "QMCWaveFunctions/JastrowBuilder.h"
#include "QMC/QMCUtilities.h"

namespace ohmmsqmc {

  MolecuApps::MolecuApps(int argc, char** argv): QMCApps(argc,argv) { 
    el.setName("e");
    int iu = el.Species.addSpecies("u");
    int id = el.Species.addSpecies("d");
    int icharge = el.Species.addAttribute("charge");
    el.Species(icharge,iu) = -1;
    el.Species(icharge,id) = -1;
  }

  ///destructor
  MolecuApps::~MolecuApps() {
    DEBUGMSG("MolecuApps::~MolecuApps")
  }

  bool MolecuApps::init() {

    if(!setParticleSets(m_root)) {
       ERRORMSG("Failed to initialize the ions and electrons. Exit now.")
       return false;
    }

    setWavefunctions(m_root);

    if(!setHamiltonian(m_root)) {
      ERRORMSG("Failed to initialize the Hamitonians. Exit now.")
	return false;
    }
    
    if(!setMCWalkers(m_root)) {
      int nup = el.last(0);
      int nions = ion.getTotalNum();
      double r=0.012;
      for (int ipart=0; ipart<el.getTotalNum(); ++ipart) {
	double costheta=2*Random()-1;
	double sintheta=sqrt(1-costheta*costheta);
	double phi=2*3.141592653*Random();
	el.R[ipart] += 
	  MCWalkerConfiguration::PosType(r*cos(phi)*sintheta,
					 r*sin(phi)*sintheta,r*costheta);
      }
    }

    cout << "Ionic configuration : " << ion.getName() << endl;
    ion.get(cout);

    cout << "Electronic configuration : " << el.getName() << endl;
    el.get(cout);

    return true;    
  }   

  bool  MolecuApps::setParticleSets(xmlNodePtr aroot) {

    bool init_els = determineNumOfElectrons(el,m_context);

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//particleset",m_context);

    xmlNodePtr el_ptr=NULL, ion_ptr=NULL;
    for(int i=0; i<result->nodesetval->nodeNr; i++) {
      xmlNodePtr cur=result->nodesetval->nodeTab[i];
      xmlChar* aname= xmlGetProp(cur,(const xmlChar*)"name");
      if(aname) {
	char fc = aname[0];
	if(fc == 'e') { el_ptr=cur;}
	else if(fc == 'i') {ion_ptr=cur;}
      }
    }

    bool donotresize = false;
    if(init_els) {
      el.setName("e");
      XMLReport("The configuration for electrons is already determined by the wave function")
      donotresize = true;
    } 
    if(el_ptr) {
      XMLParticleParser pread(el,donotresize);
      pread.put(el_ptr);
    }

    if(ion_ptr) {
      XMLParticleParser pread(ion);
      pread.put(ion_ptr);
    }

    xmlXPathFreeObject(result);

    if(!ion.getTotalNum()) {
      ion.setName("i");
      ion.create(1);
      ion.R[0] = 0.0;
    }

    return true;
  }

  /**
   *@brief Initialize the Hamiltonian
   */
  bool MolecuApps::setHamiltonian(xmlNodePtr aroot){

    string ptype("molecule");
    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//hamiltonian",m_context);

    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      return false;
    } else {
      xmlChar* att= xmlGetProp(result->nodesetval->nodeTab[0],(const xmlChar*)"type");
      if(att) {
	ptype = (const char*)att;
      }
    }
    xmlXPathFreeObject(result);

    XMLReport("Found Potential of type " << ptype)
    //always add kinetic energy first
    H.add(new BareKineticEnergy, "Kinetic");
    if(ptype == "molecule"){
      H.add(new CoulombPotentialAA(el),"ElecElec");
      H.add(new CoulombPotentialAB(ion,el),"Coulomb");
      if(ion.getTotalNum()>1) 
	H.add(new IonIonPotential(ion),"IonIon");
    } else if(ptype == "harmonic") {
      H.add(new CoulombPotentialAA(el),"ElecElec");
      H.add(new HarmonicPotential(ion,el),"Coulomb");
    } else if(ptype == "siesta") {
      H.add(new CoulombPotentialAA(el), "ElecElec");
      H.add(new LocalPPotential(ion,el), "PseudoPot");
      if(ion.getTotalNum()>1) 
	H.add(new IonIonPotential(ion),"IonIon");
    } else if(ptype == "cpp") {
      H.add(new CoulombPotentialAA(el), "ElecElec");
      H.add(new LocalPPotential(ion,el), "PseudoPot");
      //     H.add(new CPPelel(ion,el), "CPPelel");
      if(ion.getTotalNum()>1) 
	H.add(new IonIonPotential(ion),"IonIon");
    } else {
      ERRORMSG(ptype << " is not supported.")
	return false;
    }


    return true;
  }
  

  /** Find a xmlnode Wavefunction and initialize the TrialWaveFunction by adding components.
   *@param root xml node
   *
   *This function substitutes TrialWaveFunctionBuilder
   *that is intended as a Builder for any TrialWaveFunction.
   *Since MolecuApps is specialized for molecular systems,
   *the type of allowed many-body wave functions can be decided here.
   *
   *Allowed many-body wave functions for MolecuApps are
   <ul>
   <li> DeterminantSet: only one kind of SlateDetermant can be used.
   <ul>
   <li> STO-He-Optimized: HePresetHF (specialized for He, test only)
   <li> STO-Clementi-Rottie: HFAtomicSTOSet (Hartree-Fock orbitals for atoms)
   <li> MolecularOrbital: molecular orbitals with radial orbitals 
   </ul>
   <li> Jastrow: any of them can be used
   <ul>
   <li> One-Body Jastrow
   <li> Two-Body Jastrow
   </ul>
   </ul>
   The number of terms is arbitrary.
   */
  bool MolecuApps::setWavefunctions(xmlNodePtr aroot) {

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//wavefunction",m_context);

    ///make a temporary array to pass over JastrowBuilder
    vector<ParticleSet*> PtclSets;
    PtclSets.push_back(&ion);
    PtclSets.push_back(&el);

    bool foundwfs=true;
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      ERRORMSG("Wavefunction is missing. Exit." << endl)
      foundwfs=false;
    } else {
      xmlNodePtr cur = result->nodesetval->nodeTab[0]->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if (cname == OrbitalBuilderBase::detset_tag) {
	  string orbtype=(const char*)(xmlGetProp(cur, (const xmlChar *)"type"));
	  LOGMSG("Slater-determinant terms using " << orbtype)
	  if(orbtype == "STO-Clementi-Rottie") {
            HFAtomicSTOSetBuilder a(Psi,ion,el);
	    a.put(cur);
	  } else if(orbtype == "STO-He-Optimized") {
	    HePresetHFBuilder a(Psi,ion,el);
	    //a.put(cur);
	  } else if(orbtype == "MolecularOrbital") {
	    MolecularOrbitalBuilder a(Psi,ion,el);
	    a.put(cur);
	  }
	  XMLReport("Done with the initialization of SlaterDeterminant using " << orbtype)
        } // if DeterminantSet    
	else if (cname ==  OrbitalBuilderBase::jastrow_tag) {
	  JastrowBuilder a(Psi,PtclSets);
	  a.put(cur);
	}
	cur = cur->next;
      }
      LOGMSG("Completed the initialization of a many-body wave function")
    }
    xmlXPathFreeObject(result);
    return foundwfs;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
