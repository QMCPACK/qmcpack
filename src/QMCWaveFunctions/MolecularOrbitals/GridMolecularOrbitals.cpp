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
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/DetSetBuilderWithBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/RadialGridFunctorBuilder.h"

namespace ohmmsqmc {

  GridMolecularOrbitals::GridMolecularOrbitals(TrialWaveFunction& wfs, 
					       ParticleSet& ions, 
					       ParticleSet& els):
    OrbitalBuilderBase(wfs), 
    BasisSet(NULL)
  { 
    int d_ie = DistanceTable::add(ions,els);
    d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
  }   

  GridMolecularOrbitals::BasisSetType* 
  GridMolecularOrbitals::addBasisSet(xmlNodePtr cur) {

    if(!BasisSet) BasisSet = new BasisSetType;

    //map for quantum numbers
    map<string,int> nlms_id;
    nlms_id["n"] = q_n;
    nlms_id["l"] = q_l;
    nlms_id["m"] = q_m;
    nlms_id["s"] = q_s;

    QuantumNumberType nlms;
    string rnl;

    //current number of centers
    int ncenters = CenterID.size();
    int activeCenter;
    int gridmode = -1;

    //go thru the tree
    cur = cur->xmlChildrenNode;
    map<string,RGFBuilderBase*> rbuilderlist;

    while(cur!=NULL) {
      string cname((const char*)(cur->name));
      if(cname == basis_tag) {
        xmlChar* att=xmlGetProp(cur, (const xmlChar*)"species");
        if(!att) {
	  ERRORMSG("//BasisSet/Basis does not have name. Failed to initialize.")         
	  return NULL;
        }     

	string abasis((const char*)att);
	string btype((const char*)(xmlGetProp(cur, (const xmlChar *)"type")));

        bool expandlm = false;

	//choose between "STO" (Slater Type Orbitals) or the default
	//"RGOT" (Radial Grid Orbital Type)
        RGFBuilderBase* rbuilder = NULL; 
        if(btype == "STO") {
          rbuilder = new STO2GridBuilder;
        } else {
          rbuilder = new RGFBuilder;
        }

	//check, if the named center exists
	map<string,int>::iterator it = CenterID.find(abasis);

	if(it == CenterID.end()) {//add the name to the map CenterID
	  CenterID[abasis] = activeCenter = ncenters++;
	  int Lmax = 0, num=0;
	  xmlNodePtr cur1 = cur->xmlChildrenNode;
	  //find maximun angular momentum
	  while(cur1 != NULL) {
	    string cname1((const char*)(cur1->name));
	    if(cname1 == basisfunc_tag) {
    	      int l=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"l")));
	      Lmax = max(Lmax,l);
	      //expect that only Rnl is given
	      if(expandlm) 
		num += 2*l+1;
	      else		
		num++;
	    }
	    cur1 = cur1->next;
	  }

	  XMLReport("Adding a center " << abasis << " centerid " 
		    << CenterID[abasis])
          XMLReport("Maximum angular momentum    = " << Lmax)
          XMLReport("Number of centered orbitals = " << num)

	  //create a new set of atomic orbitals sharing a center
	  CenteredOrbitalType* aos = new CenteredOrbitalType(Lmax);
	  aos->LM.resize(num);
	  aos->NL.resize(num);
	  cur1 = cur->xmlChildrenNode;
	  num=0;
	  //assign radial orbitals for the new center
	  rbuilder->setOrbitalSet(aos,abasis);
	  //assign a radial grid for the new center
          rbuilder->addGrid(cur);
	  while(cur1 != NULL) {
	    string cname1((const char*)(cur1->name));
	    if(cname1 == basisfunc_tag) {
	      //check the attributes: n, l, m, s
	      xmlAttrPtr att = cur1->properties;
	      while(att != NULL) {
		string aname((const char*)(att->name));
		if(aname == "id") {
		  rnl = (const char*)(att->children->content);
		} else {
		  map<string,int>::iterator iit = nlms_id.find(aname);
		  if(iit != nlms_id.end()) {
		    nlms[(*iit).second] = atoi((const char*)(att->children->content));
		  } 
		}
		att = att->next;
	      }

	      XMLReport("")
	      XMLReport("A spherical orbital (n,l,m,s) " << nlms[0] 
			<< " " << nlms[1] << " " << nlms[2] << " " << nlms[3])
	      if(expandlm) {//add -l ..l
		map<string,int>::iterator rnl_it = RnlID.find(rnl);
		if(rnl_it == RnlID.end()) {
		  int nl = aos->Rnl.size();
		  if(rbuilder->addRadialOrbital(cur1,nlms)) {
		    RnlID[rnl] = nl;
		    int l = nlms[q_l];
		    XMLReport("Adding " << 2*l+1 << " spherical orbitals")
		    for(int tm=-l;tm<=l; tm++, num++) {
		      aos->LM[num] = aos->Ylm.index(l,tm);
		      aos->NL[num] = nl;
		    }
		  }
		}
	      } else {
		//assign the index for real Spherical Harmonic with (l,m)
		aos->LM[num] = aos->Ylm.index(nlms[q_l],nlms[q_m]);
		//radial orbitals: add only distinct orbitals
		map<string,int>::iterator rnl_it = RnlID.find(rnl);
		if(rnl_it == RnlID.end()) {
		  int nl = aos->Rnl.size();
		  if(rbuilder->addRadialOrbital(cur1,nlms)) {
		    //assign the index for radial orbital with (n,l)
		    aos->NL[num] = nl;
		    RnlID[rnl] = nl;
		  }
		} else {
		  //assign the index for radial orbital with (n,l) if repeated
		  XMLReport("Already added radial function for id: " << rnl)
		  aos->NL[num] = (*rnl_it).second;

		}
		//increment number of basis functions
		num++;
	      }
	    }
	    cur1 = cur1->next;
	  }
	  //add the new atomic basis to the basis set
	  BasisSet->add(aos);

	}else {
	  WARNMSG("Species " << abasis << " is already initialized. Ignore the input.")
	}
        if(rbuilder) delete rbuilder;
      }
      cur = cur->next;
    }

    if(BasisSet) {
      BasisSet->setTable(d_table);
      XMLReport("The total number of basis functions " << BasisSet->TotalBasis)
      return BasisSet;
    } else {
      ERRORMSG("BasisSet is not initialized.")
      return NULL;
    }
  }

  bool GridMolecularOrbitals::put(xmlNodePtr cur){
    LOGMSG("GridMolecularOrbitals::put")
    DetSetBuilderWithBasisSet<GridMolecularOrbitals> spobuilder(wfs_ref,*this);
    if(spobuilder.put(cur)) {
      BasisSet->resize(spobuilder.NumPtcl);
      return true;
    } else {
      return false;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
