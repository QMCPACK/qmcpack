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
#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
#include "QMCWaveFunctions/DetSetBuilderWithBasisSet.h"

namespace ohmmsqmc {

  STOMolecularOrbitals::STOMolecularOrbitals(TrialWaveFunction& wfs, 
					     ParticleSet& ions, 
					     ParticleSet& els):
    OrbitalBuilderBase(wfs),  BasisSet(NULL)
  { 
    int d_ie = DistanceTable::add(ions,els);
    d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
  }   
  
  STOMolecularOrbitals::BasisSetType* 
  STOMolecularOrbitals::addBasisSet(xmlNodePtr cur) {

    if(!BasisSet) BasisSet = new BasisSetType;

    int ncenters = CenterID.size();
    int activeCenter;

    //go thru the tree
    cur = cur->xmlChildrenNode;
    while(cur!=NULL) {
      string cname((const char*)(cur->name));
      if(cname == basis_tag) {
        xmlChar* att=xmlGetProp(cur, (const xmlChar*)"species");
        if(!att) {
	  ERRORMSG("//BasisSet/Basis does not have name. Failed to initialize.")         
	  return NULL;
        }     
	string abasis((const char*)att);

	//check, if the named center exists
	map<string,int>::iterator it = CenterID.find(abasis);

	if(it == CenterID.end()) {//add the name to the map CenterID
	  CenterID[abasis] = activeCenter = ncenters++;
	  int Lmax = 0, num=0;
	  xmlNodePtr cur1 = cur->xmlChildrenNode;
          vector<xmlNodePtr> basisfunc_ptr;
	  while(cur1 != NULL) {
            string cname1((const char*)(cur1->name));
	    if(cname1 == basisfunc_tag) {
              basisfunc_ptr.push_back(cur1);
    	      int l=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"l")));
	      Lmax = max(Lmax,l);num++;
	    }
	    cur1 = cur1->next;
	  }

	  XMLReport("Adding a center " << abasis 
		    << " with centerid " << CenterID[abasis])
          XMLReport("Maximum angular momentum    = " << Lmax)
          XMLReport("Number of centered orbitals = " << num)

	  //create a new set of orbitals sharing a center
	  CenteredOrbitalType* aos = new CenteredOrbitalType(Lmax);
	  aos->LM.resize(num);
	  aos->NL.resize(num);
	  num=0;
          for(int ib=0; ib<basisfunc_ptr.size(); ib++) {
            cur1=basisfunc_ptr[ib];
    	    int n=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"n")));
    	    int l=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"l")));
    	    int m=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"m")));
	    string rnl((const char*)(xmlGetProp(cur1, (const xmlChar *)"id")));

    	    //assign the index for (l,m)
	    aos->LM[num] = aos->Ylm.index(l,m);
	    double zeta = 1.0;
	    xmlNodePtr s = cur1->xmlChildrenNode;
	    while(s != NULL) {
              string cname((const char*)(s->name));
              if(cname == param_tag) {
	        putContent(zeta,s);
		break;
	      }
	      s=s->next;
	    }

	    XMLReport("A spherical orbital with (n,l,m,zeta) = (" << n << ", " << l << ", " << m << ", " << zeta << ")")

	    STONorm<RealType> anorm(n);
	    //radial orbitals: add only distinct orbitals
	    map<string,int>::iterator rnl_it = RnlID.find(rnl);
	    if(rnl_it == RnlID.end()) {
	      int nl = aos->Rnl.size();
	      aos->Rnl.push_back(new RadialOrbitalType(n-l-1,zeta,anorm(n-1,zeta)));
	      aos->NL[num] = nl;
	      RnlID[rnl] = nl;
	    } else {
	      aos->NL[num] = (*rnl_it).second;
	    }
	    num++;
	  }

	  BasisSet->add(aos);

	} else {
	  WARNMSG("Species " << 
		  abasis << " is already initialized. Ignore the input.")
	}
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

  bool STOMolecularOrbitals::put(xmlNodePtr cur){
    LOGMSG("STOMolecularOrbitals::put")
    DetSetBuilderWithBasisSet<STOMolecularOrbitals> spobuilder(wfs_ref,*this);
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
