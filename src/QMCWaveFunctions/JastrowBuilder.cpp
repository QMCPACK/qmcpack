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
// -*- C++ -*-
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/JastrowBuilder.h"
#include "QMCWaveFunctions/PadeJastrow.h"
#include "QMCWaveFunctions/NoCuspJastrow.h"
#include "QMCWaveFunctions/OneBodyJastrowFunction.h"
#include "QMCWaveFunctions/TwoBodyJastrowFunction.h"

namespace ohmmsqmc {

  /** constructor
   *@param a the wavefunction
   *@param psets a vector containing the particle
   * 
   *Jastrow wavefunctions chose the needed distance tables and the 
   *DistanceTableData objects are initialized based on the source 
   *and target particle sets.
   */
  JastrowBuilder::JastrowBuilder(TrialWaveFunction& a, 
				 vector<ParticleSet*>& psets): 
    OrbitalBuilderBase(a), PtclSets(psets),
    corr_tag("correlation") 
  { 

  }
  
  /** Create a two-body spin-dependent Jatrow function with a template class JeeType.
   *@param cur the current xmlNode 
   *@param J2 a two-body jastrow function to be created.
   *
   *The template class JeeType is a functor which handles the
   *evaluation of the function value, gradients and laplacians using
   *distance tables. This is a specialized builder function for
   *spin-dependent Jastrow function,e.g., for electrons, two functions
   *are created for uu(dd) and ud(du).
   */
  template<class JeeType>
  bool 
  JastrowBuilder::createTwoBodySpin(xmlNodePtr cur, JeeType* J2) {
    
    /**\typedef The type of a simple function,e.g., PadeJastrow<double> */
    typedef typename JeeType::FuncType FuncType;

    vector<FuncType*> jastrow;
    DistanceTableData* d_table = NULL;
    cur = cur->xmlChildrenNode;
    int ng = 0, nj = 0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
      	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
	ParticleSet* a = NULL;
	int iptcl = 0;
	//check to see that the source particle set already exists
	while(!a && iptcl<PtclSets.size()) {
	  if(PtclSets[iptcl]->getName() == source_name) { a = PtclSets[iptcl];}
	  iptcl++;
	}
      	d_table = DistanceTable::getTable(DistanceTable::add(*a));
      	ng = a->groups();
	//create a Jastrow function for each pair type
	//for spin 1/2 particles (up-up, down-up = up-down,
	//down-down) 
      	for(int i=0; i<ng*ng; i++) jastrow.push_back(NULL);
      } else if(cname ==corr_tag) {
	string spA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
	string spB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));

	int ia = d_table->origin().Species.findSpecies(spA);
	int ib = d_table->origin().Species.findSpecies(spB);
	int iab = ia*ng+ib;
	if(!(jastrow[iab])) {
	  //create the new Jastrow function
	  FuncType *j2 = new FuncType;
	  //initialize
	  j2->put(cur,wfs_ref.VarList);
	  jastrow[iab]= j2;
	  if(ia != ib) {//up-down-type pair, treat down-up the same way
	    jastrow[ib*ng+ia] = j2;
	  } else {
	    for(int iaa=0; iaa<ng; iaa++) if(iaa != ia) jastrow[iaa*ng+iaa] = j2;
	  }
	  XMLReport("Added Jastrow Correlation between "<<spA<<" and "<<spB)
	  nj++;
	} else {
	  ERRORMSG("Using an existing Jastrow Correlation "<<spA<<" and "<<spB)
	}
      }
      cur = cur->next;
    } // while cur
    
    if(!nj) {
      ERRORMSG("No valid Jastrow Correlation is provided. Do Nothing")
      return false;
    }

    J2 = new JeeType(d_table);
    int ij = 0;    
    for(int i=0; i<ng; i++) {
      for(int j=0; j<ng; j++, ij++) {
	J2->F.push_back(jastrow[ij]);
      }
    }

    wfs_ref.add(J2);
    XMLReport("Added a Two-Body Jastrow Function")
    return true;
  }


  /** Create a two-body spin-independent Jatrow function with a template class JeeType.
   *@param cur the current xmlNode 
   *@param J2 a two-body jastrow function to be created.
   *
   *The template class JeeType is a functor which handles the
   *evaluation of the function value, gradients and laplacians using
   *distance tables. This is a specialized builder function for
   *spin-dependent Jastrow function,e.g., for electrons, two functions
   *are created for uu(dd) and ud(du).
   */
  template<class JeeType>
  bool 
  JastrowBuilder::createTwoBodyNoSpin(xmlNodePtr cur, JeeType* J2) {

    /**\typedef The type of a simple function,e.g., PadeJastrow<double> */
    typedef typename JeeType::FuncType FuncType;

    DistanceTableData* d_table = NULL;
    cur = cur->xmlChildrenNode;
    int ng = 0, nj = 0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
      	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
	ParticleSet* a = NULL;
	int iptcl = 0;
	//check to see that the source particle set already exists
	while(!a && iptcl<PtclSets.size()) {
	  if(PtclSets[iptcl]->getName() == source_name) { a = PtclSets[iptcl];}
	  iptcl++;
	}
      	d_table = DistanceTable::getTable(DistanceTable::add(*a));
	if(!J2) {
	  J2 = new JeeType(d_table);
	}
      } else if(cname == corr_tag) {
	if(J2) {
	  J2->F.put(cur,wfs_ref.VarList); nj++;
	} else {
	  ERRORMSG("No distance table is given for two-body jastrow function")
	}
      }
      cur = cur->next;
    } // while cur
    
    if(!nj) {
      ERRORMSG("No valid Jastrow Correlation is provided. Do Nothing")
      return false;
    }

    if(J2) {
      wfs_ref.add(J2);
      XMLReport("Added a Two-Body Jastrow Function")
      return true;
    } else {
      XMLReport("Failed to add a Two-Body Jastrow Function")
      return false;
    }
  }

  /** Create an one-body Jatrow function with a template class JneType.
   *@param cur the current xmlNode 
   *@param J1 an one-body jastrow function to be created.
   *
   *The template class JneType is a functor which handles the
   *evaluation of the function value, gradients and laplacians using
   *distance tables. This is a specialized builder function for
   *one-body Jastrow function, e.g.,nuclei-electron Jastrow function.
   */
  template<class JneType>
  bool 
  JastrowBuilder::createOneBody(xmlNodePtr cur, JneType* J1){

    typedef typename JneType::FuncType FuncType;

    vector<FuncType*> jastrow;
    DistanceTableData* d_table = NULL;
    cur = cur->xmlChildrenNode;
    int ng = 0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
      	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
      	string target_name((const char*)(xmlGetProp(cur,(const xmlChar *)"target")));
	ParticleSet* a = NULL;
	ParticleSet* b = NULL;
	int iptcl = 0;

	while((!a || !b) && iptcl<PtclSets.size()) {
	  if(PtclSets[iptcl]->getName() == source_name) { 
	    a = PtclSets[iptcl];
	  }
	  if(PtclSets[iptcl]->getName() == target_name) { 
	    b = PtclSets[iptcl];
	  }
	  iptcl++;
	}
	d_table = DistanceTable::getTable(DistanceTable::add(*a,*b));
	ng = a->Species.getTotalNum();
	XMLReport("Number of sources " << ng)
	for(int i=0; i<ng; i++) jastrow.push_back(NULL);
      } else if(cname == corr_tag) {
	string speciesA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
	//string speciesB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
	int ia = d_table->origin().Species.findSpecies(speciesA);
	if(!(jastrow[ia])) {
	  jastrow[ia]= new FuncType;
	  jastrow[ia]->put(cur,wfs_ref.VarList);
	  XMLReport("Added Jastrow Correlation between " << speciesA)
	}
      }
      cur = cur->next;
    } // while cur
    
    J1 = new JneType(d_table);
    FuncType* func_shared = NULL;
    int ig=0;
    while(!func_shared && ig<ng) {
      if(jastrow[ig]) func_shared = jastrow[ig];
      ig++;
    }
    if(!func_shared) {
      WARNMSG("Cannot create one-body Jastrow Function. Do nothing")
      return false;
    }

    ig = 0;    
    for(; ig<ng; ig++) {
      if(jastrow[ig]) {//(i,*) center is missing. Use the first valid jastrow function
	J1->F.push_back(jastrow[ig]);
      } else {
	J1->F.push_back(func_shared);
      }
    }

    wfs_ref.add(J1);
    XMLReport("Added a One-Body Jastrow Function")
    return true;
  }

  /** Create a Jastrow function and add it to a TrialWaveFunction.
   *@param cur the current xmlNode
   *@return true if successful
   *
   *This function is called whenever a jastrow element is found.
   *The default Jastrow function is the Pade, but there is also currently a NoCusp Jastrow.
   */
  bool JastrowBuilder::put(xmlNodePtr cur) {

    string jasttype((const char*)(xmlGetProp(cur, (const xmlChar *)"type")));
    string jastname((const char*)(xmlGetProp(cur, (const xmlChar *)"name")));
    string jastfunction((const char*)(xmlGetProp(cur, (const xmlChar *)"function")));

    LOGMSG("Jastrow Factor: Name = "<< jastname <<" Type = "<<jasttype)

      /*Currently for Two-Body Jastrow only Pade function: ar/(1+br) form 
	is implemted
      */
    if(jasttype == "Two-Body-Spin") {
      TwoBodyJastrow<PadeJastrow<ValueType>,false> *J2 = NULL;
      return createTwoBodySpin(cur,J2);
    } else if(jasttype == "Two-Body") {
      TwoBodyJastrow<PadeJastrow<ValueType>,true> *J2 = NULL;
      return createTwoBodyNoSpin(cur,J2);
    }
    /*For One-Body Jastrow default is Pade function, also implements
      the nocusp function a/(1+br^2)
    */
      else if(jasttype == "One-Body") {
      if(jastfunction == "nocusp") {
	LOGMSG("Jastrow Function: nucusp.")
	OneBodyJastrow<NoCuspJastrow<ValueType> > *J1 = NULL;
	return createOneBody(cur,J1);
      } else {
	LOGMSG("Jastrow Function: pade.")
	OneBodyJastrow<PadeJastrow<ValueType> > *J1 = NULL;
	return createOneBody(cur,J1);
      }
    }
    return false;
  } 
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
