////////////////////////////////////////////////////////////////// // (c) Copyright 2003-  by Jeongnim Kim
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
#include "QMCWaveFunctions/PolarizedJastrow.h"

namespace ohmmsqmc {

  /** constructor
   * @param p target ParticleSet whose wave function is to be initialized
   *@param psi the wavefunction
   *@param psets a vector containing the ParticleSets
   * 
   *Jastrow wavefunctions chose the needed distance tables and the 
   *DistanceTableData objects are initialized based on the source 
   *and target particle sets.
   */
  JastrowBuilder::JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, 
     PtclPoolType& psets):OrbitalBuilderBase(p,psi),
                          ptclPool(psets), corr_tag("correlation") 
  { 

  }
  
  /** Create a two-body spin-dependent Jatrow function with a template 
   *class JeeType.
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

    int cur_var = targetPsi.VarList.size();

    DistanceTableData* d_table = DistanceTable::getTable(DistanceTable::add(targetPtcl));
    int	ng = targetPtcl.groups();

    map<string,FuncType*> jastrowMap;
    vector<FuncType*> jastrow(ng*ng);
    for(int i=0; i<ng*ng; i++) jastrow[i]=0;
    int nj = 0;
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      //if(cname == dtable_tag) {
      //	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
      //  //int iptcl = 0;
      //  map<string,ParticleSet*>::iterator pit(ptclPool.find(source_name));
      //  if(pit == ptclPool.end()) return false;
      //  ParticleSet* a = (*pit).second;
      //	d_table = DistanceTable::getTable(DistanceTable::add(*a));
      //	ng = a->groups();
      //  //create a Jastrow function for each pair type
      //  //for spin 1/2 particles (up-up, down-up = up-down,
      //  //down-down) 
      //	for(int i=0; i<ng*ng; i++) jastrow.push_back(NULL);
      //} else if(cname ==corr_tag) {
      if(cname ==corr_tag) {
	string spA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
	string spB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
        const xmlChar* refptr=xmlGetProp(cur,(const xmlChar *)"ref");
        const xmlChar* idptr=xmlGetProp(cur,(const xmlChar *)"id");
	int ia = targetPtcl.getSpeciesSet().findSpecies(spA);
	int ib = targetPtcl.getSpeciesSet().findSpecies(spB);
	int iab = ia*ng+ib;
	if(!(jastrow[iab])) {
	  //create the new Jastrow function
	  FuncType *j2=NULL;
          if(refptr == NULL) {
            j2 = new FuncType;
          } else {
            map<string,FuncType*>::iterator it(jastrowMap.find((const char*)refptr));
            if(it != jastrowMap.end()) {
              j2 = new FuncType((*it).second);
            } else { 
              j2 = new FuncType;
            }
          }
          if(idptr == NULL) {
            ostringstream idassigned; idassigned << "j2"<<iab;
            jastrowMap[idassigned.str()]=j2;
          } else {
            jastrowMap[(const char*)idptr]=j2;
          }

	  //initialize
	  j2->put(cur,targetPsi.VarList);
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

    J2 = new JeeType(targetPtcl,d_table);
    int ij = 0;    
    for(int i=0; i<ng; i++) {
      for(int j=0; j<ng; j++, ij++) {
	J2->F.push_back(jastrow[ij]);
      }
    }

    //set this jastrow function to be not optimizable
    if(targetPsi.VarList.size() == cur_var) {
      J2->setOptimizable(false);
    }

    targetPsi.add(J2);
    XMLReport("Added a Two-Body Jastrow Function")
    return true;
  }


  /** Create a two-body spin-independent Jatrow function with a template
   *class JeeType.
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

    int cur_var = targetPsi.VarList.size();

    DistanceTableData* d_table = DistanceTable::getTable(DistanceTable::add(targetPtcl));
    J2 = new JeeType(targetPtcl,d_table);
    int nj = 0;
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      //if(cname == dtable_tag) {
      //	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
      //  map<string,ParticleSet*>::iterator pit(ptclPool.find(source_name));
      //  if(pit == ptclPool.end()) return false;
      //  ParticleSet* a = (*pit).second;
      //	d_table = DistanceTable::getTable(DistanceTable::add(*a));
      //  if(!J2) {
      //    LOGMSG("Creating spin-independent two-body Jastrow " << a->getName())
      //    J2 = new JeeType(targetPtcl,d_table);
      //  }
      //} else if(cname == corr_tag) {
      if(cname == corr_tag) {
	if(J2) {
	  J2->F.put(cur,targetPsi.VarList); nj++;
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

    //set this jastrow function to be not optimizable
    if(targetPsi.VarList.size() == cur_var) {
      J2->setOptimizable(false);
    }
    targetPsi.add(J2);
    XMLReport("Added a Two-Body Jastrow Function")
    return true;
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
    int cur_var = targetPsi.VarList.size();

    vector<FuncType*> jastrow;
    DistanceTableData* d_table = NULL;
    cur = cur->xmlChildrenNode;
    int ng = 0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
      	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
      	string target_name((const char*)(xmlGetProp(cur,(const xmlChar *)"target")));
        map<string,ParticleSet*>::iterator pa_it(ptclPool.find(source_name));
        map<string,ParticleSet*>::iterator pb_it(ptclPool.find(target_name));

        if(pa_it == ptclPool.end()  || pb_it == ptclPool.end()) return false;

	ParticleSet* a = (*pa_it).second;
	ParticleSet* b = (*pb_it).second;
	d_table = DistanceTable::getTable(DistanceTable::add(*a,*b));
	ng = a->getSpeciesSet().getTotalNum();
	XMLReport("Number of sources " << ng)
	for(int i=0; i<ng; i++) jastrow.push_back(NULL);
      } else if(cname == corr_tag) {
	string speciesA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
	//string speciesB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
	int ia = d_table->origin().getSpeciesSet().findSpecies(speciesA);
	if(!(jastrow[ia])) {
	  jastrow[ia]= new FuncType;
	  jastrow[ia]->put(cur,targetPsi.VarList);
	  XMLReport("Added Jastrow Correlation between " << speciesA)
	}
      }
      cur = cur->next;
    } // while cur
    
    J1 = new JneType(targetPtcl, d_table);
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

    //set this jastrow function to be not optimizable
    if(targetPsi.VarList.size() == cur_var) {
      J1->setOptimizable(false);
    }
    targetPsi.add(J1);
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
    string jastfunction;
    xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
    if(ftype) jastfunction = (const char*) ftype;

    LOGMSG("")
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
    } else if(jasttype == "One-Body") {
      if(jastfunction == "nocusp") {
	LOGMSG("Jastrow Function: nucusp")
	OneBodyJastrow<NoCuspJastrow<ValueType> > *J1 = NULL;
	return createOneBody(cur,J1);
      } else if(jastfunction == "pade2") {
	LOGMSG("Jastrow Function: pade form (a*r+c*r^2)/(1+br).")
	OneBodyJastrow<PadeJastrow2<ValueType> > *J1 = NULL;
	return createOneBody(cur,J1);
      } else {
	LOGMSG("Jastrow Function: pade.")
	OneBodyJastrow<PadeJastrow<ValueType> > *J1 = NULL;
	return createOneBody(cur,J1);
      }
    } else if(jasttype == "Polarization") {
      PolarizedJastrow *jp=new PolarizedJastrow;
      if(jp) {
        //compare the size of VarList before/after reading in
        int cur_var = targetPsi.VarList.size();
        jp->put(cur,targetPsi.VarList);
        if(targetPsi.VarList.size() == cur_var) {
          jp->setOptimizable(false);
        }
        targetPsi.add(jp);
        return true;
      } else {
        return false;
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
