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
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/JAABuilder.h"
#include "QMCWaveFunctions/ModPadeJastrow.h"
#include "QMCWaveFunctions/TwoBodyJastrowOrbital.h"

namespace qmcplusplus {

   JAABuilder::JAABuilder(ParticleSet& p, TrialWaveFunction& psi):OrbitalBuilderBase(p,psi),
    IgnoreSpin(true){
    }
  /** Create a two-body Jatrow function with a template 
   *@param cur the current xmlNode 
   *@param dummy null pointer used to identify FN
   *
   *The template class JeeType is a functor which handles the
   *evaluation of the function value, gradients and laplacians using
   *distance tables. This is a specialized builder function for
   *spin-dependent Jastrow function,e.g., for electrons, two functions
   *are created for uu(dd) and ud(du).
   */
  template<class FN>
  bool JAABuilder::createJAA(xmlNodePtr cur, FN* dummy) {

    string corr_tag("correlation");
    int ng = targetPtcl.groups();
    //temporary storage for counting only
    map<string,FN*> j2Unique;
    vector<FN*> jastrow(ng*ng);
    for(int i=0; i<ng*ng; i++) jastrow[i]=0;

    int ia=0, ib=0, iab=0;
    int cur_var = targetPsi.VarList.size();
    xmlNodePtr gridPtr=NULL;
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == corr_tag) {
        string spA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
        string spB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
        const xmlChar* refptr=xmlGetProp(cur,(const xmlChar *)"ref");
        const xmlChar* idptr=xmlGetProp(cur,(const xmlChar *)"id");

        if(!IgnoreSpin) { //get the species
          ia = targetPtcl.getSpeciesSet().findSpecies(spA);
          ib = targetPtcl.getSpeciesSet().findSpecies(spB);
          iab = ia*ng+ib;
        }
        if(!(jastrow[iab])) {
          //create the new Jastrow function
          FN *j2= new FN(ia==ib);
          if(targetPtcl.Lattice.BoxBConds[0])
            j2->setDensity(targetPtcl.getTotalNum()/targetPtcl.Lattice.Volume);
          if(idptr == NULL) {
            ostringstream idassigned; idassigned << "j2"<<iab;
            j2Unique[idassigned.str()] = j2;
          } else {
            j2Unique[(const char*)idptr]=j2;
          }
          //initialize
          j2->put(cur,targetPsi.VarList);
          jastrow[iab]= j2;
          if(ia != ib) {//up-down-type pair, treat down-up the same way
            jastrow[ib*ng+ia] = j2;
          } else {
            for(int iaa=0; iaa<ng; iaa++) if(iaa != ia && jastrow[iaa*ng+iaa]==0) jastrow[iaa*ng+iaa] = j2;
          }
          XMLReport("Added Jastrow Correlation between "<<spA<<" and "<<spB)
        } else {
          ERRORMSG("Using an existing Jastrow Correlation "<<spA<<" and "<<spB)
        }
      }
      cur = cur->next;
    } // while cur

    DistanceTableData* d_table = DistanceTable::add(targetPtcl);
    /*
    if(IgnoreSpin) {
      TwoBodyJastrow<FN,true>* J2 = new TwoBodyJastrow<FN,true>(targetPtcl,d_table);
      J2->F=(*jastrow[0]);
      J2->setOptimizable(true);
      targetPsi.add(J2);
    } else {
      TwoBodyJastrow<FN,false>* J2 = new TwoBodyJastrow<FN,false>(targetPtcl,d_table);
      int ij = 0;    
      for(int i=0; i<ng; i++) {
        for(int j=0; j<ng; j++, ij++) {
          J2->F.push_back(jastrow[ij]);
        }
      }
      J2->setOptimizable(true);
      targetPsi.add(J2);
    }
    */
    typedef TwoBodyJastrowOrbital<FN> JeeType;
    JeeType *J2 = new JeeType(targetPtcl,d_table);
    J2->insert(j2Unique);
    if(IgnoreSpin) {
      LOGMSG("  Spin-indepdenent two-body jastrow ")
        for(int i=0; i<ng; i++) {
          for(int j=0; j<ng; j++) {
            J2->addFunc(jastrow[0]);
          }
        }
    } else {
      LOGMSG("  Spin-depdenent two-body jastrow. Unique functors = "<< j2Unique.size())
        int ij = 0;    
      for(int i=0; i<ng; i++) {
        for(int j=0; j<ng; j++, ij++) {
          J2->addFunc(jastrow[ij]);
        }
      }
    }

    //set this jastrow function to be not optimizable
    if(targetPsi.VarList.size() == cur_var) {
      J2->setOptimizable(false);
    }

    targetPsi.addOrbital(J2);
    XMLReport("Added a Two-Body Jastrow Function")
    return true;
  }

  bool JAABuilder::put(xmlNodePtr cur) {

    const xmlChar* spin=xmlGetProp(cur,(const xmlChar*)"spin");
    if(spin != NULL) {
      string a((const char*)spin);
      if(a == "yes") IgnoreSpin=false;
    }

    string jasttype((const char*)(xmlGetProp(cur, (const xmlChar *)"type")));
    if(jasttype == "Two-Body-Spin") {
      //only for the backward compability
      IgnoreSpin=false;
    }

    string jastfunction("pade");
    const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
    if(ftype) jastfunction = (const char*) ftype;

    bool success=false;
    //if(jastfunction == "pade") {
    //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
    //  PadeJastrow<RealType> *dummy = 0;
    //  success = createJAA(cur,dummy);
    //} else 
    if(jastfunction == "short") {
      app_log() << "  Modified Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
      IgnoreSpin=true;
      ModPadeJastrow<RealType> *dummy = 0;
      success = createJAA(cur,dummy);
    }
    //} else if(jastfunction == "rpa") {
    //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
    //  RPAJastrow<RealType> *dummy = 0;
    //  success = createJAA(cur,dummy);
    //}
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
