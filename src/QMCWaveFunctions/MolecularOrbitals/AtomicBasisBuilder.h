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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ATOMICORBITALBUILDER_H
#define QMCPLUSPLUS_ATOMICORBITALBUILDER_H

#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  template<class RFB>
  struct AtomicBasisBuilder: public BasisSetBuilder {

    typedef typename RFB::CenteredOrbitalType COT;

    enum {DONOT_EXPAND=0, GAUSSIAN_EXPAND=1, NATURAL_EXPAND};

    RFB radFuncBuilder;

    bool addsignforM;
    int expandlm;
    string Morder;
    string sph;
    string basisType;
    string elementType;

    ///map for the radial orbitals
    map<string,int>  RnlID;

    ///map for (n,l,m,s) to its quantum number index
    map<string,int> nlms_id;

    AtomicBasisBuilder(const string& eName);

    bool put(xmlNodePtr cur);

    SPOSetBase* createSPOSet(xmlNodePtr cur) { return 0;}

    COT* createAOSet(xmlNodePtr cur);

    int expandYlm(const string& rnl, const QuantumNumberType& nlms, 
        int num, COT* aos, xmlNodePtr cur1, int expandlm=DONOT_EXPAND);

  };

  template<class RFB>
  AtomicBasisBuilder<RFB>::AtomicBasisBuilder(const string& eName):
    addsignforM(false), expandlm(GAUSSIAN_EXPAND), Morder("gaussian"),
    sph("default"), basisType("Numerical"), elementType(eName)
  {
    nlms_id["n"] = q_n;
    nlms_id["l"] = q_l;
    nlms_id["m"] = q_m;
    nlms_id["s"] = q_s;
  }

  template<class RFB>
  bool AtomicBasisBuilder<RFB>::put(xmlNodePtr cur) {
    //Register valid attributes attributes
    OhmmsAttributeSet aAttrib;
    //aAttrib.add(elementType,"elementType"); aAttrib.add(elementType,"species");
    aAttrib.add(basisType,"type");
    aAttrib.add(sph,"angular"); 
    aAttrib.add(addsignforM,"expM"); 
    aAttrib.add(Morder,"expandYlm"); 
    aAttrib.put(cur);

    if(sph == "spherical") addsignforM=1; //include (-1)^m
    if(Morder == "gaussian") {
      expandlm = GAUSSIAN_EXPAND;
    } else if(Morder == "natural"){
      expandlm = NATURAL_EXPAND;
    } else if(Morder == "no") {
      expandlm = DONOT_EXPAND;
    }

    return radFuncBuilder.putCommon(cur);
  }

  template<class RFB>
  typename AtomicBasisBuilder<RFB>::COT*
  AtomicBasisBuilder<RFB>::createAOSet(xmlNodePtr cur) {

    app_log() << "  AO BasisSet for " << elementType << "\n";
    if(addsignforM) 
      app_log() << "   Spherical Harmonics contain (-1)^m factor" << endl;
    else
      app_log() << "   Spherical Harmonics  DO NOT contain (-1)^m factor" << endl;

    switch(expandlm) {
      case(GAUSSIAN_EXPAND):
      app_log() << "   Angular momentum m expanded addocrding to Gaussian" << endl;
      break;
      case(NATURAL_EXPAND):
      app_log() << "   Angular momentum m expanded as -l, ... ,l" << endl;
      break;
      default:
      app_log() << "   Angular momentum m is explicitly given." << endl;
    }

    QuantumNumberType nlms;
    string rnl;

    int Lmax(0); //maxmimum angular momentum of this center
    int num(0);//the number of localized basis functions of this center

    //process the basic property: maximun angular momentum, the number of basis functions to be added
    vector<xmlNodePtr> radGroup;
    xmlNodePtr cur1 = cur->xmlChildrenNode;
    xmlNodePtr gptr=0;
    while(cur1 != NULL) {
      string cname1((const char*)(cur1->name));
      if(cname1 == "basisGroup") {
        radGroup.push_back(cur1);
        int l=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"l")));
        Lmax = max(Lmax,l);
        //expect that only Rnl is given
        if(expandlm) 
          num += 2*l+1;
        else
          num++;
      } else if(cname1 == "grid") {
        gptr = cur1;
      }
      cur1 = cur1->next;
    }
    //create a new set of atomic orbitals sharing a center with (Lmax, num)
    //if(addsignforM) the basis function has (-1)^m sqrt(2)Re(Ylm)
    COT* aos = new COT(Lmax,addsignforM);
    aos->LM.resize(num);
    aos->NL.resize(num);

    //Now, add distinct Radial Orbitals and (l,m) channels
    num=0;
    radFuncBuilder.setOrbitalSet(aos,elementType); //assign radial orbitals for the new center
    radFuncBuilder.addGrid(gptr); //assign a radial grid for the new center

    vector<xmlNodePtr>::iterator it(radGroup.begin());
    vector<xmlNodePtr>::iterator it_end(radGroup.end());
    while(it != it_end) {
      cur1 = (*it);
      xmlAttrPtr att = cur1->properties;
      while(att != NULL) {
        string aname((const char*)(att->name));
        if(aname == "rid" || aname == "id") { //accept id/rid
          rnl = (const char*)(att->children->content);
        } else { 
          map<string,int>::iterator iit = nlms_id.find(aname);
          if(iit != nlms_id.end()) { //valid for n,l,m,s
            nlms[(*iit).second] = atoi((const char*)(att->children->content));
          } 
        }
        att = att->next;
      }
      //add Ylm channels
      app_log() << "   R(n,l,m,s) " << nlms[0] << " " << nlms[1] << " " << nlms[2] << " " << nlms[3] << endl;
      num = expandYlm(rnl,nlms,num,aos,cur1,expandlm);
      ++it;
    }

    aos->setBasisSetSize(-1);
    app_log() << "   Maximu Angular Momentum   = " << aos->Ylm.Lmax << endl
              << "   Number of Radial functors = " << aos->Rnl.size() << endl
              << "   Basis size                = " << aos->getBasisSetSize() << "\n\n";
    return aos;
  }


  template<class RFB>
  int AtomicBasisBuilder<RFB>::expandYlm(const string& rnl, const QuantumNumberType& nlms, int num, 
                  COT* aos, xmlNodePtr cur1, int expandlm) {
      if(expandlm == GAUSSIAN_EXPAND) {
        DEBUGMSG("Expanding Ylm according to Gaussian98")
        map<string,int>::iterator rnl_it = RnlID.find(rnl);
        if(rnl_it == RnlID.end()) {
          int nl = aos->Rnl.size();
          if(radFuncBuilder.addRadialOrbital(cur1,nlms)) {
            RnlID[rnl] = nl;
            int l = nlms[q_l];
            DEBUGMSG("Adding " << 2*l+1 << " spherical orbitals for l= " << l)
            switch (l) 
            {
              case(0):
                aos->LM[num] = aos->Ylm.index(0,0);  aos->NL[num] = nl; num++;
                break;
              case(1)://px(1),py(-1),pz(0)            
                aos->LM[num] = aos->Ylm.index(1,1);  aos->NL[num] = nl; num++;
                aos->LM[num] = aos->Ylm.index(1,-1); aos->NL[num] = nl; num++;
                aos->LM[num] = aos->Ylm.index(1,0);  aos->NL[num] = nl; num++;
                break; 
              default://0,1,-1,2,-2,...,l,-l
                aos->LM[num] = aos->Ylm.index(l,0);  aos->NL[num] = nl; num++;
                for(int tm=1; tm<=l; tm++) {
                  aos->LM[num] = aos->Ylm.index(l,tm);  aos->NL[num] = nl; num++;
                  aos->LM[num] = aos->Ylm.index(l,-tm); aos->NL[num] = nl; num++;
                }
                break;
            }
          }
        }
      } else if(expandlm == NATURAL_EXPAND) {
        DEBUGMSG("Expanding Ylm as -l,-l+1,...,l-1,l")
        map<string,int>::iterator rnl_it = RnlID.find(rnl);
        if(rnl_it == RnlID.end()) {
          int nl = aos->Rnl.size();
          if(radFuncBuilder.addRadialOrbital(cur1,nlms)) {
            RnlID[rnl] = nl;
            int l = nlms[q_l];
            LOGMSG("Adding " << 2*l+1 << " spherical orbitals")
            for(int tm=-l; tm<=l; tm++,num++) {
              aos->LM[num] = aos->Ylm.index(l,tm);  aos->NL[num] = nl;
            }
          }
        }
      } else {
        DEBUGMSG("Ylm is NOT expanded.")
        //assign the index for real Spherical Harmonic with (l,m)
        aos->LM[num] = aos->Ylm.index(nlms[q_l],nlms[q_m]);
        //radial orbitals: add only distinct orbitals
        map<string,int>::iterator rnl_it = RnlID.find(rnl);
        if(rnl_it == RnlID.end()) {
          int nl = aos->Rnl.size();
          if(radFuncBuilder.addRadialOrbital(cur1,nlms)) { //assign the index for radial orbital with (n,l)
            aos->NL[num] = nl;
            RnlID[rnl] = nl;
          }
        } else {
          //assign the index for radial orbital with (n,l) if repeated
          DEBUGMSG("Already added radial function for id: " << rnl)
          aos->NL[num] = (*rnl_it).second;
        }
        //increment number of basis functions
        num++;
      }
      return num;

    }

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
