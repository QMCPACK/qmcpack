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
#include "QMCWaveFunctions/Jastrow/WMConstraints.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  template<class T> bool WMFunctor<T>::put(xmlNodePtr cur)
  {
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID,"id"); //rAttrib.add(a->B0,"b");
    rAttrib.put(cur);
    return putContent(B0,cur);
  }

  WMConstraints::~WMConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool WMConstraints::put(xmlNodePtr cur) {
    //get generic parameters and grid information
    bool success=getVariables(cur);
    return success;
  }

  void WMConstraints::apply() {
    for(int i=0; i<FuncList.size(); i++) {
      FuncList[i]->reset();
    }
  }

  void WMConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
    map<string,BasisSetType*>::iterator it(myBasisSet.begin()),it_end(myBasisSet.end());
    while(it != it_end) {
      BasisSetType& basis(*((*it).second));
      for(int ib=0; ib<basis.size(); ib++) {
        basis[ib]->addOptimizables(outVars);
      }
      ++it;
    }
    for(int i=0; i<InFuncList.size(); i++) {
      InFuncList[i]->addOptimizables(outVars);
    }
  }

  void WMConstraints::addBasisGroup(xmlNodePtr cur) {
    string sourceOpt("e");
    string elementType("e");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(sourceOpt,"source");
    aAttrib.add(elementType,"elementType");
    aAttrib.put(cur);

    RealType rcut(myGrid->rmax());
    map<string,BasisSetType*>::iterator it(myBasisSet.find(elementType));
    if(it == myBasisSet.end()) {
       BasisSetType* newBasis=new BasisSetType;
       cur=cur->children;
       while(cur != NULL) {
         string cname((const char*)(cur->name));
         if(cname == "parameter") {
           //BasisType* a=new BasisType(1.0,rcut);
           WMFunctor<RealType>* a = new WMFunctor<RealType>(1.0,rcut);
           a->put(cur);
           newBasis->push_back(a);
         }
         cur=cur->next;
       }
       //add a new BasisSet
       myBasisSet[elementType]=newBasis;
    }
  }


  OrbitalBase* WMConstraints::createTwoBody(ParticleSet& target) {

    InFuncType* infunc=0;
    setRadialGrid(target);
    xmlNodePtr cur=myNode->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "basisGroup") {
        addBasisGroup(cur);
      } else if(cname =="correlation") {
        string speciesA(target.getName());
        string speciesB(target.getName());
        OhmmsAttributeSet aAttrib;
        aAttrib.add(speciesA,"speciesA");
        aAttrib.add(speciesB,"speciesB");
        aAttrib.put(cur);
        if(speciesA==speciesB) {
          map<string,BasisSetType*>::iterator it(myBasisSet.find(speciesA));
          if(it == myBasisSet.end()) {
            app_error() <<  "  WMBasisSet for " << speciesA << " does not exist." << endl;
            continue;
          } 
          app_log() << "    Creating a correlation function = " << speciesA << "-" << speciesB << endl;
          infunc = createCorrelation(cur,(*it).second);
        }
      }
      cur=cur->next;
    }

    if(infunc==0) return 0;

    InFuncList.push_back(infunc);
    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(target);
    FuncType* nfunc= new FuncType(infunc,myGrid);
    for(int i=0; i<4; i++) J2->addFunc(nfunc);
    FuncList.push_back(nfunc);
    return J2;
  }

  OrbitalBase* WMConstraints::createOneBody(ParticleSet& target, ParticleSet& source) {
    vector<InFuncType*> jnSet;
    jnSet.resize(source.getSpeciesSet().getTotalNum(),0);
    xmlNodePtr cur=myNode->children;
    bool noOneBody=true;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "basisGroup") {
        addBasisGroup(cur);
      } else if(cname =="correlation") {
        string speciesA("e");
        string speciesB("e");
        OhmmsAttributeSet aAttrib;
        aAttrib.add(speciesA,"speciesA");
        aAttrib.add(speciesB,"speciesB");
        aAttrib.put(cur);
        if(speciesA != speciesB) {
          map<string,BasisSetType*>::iterator it(myBasisSet.find(speciesA));
          if(it == myBasisSet.end()) {
            app_error() <<  "  WMBasisSet for " << speciesA << " does not exist." << endl;
            continue;
          }
          app_log() << "    Creating a correlation function = " << speciesA << "-" << speciesB << endl;
          int gid=source.getSpeciesSet().addSpecies(speciesA);
          jnSet[gid] = createCorrelation(cur,(*it).second);
          noOneBody=false;
        }
      }
      cur=cur->next;
    }

    if(noOneBody) return 0;

    typedef OneBodyJastrow<FuncType> JneType;
    JneType* jne=new JneType(source,target);
    for(int ig=0; ig<jnSet.size(); ig++) {
      if(jnSet[ig]) {
        FuncType* nfunc= new FuncType(jnSet[ig],myGrid);
        jne->addFunc(ig,nfunc);
        FuncList.push_back(nfunc);
        InFuncList.push_back(jnSet[ig]);
      }
    }
    return jne;
  }

  WMConstraints::InFuncType*
    WMConstraints::createCorrelation(xmlNodePtr cur,BasisSetType* basis) {
      int nc=0;
      InFuncType* acombo=new InFuncType;
      cur=cur->children;
      while(cur != NULL) {
        string cname((const char*)(cur->name));
        if(cname == "parameter") {
          string id("0");
          string ref("0");
          RealType c=1.0;
          OhmmsAttributeSet aAttrib;
          aAttrib.add(id,"id");
          aAttrib.add(ref,"ref");
          aAttrib.put(cur);
          putContent(c,cur);
          if(nc <basis->size()) acombo->add((*basis)[nc++],c,id);
        }
        cur=cur->next;
      }

      if(nc)
        return acombo;
      else {
        delete acombo; 
        return 0;
      }
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
