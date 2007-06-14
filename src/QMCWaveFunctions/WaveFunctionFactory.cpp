/////////////////////////////////////////////////////////////////
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
/**@file WaveFunctionFactory.cpp
 *@brief Definition of a WaveFunctionFactory 
 */
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#else
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
#endif
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus {
  WaveFunctionFactory::WaveFunctionFactory(ParticleSet* qp, PtclPoolType& pset): 
    targetPtcl(qp),ptclPool(pset),targetPsi(0), myNode(NULL) {
      myName="psi0";
    }

  bool WaveFunctionFactory::build(xmlNodePtr cur, bool buildtree) {

    if(cur == NULL) return false;

    bool attach2Node=false;
    if(buildtree) {
      if(myNode == NULL) {
        myNode = xmlCopyNode(cur,1);
      } else {
        attach2Node=true;
      }
    }


    if(targetPsi==0) {//allocate targetPsi and set the name
      targetPsi  = new TrialWaveFunction;
      targetPsi->setName(myName);
      app_log() << "  Creating a trial wavefunction " << myName << endl;
    }

    cur = cur->children;
    bool success=true;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if (cname == OrbitalBuilderBase::detset_tag) 
      {
        success = addFermionTerm(cur);
      } 
      else if (cname ==  OrbitalBuilderBase::jastrow_tag) 
      {
        OrbitalBuilderBase *jbuilder = new JastrowBuilder(*targetPtcl,*targetPsi,ptclPool);
        success = jbuilder->put(cur);
        addNode(jbuilder,cur);
      }
      else if(cname == "agp") 
      {
#if defined(QMC_COMPLEX)
        app_error() << "  AGPDeterminant cannot be used with QMC_COMPLEX=1" << endl;
        return false;
#else
        AGPDeterminantBuilder* agpbuilder = new AGPDeterminantBuilder(*targetPtcl,*targetPsi,ptclPool);
        success = agpbuilder->put(cur);
        addNode(agpbuilder,cur);
#endif
      } 
      if(attach2Node) xmlAddChild(myNode,xmlCopyNode(cur,1));
      cur = cur->next;
    }

    if(OrbitalBuilderBase::print_level>0)
    {
      app_log() << "  List of optimizable variables " << endl;
      targetPsi->VarList.print(app_log());

      //set to zero so that nothing is written again
      OrbitalBuilderBase::print_level=0;
    }

    return success;
  }

  bool WaveFunctionFactory::addFermionTerm(xmlNodePtr cur) {

    OrbitalBuilderBase* detbuilder=0;

    string orbtype("MolecularOrbital");
    string nuclei("i");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(orbtype,"type");
    oAttrib.add(nuclei,"source");
    oAttrib.put(cur);

    app_log() << "\n  Slater determinant terms using " << orbtype << endl;

#if defined(QMC_COMPLEX)
    if(orbtype == "electron-gas") 
    {
      detbuilder = new ElectronGasComplexOrbitalBuilder(*targetPtcl,*targetPsi);
    } 
#else
    if(orbtype == "electron-gas") 
    {
      detbuilder = new ElectronGasOrbitalBuilder(*targetPtcl,*targetPsi);
    } 
#endif
    else if(orbtype == "PWBasis" || orbtype == "PW" || orbtype == "pw") 
    {
      detbuilder = new PWOrbitalBuilder(*targetPtcl,*targetPsi);
    } 
    //else if(orbtype == "MolecularOrbital") 
    //{
    //  detbuilder = new MolecularOrbitalBuilder(*targetPtcl,*targetPsi,ptclPool);
    //} 
    else 
    {
      app_log() << "  Creating concrete SlaterDeterminant class with SlaterDetBuilder." << endl;
      detbuilder = new SlaterDetBuilder(*targetPtcl,*targetPsi,ptclPool);
    }

    if(detbuilder) {//valid determinant set
      detbuilder->put(cur);
      addNode(detbuilder,cur);
      return true;
    } else {
      return false;
    }
  }


  bool WaveFunctionFactory::addNode(OrbitalBuilderBase* b, xmlNodePtr cur) {
    psiBuilder.push_back(b);
    ///if(myNode != NULL) {
    ///  cout << ">>>> Adding " << (const char*)cur->name << endl;
    ///  xmlAddChild(myNode,xmlCopyNode(cur,1));
    ///}
    return true;
  }

  void WaveFunctionFactory::setCloneSize(int np) {
    myClones.resize(np,0);
  }

  WaveFunctionFactory*
  WaveFunctionFactory::clone(ParticleSet* qp, int ip, const string& aname) {
    WaveFunctionFactory* aCopy= new WaveFunctionFactory(qp,ptclPool);
    aCopy->setName(aname);
    aCopy->build(myNode,false);
    myClones[ip]=aCopy;
    return aCopy;
  }

  WaveFunctionFactory::~WaveFunctionFactory() {
    //clean up everything
  }

  bool WaveFunctionFactory::get(std::ostream& ) const {
    return true;
  }

  bool WaveFunctionFactory::put(std::istream& ) {
    return true;
  }

  bool WaveFunctionFactory::put(xmlNodePtr cur) {
    return build(cur,true);
  }

  void WaveFunctionFactory::reset() { }

//  bool WaveFunctionFactory::addJastrowTerm(xmlNodePtr cur) {
//    string jasttype("0");
//    string jastname("0");
//    string funcname("0");
//
//    OhmmsAttributeSet oAttrib;
//    oAttrib.add(jasttype,"type");
//    oAttrib.add(jastname,"name");
//    oAttrib.add(funcname,"function");
//    oAttrib.put(cur);
//
//    if(jasttype[0] == '0')
//    {
//      app_warning() << "  WaveFunctionFactory::addJastrowTerm missing type. Ignore " << jastname << endl;
//      return false;
//    }
//
//    //string jasttype((const char*)(xmlGetProp(cur, (const xmlChar *)"type")));
//    //string jastname((const char*)(xmlGetProp(cur, (const xmlChar *)"name")));
//    //string funcname((const char*)(xmlGetProp(cur, (const xmlChar *)"function")));
//    bool useSpline=false;
//    const xmlChar* gptr=xmlGetProp(cur,(const xmlChar*)"transform");
//    if(gptr != NULL) {
//      if(xmlStrEqual(gptr,(const xmlChar*)"yes")) {
//        useSpline=true;
//      } 
//    }
//
//    OrbitalBuilderBase* jbuilder=0;
//    if(jasttype.find("Two") < jasttype.size()) 
//    {
//      jbuilder=new TwoBodyJastrowBuilder(*targetPtcl,*targetPsi,ptclPool);
//    } 
//    else if(jasttype == "TEST")
//    {
//      app_log() << "\n  Using JastrowBasisBuilder for TESTING ONLY" << endl;
//      jbuilder=new JastrowBuilder(*targetPtcl,*targetPsi,ptclPool);
//    }
//    else if(jasttype == "Long-Range") 
//    {
//      app_log() << "\n  Using JAAPBCBuilder for two-body jatrow TESTING ONLY" << endl;
//      jbuilder = new JAAPBCBuilder(*targetPtcl,*targetPsi);
//    } 
//    else if(jasttype == "One-Body") 
//    {
//      if(useSpline) {
//        app_log() << "\n  Using NJABBuilder for one-body jatrow with spline functions" << endl;
//        jbuilder = new NJABBuilder(*targetPtcl,*targetPsi,ptclPool);
//      } else {
//        app_log() << "\n  Using JABBuilder for one-body jatrow with analytic functions" << endl;
//        jbuilder = new JABBuilder(*targetPtcl,*targetPsi,ptclPool);
//      }
//    } 
//#if !defined(QMC_COMPLEX)
//    else if(jasttype == "Three-Body-Geminal") {
//      app_log() << "\n  creating Three-Body-Germinal Jastrow function " << endl;
//      string source_name("i");
//      const xmlChar* iptr = xmlGetProp(cur, (const xmlChar *)"source");
//      if(iptr != NULL) source_name=(const char*)iptr;
//      PtclPoolType::iterator pit(ptclPool.find(source_name));
//      if(pit != ptclPool.end()) {
//        jbuilder = new ThreeBodyGeminalBuilder(*targetPtcl,*targetPsi,*((*pit).second));
//      }
//    } else if (jasttype == "Three-Body-Pade") {
//      app_log() << "\n  creating Three-Body-Pade Jastrow function " << endl;
//      string source_name("i");
//      const xmlChar* iptr = xmlGetProp(cur, (const xmlChar *)"source");
//      //if(iptr != NULL) source_name=(const char*)iptr;
//      PtclPoolType::iterator pit(ptclPool.find(source_name));
//      if(pit != ptclPool.end()) {
//        jbuilder = new ThreeBodyPadeBuilder(*targetPtcl,*targetPsi,*((*pit).second));
//      }
//    }
//#endif
//
//    if(jbuilder) {
//      jbuilder->put(cur);
//      addNode(jbuilder,cur);
//      return true;
//    } else {
//      app_warning() << "    " << jasttype << " is not valid." << endl;
//      return false;
//    }
//  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
