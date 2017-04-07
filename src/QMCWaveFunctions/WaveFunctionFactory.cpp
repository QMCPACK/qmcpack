//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file WaveFunctionFactory.cpp
 *@brief Definition of a WaveFunctionFactory
 */
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/IonOrbitalBuilder.h"

#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#else
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#if OHMMS_DIM==3 && QMC_BUILD_LEVEL>1 && !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
#endif

#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{

WaveFunctionFactory::WaveFunctionFactory(ParticleSet* qp, PtclPoolType& pset, Communicate* c)
  : MPIObjectBase(c)
  , targetPtcl(qp),ptclPool(pset),targetPsi(0), myNode(NULL)
{
  ClassName="WaveFunctionFactory";
  myName="psi0";
}

void WaveFunctionFactory::setPsi(TrialWaveFunction* psi)
{
  this->setName(psi->getName());
  targetPsi=psi;
}

bool WaveFunctionFactory::build(xmlNodePtr cur, bool buildtree)
{
  ReportEngine PRE(ClassName,"build");
  if(cur == NULL)
    return false;
  bool attach2Node=false;
  if(buildtree)
  {
    if(myNode == NULL)
    {
      myNode = xmlCopyNode(cur,1);
    }
    else
    {
      attach2Node=true;
    }
  }
  if(targetPsi==0) //allocate targetPsi and set the name
  {
    targetPsi  = new TrialWaveFunction(myComm);
    targetPsi->setName(myName);
    targetPsi->setMassTerm(*targetPtcl);
  }
  cur = cur->children;
  bool success=true;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname =="sposet_builder")
    {
      BasisSetFactory basisFactory(*targetPtcl,*targetPsi,ptclPool);
      basisFactory.build_sposet_collection(cur);
    }
    else if (cname == OrbitalBuilderBase::detset_tag)
    {
      success = addFermionTerm(cur);
      bool foundtwist(false);
      xmlNodePtr kcur = cur->children;
      while(kcur != NULL)
      {
        std::string kname((const char*)(kcur->name));
        if (kname=="h5tag")
        {
          std::string hdfName;
          OhmmsAttributeSet attribs;
          attribs.add (hdfName, "name");
          if (hdfName=="twistAngle")
          {
            std::vector<ParticleSet::RealType> tsts(3,0);
            putContent(tsts,kcur);
            targetPsi->setTwist(tsts);
            foundtwist=true;
          }
        }
        kcur=kcur->next;
      }
      if(!foundtwist)
      {
        //default twist is [0 0 0]
        std::vector<ParticleSet::RealType> tsts(3,0);
        targetPsi->setTwist(tsts);
      }
    }
    else if (cname ==  OrbitalBuilderBase::jastrow_tag)
    {
      OrbitalBuilderBase *jbuilder = new JastrowBuilder(*targetPtcl,*targetPsi,ptclPool);
      jbuilder->setReportLevel(ReportLevel);
      success = jbuilder->put(cur);
      addNode(jbuilder,cur);
    }
    else if (cname == OrbitalBuilderBase::ionorb_tag)
    {
      IonOrbitalBuilder *builder = new IonOrbitalBuilder
        (*targetPtcl, *targetPsi, ptclPool);
      success = builder->put(cur);
      addNode(builder,cur);
    }
    else if ((cname ==  "Molecular") || (cname =="molecular"))
    {
      APP_ABORT("  Removed Helium Molecular terms from qmcpack ");
    }
#if QMC_BUILD_LEVEL>2 && !defined(QMC_COMPLEX) && OHMMS_DIM==3
    else if(cname == "agp")
    {
      AGPDeterminantBuilder* agpbuilder = new AGPDeterminantBuilder(*targetPtcl,*targetPsi,ptclPool);
      success = agpbuilder->put(cur);
      addNode(agpbuilder,cur);
    }
#endif
    if(attach2Node)
      xmlAddChild(myNode,xmlCopyNode(cur,1));
    cur = cur->next;
  }
  //{
  //  ReportEngine PREA("TrialWaveFunction","print");
  //  targetPsi->VarList.print(app_log());
  //}
// synch all parameters. You don't want some being different if same name.
  opt_variables_type dummy;
  targetPsi->checkInVariables(dummy);
  dummy.resetIndex();
  targetPsi->checkOutVariables(dummy);
  targetPsi->resetParameters(dummy);
  return success;
}


bool WaveFunctionFactory::addFermionTerm(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"addFermionTerm");
  std::string orbtype("MolecularOrbital");
  std::string nuclei("i");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(orbtype,"type");
  oAttrib.add(nuclei,"source");
  oAttrib.put(cur);
  OrbitalBuilderBase* detbuilder=0;
  if(orbtype == "electron-gas")
  {
#if defined(QMC_COMPLEX)
    detbuilder = new ElectronGasComplexOrbitalBuilder(*targetPtcl,*targetPsi);
#else
    detbuilder = new ElectronGasOrbitalBuilder(*targetPtcl,*targetPsi);
#endif
  }
//#if OHMMS_DIM == 3 && QMC_BUILD_LEVEL>1
    else if(orbtype == "PWBasis" || orbtype == "PW" || orbtype == "pw")
    {
      detbuilder = new PWOrbitalBuilder(*targetPtcl,*targetPsi,ptclPool);
    }
//#endif /* QMC_BUILD_LEVEL>1 */
  else
    detbuilder = new SlaterDetBuilder(*targetPtcl,*targetPsi,ptclPool);
  detbuilder->setReportLevel(ReportLevel);
  detbuilder->put(cur);
  addNode(detbuilder,cur);
  return true;
}


bool WaveFunctionFactory::addNode(OrbitalBuilderBase* b, xmlNodePtr cur)
{
  psiBuilder.push_back(b);
  ///if(myNode != NULL) {
  ///  std::cout << ">>>> Adding " << (const char*)cur->name << std::endl;
  ///  xmlAddChild(myNode,xmlCopyNode(cur,1));
  ///}
  return true;
}

void WaveFunctionFactory::setCloneSize(int np)
{
  myClones.resize(np,0);
}

WaveFunctionFactory*
WaveFunctionFactory::clone(ParticleSet* qp, int ip, const std::string& aname)
{
  WaveFunctionFactory* aCopy= new WaveFunctionFactory(qp,ptclPool,myComm);
  //turn off the report for the clones
  aCopy->setReportLevel(0);
  aCopy->setName(aname);
  aCopy->build(myNode,false);
  myClones[ip]=aCopy;
  return aCopy;
}

WaveFunctionFactory::~WaveFunctionFactory()
{
  DEBUG_MEMORY("WaveFunctionFactory::~WaveFunctionFactory");
  delete_iter(psiBuilder.begin(),psiBuilder.end());
}

bool WaveFunctionFactory::put(xmlNodePtr cur)
{
  return build(cur,true);
}

void WaveFunctionFactory::reset() { }

//  bool WaveFunctionFactory::addJastrowTerm(xmlNodePtr cur) {
//    std::string jasttype("0");
//    std::string jastname("0");
//    std::string funcname("0");
//
//    OhmmsAttributeSet oAttrib;
//    oAttrib.add(jasttype,"type");
//    oAttrib.add(jastname,"name");
//    oAttrib.add(funcname,"function");
//    oAttrib.put(cur);
//
//    if(jasttype[0] == '0')
//    {
//      app_warning() << "  WaveFunctionFactory::addJastrowTerm missing type. Ignore " << jastname << std::endl;
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
//      app_log() << "\n  Using JastrowBasisBuilder for TESTING ONLY" << std::endl;
//      jbuilder=new JastrowBuilder(*targetPtcl,*targetPsi,ptclPool);
//    }
//    else if(jasttype == "Long-Range")
//    {
//      app_log() << "\n  Using JAAPBCBuilder for two-body jatrow TESTING ONLY" << std::endl;
//      jbuilder = new JAAPBCBuilder(*targetPtcl,*targetPsi);
//    }
//    else if(jasttype == "One-Body")
//    {
//      if(useSpline) {
//        app_log() << "\n  Using NJABBuilder for one-body jatrow with spline functions" << std::endl;
//        jbuilder = new NJABBuilder(*targetPtcl,*targetPsi,ptclPool);
//      } else {
//        app_log() << "\n  Using JABBuilder for one-body jatrow with analytic functions" << std::endl;
//        jbuilder = new JABBuilder(*targetPtcl,*targetPsi,ptclPool);
//      }
//    }
//#if !defined(QMC_COMPLEX)
//    else if(jasttype == "Three-Body-Geminal") {
//      app_log() << "\n  creating Three-Body-Germinal Jastrow function " << std::endl;
//      std::string source_name("i");
//      const xmlChar* iptr = xmlGetProp(cur, (const xmlChar *)"source");
//      if(iptr != NULL) source_name=(const char*)iptr;
//      PtclPoolType::iterator pit(ptclPool.find(source_name));
//      if(pit != ptclPool.end()) {
//        jbuilder = new ThreeBodyGeminalBuilder(*targetPtcl,*targetPsi,*((*pit).second));
//      }
//    } else if (jasttype == "Three-Body-Pade") {
//      app_log() << "\n  creating Three-Body-Pade Jastrow function " << std::endl;
//      std::string source_name("i");
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
//      app_warning() << "    " << jasttype << " is not valid." << std::endl;
//      return false;
//    }
//  }
}
