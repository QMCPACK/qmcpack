//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/PadeJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"

namespace qmcplusplus
{

PadeJastrowBuilder::PadeJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                                       PtclPoolType& psets):
  OrbitalBuilderBase(target,psi),ptclPool(psets)
{
  ClassName="PadeJastrowBuilder";
}

bool PadeJastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put()");
  std::string sourceOpt=targetPtcl.getName();
  std::string jname="PadeJastrow";
  std::string spin="yes";
  std::string id_b="jee_b";
  RealType pade_b=1.0;
  OhmmsAttributeSet pattrib;
  pattrib.add(jname,"name");
  pattrib.add(spin,"spin");
  pattrib.add(sourceOpt,"source");
  pattrib.put(cur);
  //bool spindep=(spin=="yes");
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd=species.addAttribute("charge");
  typedef PadeFunctor<RealType> RadFuncType;
  if(sourceOpt == targetPtcl.getName())
  {
    //two-body
#if defined(ENABLE_SOA)
    typedef J2OrbitalSoA<RadFuncType> J2Type;
#else
    typedef TwoBodyJastrowOrbital<RadFuncType> J2Type;
#endif
    typedef DiffTwoBodyJastrowOrbital<RadFuncType> dJ2Type;
    int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
    J2Type *J2 = new J2Type(targetPtcl,taskid);
    dJ2Type *dJ2 = new dJ2Type(targetPtcl);
    cur= cur->xmlChildrenNode;
    while(cur!=NULL)
    {
      std::string cname((const char*)cur->name);
      if (cname == "correlation")
      {
        OhmmsAttributeSet rAttrib;
        RealType cusp=-1e10;
        std::string spA(species.speciesName[0]);
        std::string spB(species.speciesName[0]);
        rAttrib.add(spA,"speciesA");
        rAttrib.add(spB,"speciesB");
        rAttrib.add(cusp,"cusp");
        rAttrib.put(cur);
        int ia = species.findSpecies(spA);
        int ib = species.findSpecies(spB);
        if(ia==species.size() || ib == species.size())
        {
          PRE.error("Failed. Species are incorrect.",true);
        }
        if(cusp<-1e6)
        {
          RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
          cusp = (ia==ib)? -0.25*qq:-0.5*qq;
        }
        std::ostringstream o;
        o<<"j2"<<ia<<ib;
        RadFuncType *functor = new RadFuncType(cusp,o.str());
        functor->put(cur);
        J2->addFunc(ia,ib,functor);
        dJ2->addFunc(ia,ib,functor);
        //{
        //  std::ostringstream o;
        //  o<< "pade"<<ia<<"-"<<ib;
        //  std::ofstream fstream(o.str().c_str());
        //  int n=100;
        //  RealType d=10/100.,r=0.001;
        //  RealType u,du,d2du;
        //  for (int i=0; i<n; ++i)
        //  {
        //    u=functor->evaluate(r,du,d2du);
        //    fstream << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du
        //      << std::setw(22) << d2du << std::endl;
        //    r+=d;
        //  }
        //}
      }
      cur=cur->next;
    }
    J2->dPsi=dJ2;
    targetPsi.addOrbital(J2,"J2_pade");
    J2->setOptimizable(true);
  }
  return true;
}
//  bool PadeJastrowBuilder::put(xmlNodePtr cur)
//  {
//
//    ReportEngine PRE(ClassName,"put()");
//
//    std::string sourceOpt=targetPtcl.getName();
//    std::string jname="PadeJastrow";
//    std::string spin="no";
//    std::string id_b="jee_b";
//    RealType pade_b=1.0;
//    OhmmsAttributeSet pattrib;
//    pattrib.add(jname,"name");
//    pattrib.add(spin,"spin");
//    pattrib.add(sourceOpt,"source");
//    pattrib.put(cur);
//
//    cur=cur->children;
//    while(cur != NULL)
//    {
//      {//just to hide this
//        std::string pname="0";
//        OhmmsAttributeSet aa;
//        aa.add(pname,"name");
//        aa.add(id_b,"id");
//        aa.put(cur);
//        if(pname[0]=='B') putContent(pade_b,cur);
//      }
//
//      xmlNodePtr cur1=cur->children;
//      while(cur1!= NULL)
//      {
//        std::string pname="0";
//        OhmmsAttributeSet aa;
//        aa.add(pname,"name");
//        aa.add(id_b,"id");
//        aa.put(cur1);
//        if(pname[0]=='B') putContent(pade_b,cur1);
//        cur1=cur1->next;
//        std::cout << "TESTING " << pade_b << std::endl;
//      }
//      cur=cur->next;
//    }
//
//    app_log() << "PadeJastrowBuilder " << id_b << " = " << pade_b << std::endl;
//
//    typedef PadeFunctor<RealType> FuncType;
//
//    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
//    JeeType *J2 = new JeeType(targetPtcl,targetPsi.is_manager());
//
//    SpeciesSet& species(targetPtcl.getSpeciesSet());
//    RealType q=species(0,species.addAttribute("charge"));
//
//    if(spin == "no")
//    {//
//      RealType cusp=-0.5*q*q;
//      FuncType *func=new FuncType(cusp,pade_b);
//      func->setIDs("jee_cusp",id_b,false,true);//set the ID's, fixed cusp
//      J2->addFunc(0,0,func);
//      //DerivFuncType *dfunc=new DerivFuncType(cusp,B);
//      //dJ2->addFunc("pade_uu",0,0,dfunc);
//      //dFuncList.push_back(dfunc);
//      app_log() << "    Adding Spin-independent Pade Two-Body Jastrow Cusp " << cusp<< "\n";
//    }
//    else
//    {
//      //build uu functor
//      RealType cusp_uu=-0.25*q*q;
//      FuncType *funcUU=new FuncType(cusp_uu,pade_b);
//      funcUU->setIDs("pade_uu",id_b);//set the ID's
//
//      //build ud functor
//      RealType cusp_ud=-0.5*q*q;
//      FuncType *funcUD=new FuncType(cusp_ud,pade_b);
//      funcUD->setIDs("pade_ud",id_b);//set the ID's
//
//      J2->addFunc(0,0,funcUU);
//      J2->addFunc(0,1,funcUD);
//
//      //DerivFuncType *dfuncUU=new DerivFuncType(cusp_uu,B);
//      //DerivFuncType *dfuncUD=new DerivFuncType(cusp_ud,B);
//      //dJ2->addFunc("pade_uu",0,0,dfuncUU);
//      //dJ2->addFunc("pade_ud",0,1,dfuncUD);
//      app_log() << "    Adding Spin-dependent Pade Two-Body Jastrow " << "\n";
//      app_log() << "      parallel spin     " << cusp_uu << "\n";
//      app_log() << "      antiparallel spin " << cusp_ud << "\n";
//    }
//
//    targetPsi.addOrbital(J2,"J2_pade");
//
//    if(sourceOpt != targetPtcl.getName())
//    {
//      std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
//      if(pa_it == ptclPool.end())
//      {
//        PRE.warning("PadeJastrowBuilder::put failed. "+sourceOpt+" does not exist.");
//        return true;
//      }
//      ParticleSet& sourcePtcl= (*(*pa_it).second);
//
//      app_log() << "  PadeBuilder::Adding Pade One-Body Jastrow with effective ionic charges." << std::endl;
//      typedef OneBodyJastrowOrbital<FuncType> JneType;
//      JneType* J1 = new JneType(sourcePtcl,targetPtcl);
//
//      //typedef OneBodyJastrowOrbital<DerivFuncType> DerivJneType;
//      //DerivJneType* dJ1=new DerivJneType(sourcePtcl,targetPtcl);
//
//      SpeciesSet& Species(sourcePtcl.getSpeciesSet());
//      int ng=Species.getTotalNum();
//      int icharge = Species.addAttribute("charge");
//      for(int ig=0; ig<ng; ++ig)
//      {
//        RealType zeff=Species(icharge,ig);
//        std::ostringstream j1id;
//        j1id<<"pade_"<<Species.speciesName[ig];
//
//        RealType sc=std::pow(2*zeff,0.25);
//        FuncType *func=new FuncType(-zeff,pade_b,sc);
//        func->setIDs(j1id.str(),id_b);
//
//        J1->addFunc(ig,func);
//
//        //DerivFuncType *dfunc=new DerivFuncType(-zeff,B,sc);
//        //dJ1->addFunc(ig,dfunc);
//        //dFuncList.push_back(dfunc);
//
//        app_log() << "    " << Species.speciesName[ig] <<  " Zeff = " << zeff << " B= " << pade_b*sc << std::endl;
//      }
//      targetPsi.addOrbital(J1,"J1_pade");
//    }
//    return true;
//  }

}
