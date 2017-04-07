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
    
    
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
//#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
//#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"

namespace qmcplusplus
{

////////////////////////////////////////
//PadeConstraints definitions
////////////////////////////////////////
PadeConstraints::PadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
  OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
{
  ClassName="PadeConstraints";
  JComponent.set(MULTIPLE);
  JComponent.set(ONEBODY);
  JComponent.set(TWOBODY);
  //dPsi=new AnalyticDiffOrbital(p.getTotalNum(),0);
}

PadeConstraints::~PadeConstraints()
{
  delete_iter(FuncList.begin(), FuncList.end());
  //delete_iter(dFuncList.begin(), dFuncList.end());
}

bool PadeConstraints::put(xmlNodePtr cur)
{
  //always copy the node
  bool success=getVariables(cur);
  std::map<std::string,std::pair<std::string,RealType> >::iterator vit(inVars.find("B"));
  if(vit == inVars.end())
    return false; //disaster, need to abort
  ID_B=(*vit).second.first;
  B=(*vit).second.second;
  node_B=cur; //save the last B
  return true;
}

void PadeConstraints::resetParameters(const opt_variables_type& active)
{
  //nothing to do
}

OrbitalBase* PadeConstraints::createTwoBody()
{
  ReportEngine PRE(ClassName,"createTwoBody()");
  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
  JeeType *J2 = new JeeType(targetPtcl);
  //typedef DPadeDBFunctor<RealType> DerivFuncType;
  //typedef TwoBodyJastrowOrbital<DerivFuncType> DerivJeeType;
  //DerivJeeType *dJ2= new DerivJeeType(targetPtcl);
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  RealType q=species(0,species.addAttribute("charge"));
  if(IgnoreSpin)
  {
    RealType cusp=-0.5*q*q;
    FuncType *func=new FuncType(cusp,B);
    func->setIDs("jee_cusp",ID_B);//set the ID's
    J2->addFunc("pade_uu",0,0,func);
    FuncList.push_back(func);
    //DerivFuncType *dfunc=new DerivFuncType(cusp,B);
    //dJ2->addFunc("pade_uu",0,0,dfunc);
    //dFuncList.push_back(dfunc);
    app_log() << "    Adding Spin-independent Pade Two-Body Jastrow Cusp " << cusp<< "\n";
  }
  else
  {
    //build uu functor
    RealType cusp_uu=-0.25*q*q;
    FuncType *funcUU=new FuncType(cusp_uu,B);
    funcUU->setIDs("pade_uu",ID_B);//set the ID's
    //build ud functor
    RealType cusp_ud=-0.5*q*q;
    FuncType *funcUD=new FuncType(cusp_ud,B);
    funcUD->setIDs("pade_ud",ID_B);//set the ID's
    J2->addFunc("pade_uu",0,0,funcUU);
    J2->addFunc("pade_ud",0,1,funcUD);
    FuncList.push_back(funcUU);
    FuncList.push_back(funcUD);
    //DerivFuncType *dfuncUU=new DerivFuncType(cusp_uu,B);
    //DerivFuncType *dfuncUD=new DerivFuncType(cusp_ud,B);
    //dJ2->addFunc("pade_uu",0,0,dfuncUU);
    //dJ2->addFunc("pade_ud",0,1,dfuncUD);
    //dFuncList.push_back(dfuncUU);
    //dFuncList.push_back(dfuncUD);
    app_log() << "    Adding Spin-dependent Pade Two-Body Jastrow " << "\n";
    app_log() << "      parallel spin     " << cusp_uu << "\n";
    app_log() << "      antiparallel spin " << cusp_ud << "\n";
  }
//#if defined(ENABLE_SMARTPOINTER)
//    if(!dPsi.get())
//      dPsi = DiffOrbitalBasePtr(new AnalyticDiffOrbital(targetPtcl.getTotalNum(),dJ2));
//#else
//    if(dPsi==0)
//      dPsi=new AnalyticDiffOrbital(targetPtcl.getTotalNum(),dJ2);
//#endif
//    else
//      dPsi->addOrbital(dJ2);
//
//    dPsi->resetTargetParticleSet(targetPtcl);
  app_log() << "  PadeConstraints:: B = " << B <<"\n";
  return J2;
}

OrbitalBase* PadeConstraints::createOneBody(ParticleSet& source)
{
  app_log() << "  PadeBuilder::Adding Pade One-Body Jastrow with effective ionic charges." << std::endl;
  typedef OneBodyJastrowOrbital<FuncType> JneType;
  JneType* J1 = new JneType(source,targetPtcl);
  //typedef OneBodyJastrowOrbital<DerivFuncType> DerivJneType;
  //DerivJneType* dJ1=new DerivJneType(source,targetPtcl);
  SpeciesSet& Species(source.getSpeciesSet());
  int ng=Species.getTotalNum();
  int icharge = Species.addAttribute("charge");
  for(int ig=0; ig<ng; ig++)
  {
    RealType zeff=Species(icharge,ig);
    std::ostringstream j1id;
    j1id<<"pade_"<<Species.speciesName[ig];
    RealType sc=std::pow(2*zeff,0.25);
    FuncType *func=new FuncType(-zeff,B,sc);
    func->setIDs(j1id.str(),ID_B);
    J1->addFunc(ig,func);
    FuncList.push_back(func);
    //DerivFuncType *dfunc=new DerivFuncType(-zeff,B,sc);
    //dJ1->addFunc(ig,dfunc);
    //dFuncList.push_back(dfunc);
    app_log() << "    " << Species.speciesName[ig] <<  " Zeff = " << zeff << " B= " << B*sc << std::endl;
  }
  //add dJ1 to dPsi
  //dPsi->addOrbital(dJ1);
  return J1;
}

////////////////////////////////////////
//ScaledPadeConstraints definitions
////////////////////////////////////////
ScaledPadeConstraints::ScaledPadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
  OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
{
  JComponent.set(TWOBODY);
}
ScaledPadeConstraints::~ScaledPadeConstraints()
{
  delete_iter(FuncList.begin(), FuncList.end());
}

void ScaledPadeConstraints::resetParameters(const opt_variables_type& active)
{
  APP_ABORT("ScaledPadeConstraints::resetParameters is broken. Fix it!");
  // bool update=false;
  // OptimizableSetType::iterator it(optVariables.find(ID_B));
  // if(it != optVariables.end())
  // {
  //   B=(*it).second;
  //   update=true;
  // }
  // OptimizableSetType::iterator it_c(optVariables.find(ID_C));
  // if(it_c != optVariables.end())
  // {
  //   C=(*it_c).second;
  //   update=true;
  // }
  // if(update)
  //   for(int i=0; i<FuncList.size(); i++) {
  //     FuncList[i]->B=B;
  //     FuncList[i]->C=C;
  //     FuncList[i]->resetParameters(optVariables);
  //   }
}


bool ScaledPadeConstraints::put(xmlNodePtr cur)
{
  bool success=getVariables(cur);
  std::map<std::string,std::pair<std::string,RealType> >::iterator bit(inVars.find("B"));
  std::map<std::string,std::pair<std::string,RealType> >::iterator cit(inVars.find("C"));
  if(bit == inVars.end() || cit == inVars.end())
    return false;
  ID_B=(*bit).second.first;
  B=(*bit).second.second;
  ID_C=(*cit).second.first;
  C=(*cit).second.second;
  return true;
}

OrbitalBase* ScaledPadeConstraints::createTwoBody()
{
  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
  JeeType *J2 = new JeeType(targetPtcl);
  if(IgnoreSpin)
  {
    app_log() << "  ScaledPadeConstraints::Adding Spin-independent Pade Two-Body Jastrow " << std::endl;
    FuncType *func=new FuncType(-0.5,B,C);
    J2->addFunc("pade_uu",0,0,func);
    //dJ2->addFunc("pade_uu",0,0,dfunc);
    FuncList.push_back(func);
  }
  else
  {
    app_log() << "  ScaledPadeConstraints::Adding Spin-dependent Pade Two-Body Jastrow " << std::endl;
    FuncType *funcUU=new FuncType(-0.25,B,C);
    FuncType *funcUD=new FuncType(-0.5,B,C);
    J2->addFunc("pade_uu",0,0,funcUU);
    J2->addFunc("pade_ud",0,1,funcUD);
    FuncList.push_back(funcUU);
    FuncList.push_back(funcUD);
  }
  app_log() << "  ScaledPadeConstraints:: B = " << B << " C = " << C << std::endl;
  return J2;
}

OrbitalBase*
ScaledPadeConstraints::createOneBody(ParticleSet& source)
{
  //return 0 for now
  return 0;
}

//////////////////////////////////////////
////PadeOnGridConstraints definitions
//////////////////////////////////////////
//PadeOnGridConstraints::~PadeOnGridConstraints() {
//  delete_iter(FuncList.begin(), FuncList.end());
//}

//bool PadeOnGridConstraints::put(xmlNodePtr cur) {
//  bool success=getVariables(cur);
//  std::map<std::string,std::pair<std::string,RealType> >::iterator vit(inVars.find("B"));
//  if(vit == inVars.end()) return false; //disaster, need to abort
//  ID=(*vit).second.first; B=(*vit).second.second;
//  return true;
//}

//void PadeOnGridConstraints::apply() {
//  for(int i=0; i<FuncList.size(); i++) {
//    InFuncList[i]->B0=B;
//    FuncList[i]->reset();
//  }
//}

//void PadeOnGridConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
//  //potentially add Rcut
//  outVars.add(ID,&B,1);
//}

//OrbitalBase* PadeOnGridConstraints::createTwoBody() {

//  setRadialGrid(targetPtcl);

//  typedef FuncType::FNOUT OutFuncType;
//  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
//  JeeType *J2 = new JeeType(target);
//  if(IgnoreSpin) {
//    app_log() << "  PadeOnGridConstraints::Adding Spin-independent Pade Two-Body Jastrow B=" << B << std::endl;
//    //create an analytic input functor
//    InFuncType *infunc=new InFuncType(-0.5,B);
//    //create a numerical functor
//    FuncType* nfunc= new FuncType;
//    //initialize the numerical functor
//    nfunc->initialize(infunc,myGrid,Rcut);

//    InFuncList.push_back(infunc);
//    FuncList.push_back(nfunc);
//    for(int i=0; i<4; i++) J2->addFunc(nfunc);
//  } else {
//    app_log() << "  PadeOnGridConstraints::Adding Spin-dependent Pade Two-Body Jastrow B= " << B << std::endl;
//    InFuncType *uu=new InFuncType(-0.25,B);
//    InFuncType *ud=new InFuncType(-0.5,B);

//    FuncType *funcUU= new FuncType;
//    funcUU->initialize(uu,myGrid,Rcut);

//    FuncType *funcUD= new FuncType;
//    funcUD->initialize(ud,myGrid,Rcut);

//    InFuncList.push_back(uu);
//    InFuncList.push_back(ud);

//    FuncList.push_back(funcUU);
//    FuncList.push_back(funcUD);

//    J2->addFunc(funcUU);//uu
//    J2->addFunc(funcUD);//ud
//    J2->addFunc(funcUD);//du
//    J2->addFunc(funcUU);//dd
//  }
//  return J2;
//}

//OrbitalBase* PadeOnGridConstraints::createOneBody(ParticleSet& source) {
//  //return 0 for the moment
//  return 0;
//}

}
