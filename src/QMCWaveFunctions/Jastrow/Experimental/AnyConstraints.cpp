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
    
    
#include "QMCWaveFunctions/Jastrow/AnyConstraints.h"
#include "QMCWaveFunctions/Jastrow/LinearCombinationFunctor.h"
#include "QMCWaveFunctions/Jastrow/CompositeFunctor.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
AnyConstraints::AnyConstraints(ParticleSet& p, TrialWaveFunction& psi):
  OrbitalConstraintsBase(p,psi)
{
  ClassName="AnyConstraints";
  //note MULTIPLE is false
  JComponent.set(ONEBODY);
  JComponent.set(TWOBODY);
}

AnyConstraints::~AnyConstraints()
{
  //delete_iter(FuncList.begin(), FuncList.end());
}

bool AnyConstraints::put(xmlNodePtr cur)
{
  //get generic parameters and grid information
  bool success=getVariables(cur);
  return success;
}

void AnyConstraints::resetParameters(const opt_variables_type& active)
{
  //map<std::string,BasisGroupType*>::iterator it(BasisGroups.begin()),it_end(BasisGroups.end());
  //while(it != it_end)
  //{
  //  (*it).second->Out_->resetParameters(active);
  //  ++it;
  //}
}

//void AnyConstraints::addOptimizables(OptimizableSetType& outVars)
//{
//  //createBasis will take care of this. Do not add again
//  std::map<std::string,BasisGroupType*>::iterator it(BasisGroups.begin()),it_end(BasisGroups.end());
//  while(it != it_end)
//  {
//    (*it).second->In_->addOptimizables(outVars);
//    ++it;
//  }
//}

void AnyConstraints::addSingleBasisPerSpecies(xmlNodePtr cur)
{
  RealType rcut=10.0;
  int npts=101;
  RealType step=-1.0;
  if(myGrid)
  {
    rcut=myGrid->rmax();
    npts = myGrid->size();
  }
  OhmmsAttributeSet gAttrib;
  gAttrib.add(rcut,"rf");
  gAttrib.add(npts,"npts");
  gAttrib.add(step,"step");
  std::string tpname(targetPtcl.getName());
  BasisGroupType* curBG=0;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    std::string elementType("e");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(elementType,"elementType");
    aAttrib.put(cur);
    if(cname == "atomicBasisSet")
    {
      //replace elementType for clones
      if(tpname.find(elementType)<tpname.size())
        elementType=tpname;
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)(cur1->name));
        if(cname1 == "basisGroup")
        {
          curBG=createBasisGroup(cur1,elementType);
          curBG->setGrid(rcut,npts);
          add2BasisGroup(curBG,cur1);
        }
        else
          if(cname1 == "grid")
            gAttrib.put(cur1);
        cur1=cur1->next;
      }
    }
    else
      if(cname == "basisGroup")
      {
        //replace elementType for clones
        if(tpname.find(elementType)<tpname.size())
          elementType=tpname;
        curBG=createBasisGroup(cur,elementType);
        curBG->setGrid(rcut,npts);
        add2BasisGroup(curBG,cur);
      }
      else
        if(cname == "grid")
          gAttrib.put(cur);
    cur=cur->next;
  }
}


AnyConstraints::BasisGroupType*
AnyConstraints::createBasisGroup(xmlNodePtr cur, const std::string& elementType)
{
  ReportEngine PRE(ClassName,"createBasisGroup(...)");
  std::string type("WM");
  RealType cusp=0.0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(type,"type");
  aAttrib.add(cusp,"cusp");
  aAttrib.put(cur);
  BGContainerType::iterator it(BasisGroups.find(elementType));
  BasisGroupType* curBG=0;
  if(it == BasisGroups.end())
  {
    curBG=new BasisGroupType;
    BasisGroups[elementType]=curBG;
  }
  else
  {
    curBG=(*it).second;
  }
  //overwrite a cusp
  curBG->Cusp=cusp;
  return curBG;
}

void AnyConstraints::add2BasisGroup(BasisGroupType* curBG, xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"add2BasisGroup(...)");
  typedef LinearCombinationFunctor<RealType> ComboFuncType;
  ComboFuncType* acombo=new ComboFuncType;
  curBG->In_ = acombo;
  //DerivFuncType* aderiv=new DerivFuncType(curBG->Rcut);
  std::string radID("0");
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "radfunc")
    {
      OhmmsAttributeSet rAttrib;
      std::string rfuncName("WM");
      RealType exponent=1.0;
      RealType contraction=1.0;
      RealType rcut(curBG->Rcut);
      int rpower=0;
      rAttrib.add(radID,"id"); //rAttrib.add(a->B0,"b");
      rAttrib.add(exponent,"exponent");
      rAttrib.add(contraction,"contraction");
      //rAttrib.add(rpower,"node");
      rAttrib.add(rcut,"rcut");
      rAttrib.add(rfuncName,"type");
      rAttrib.put(cur);
      acombo->cutoff_radius=rcut;
      OptimizableFunctorBase *a=0;
      OptimizableFunctorBase *da=0;
      if(rfuncName == "cusp")
      {
        //a= new CuspCorrectionFunctor<RealType>(exponent,rcut);
        curBG->Cusp=contraction;
        rpower=0;//overwrite the power
      }
      else
      {
        a= new WMFunctor<RealType>(exponent,rcut);
        a->put(cur);
        std::string id_c=radID+"_C";
        acombo->addComponent(a,contraction,id_c);
        //add a component to the derivative
        //aderiv->addComponent(contraction,exponent,radID);
      }
      //else
      //{//this is useless
      //  AnyTimesRnFunctor<RealType>* awrap=new AnyTimesRnFunctor<RealType>(a,rpower);
      //  acombo->addComponent(awrap,contraction,id_c,cur);
      //}
      app_log()  << "<radfunc id=\"" << radID << "\" exponent=\""<< exponent
                 << "\" contraction=\"" << contraction
                 << "\" node=\"" << rpower << "\" rcut=\"" << rcut << "\"/>" << std::endl;
    }
    cur=cur->next;
  }
  //add variables to the optimization list
  //aderiv->addOptimizables(targetPsi.VarList);
  //for(int i=0; i<targetPsi.VarList.size(); ++i)
  //{
  //  std::cout << targetPsi.VarList.getName(i) << " " << targetPsi.VarList.getValue(i) << std::endl;
  //}
  //non-zero cusp
  if(std::abs(curBG->Cusp)>std::numeric_limits<RealType>::epsilon())
  {
    CuspCorrectionFunctor<RealType> *a  = new CuspCorrectionFunctor<RealType>(2.0,curBG->Rcut);
    app_log() << "  Adding a cusp term: " << curBG->Cusp << "* (-1/b exp(-br)), b=" << a->E << std::endl;
    std::string cusp_tag=radID+"_cusp"; //non optimizable
    acombo->addComponent(a,curBG->Cusp,cusp_tag,true);//this cannot be modified
  }
  //curBG->Deriv_ = aderiv;
  //add optimizable values now
  //aderiv->addOptimizables(targetPsi.VarList);
}

OrbitalBase* AnyConstraints::createTwoBody()
{
  xmlNodePtr cur=myNode->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "basisset")
    {
      addSingleBasisPerSpecies(cur);
    }
    cur=cur->next;
  }
  if(BasisGroups.empty())
  {
    app_error() << "  AnyConstraints::createTwoBody fails to create a TwoBodyJastrow "
                << " due to missing <basisset/> " << std::endl;
    return 0;
  }
  BasisGroupType* curGroup=0;
  BGContainerType::iterator it(BasisGroups.find(targetPtcl.getName()));
  if(it == BasisGroups.end())
  {
    return 0;
  }
  else
  {
    curGroup=(*it).second;
  }
  InFuncType* infunc=curGroup->In_ ;
  OutFuncType* nfunc= new OutFuncType(infunc, curGroup->Rcut, curGroup->NumGridPoints);
  curGroup->Out_ = nfunc;
#if !defined(HAVE_MPI)
  if(PrintTables)
  {
    std::ofstream fout("J2.dat");
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout << "# Two-body Jastrow generated by AnyContraints::createTwoBody" << std::endl;
    nfunc->print(fout);
  }
#endif
  //create a Jastrow function
  typedef TwoBodyJastrowOrbital<OutFuncType> JeeType;
  JeeType *J2 = new JeeType(targetPtcl);
  J2->addFunc("Jee",0,0,nfunc);
  //2008-04-07 derivatives not complete
  //typedef DiffTwoBodyJastrowOrbital<DerivFuncType> dJeeType;
  //dJeeType *dJ2 = new dJeeType(targetPtcl);
  //dJ2->addFunc("Jee",0,0,curGroup->Deriv_);
  //dJ2->initialize();
  ////add a derivative function
  //J2->setDiffOrbital(dJ2);
  return J2;
}

OrbitalBase* AnyConstraints::createOneBody(ParticleSet& source)
{
  std::map<std::string,InFuncType*> jnSet;
  xmlNodePtr cur=myNode->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "basisset")
    {
      addSingleBasisPerSpecies(cur);
    }
    cur=cur->next;
  }
  int nSpecies = source.getSpeciesSet().getTotalNum();
  typedef OneBodyJastrowOrbital<OutFuncType> JneType;
  JneType* jne=new JneType(source,targetPtcl);
  //typedef DiffOneBodyJastrowOrbital<DerivFuncType> dJneType;
  //dJneType *djne = new dJneType(source,targetPtcl);
  bool foundit=false;
  BGContainerType::iterator jit(BasisGroups.begin()), jit_end(BasisGroups.end());
  while(jit != jit_end)
  {
    int ig=source.getSpeciesSet().findSpecies((*jit).first);
    if(ig < nSpecies) //should not add any species here
    {
      foundit=true;
      BasisGroupType* curG((*jit).second);
      OutFuncType* nfunc= new OutFuncType(curG->In_,curG->Rcut, curG->NumGridPoints);
      jne->addFunc(ig,nfunc);
      curG->Out_ = nfunc;
      //add derivatives
      //djne->addFunc(ig,curG->Deriv_);
#if !defined(HAVE_MPI)
      if(PrintTables)
      {
        char fname[16];
        sprintf(fname,"J1.%s.dat",source.getSpeciesSet().speciesName[ig].c_str());
        std::ofstream fout(fname);
        fout.setf(std::ios::scientific, std::ios::floatfield);
        fout << "# One-body Jastrow " << source.getSpeciesSet().speciesName[ig]
             << " generated by AnyContraints::createOneBody" << std::endl;
        nfunc->print(fout);
      }
#endif
    }
    ++jit;
  }
  //if(foundit)
  //{
  //  djne->initialize();
  //  jne->setDiffOrbital(djne);
  //}
  //else
  //{
  //  delete jne; jne=0;
  //  delete djne;
  //}
  return jne;
}


// void
//   AnyConstraints::createDefaultTwoBody(xmlNodePtr cur, const std::string& tname) {
//    //xmlNodePtr basisNode = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
//    //xmlNewProp(basisNode,(const xmlChar*)"source",(const xmlChar*)tname.c_str());
//    //xmlNodePtr b0 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"4.55682");
//    //xmlNewProp(b0,(const xmlChar*)"id",(const xmlChar*)"eeB0");
//    //xmlNewProp(b0,(const xmlChar*)"name",(const xmlChar*)"B");
//    //xmlNodePtr b1 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"1.83679");
//    //xmlNewProp(b1,(const xmlChar*)"id",(const xmlChar*)"eeB1");
//    //xmlNewProp(b1,(const xmlChar*)"name",(const xmlChar*)"B");
//    //xmlNodePtr b2 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"4.59201");
//    //xmlNewProp(b2,(const xmlChar*)"id",(const xmlChar*)"eeB2");
//    //xmlNewProp(b2,(const xmlChar*)"name",(const xmlChar*)"B");
//    //xmlNodePtr b3 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0.226152");
//    //xmlNewProp(b3,(const xmlChar*)"id",(const xmlChar*)"eeB3");
//    //xmlNewProp(b3,(const xmlChar*)"name",(const xmlChar*)"B");
//    //xmlAddChild(cur,basisNode);

//    //addBasisGroup(basisNode);

//    //xmlNodePtr corrNode = xmlNewNode(NULL,(const xmlChar*)"correlation");
//    //xmlNewProp(corrNode,(const xmlChar*)"speciesA",(const xmlChar*)tname.c_str());
//    //xmlNewProp(corrNode,(const xmlChar*)"speciesB",(const xmlChar*)tname.c_str());
//    //xmlNodePtr c0 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"-0.137958");
//    //xmlNewProp(c0,(const xmlChar*)"id",(const xmlChar*)"eeC0");
//    //xmlNewProp(c0,(const xmlChar*)"name",(const xmlChar*)"C");
//    //xmlNodePtr c1 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"-1.09405");
//    //xmlNewProp(c1,(const xmlChar*)"id",(const xmlChar*)"eeC1");
//    //xmlNewProp(c1,(const xmlChar*)"name",(const xmlChar*)"C");
//    //xmlNodePtr c2 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"1.06578");
//    //xmlNewProp(c2,(const xmlChar*)"id",(const xmlChar*)"eeC2");
//    //xmlNewProp(c2,(const xmlChar*)"name",(const xmlChar*)"C");
//    //xmlNodePtr c3 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0.96602");
//    //xmlNewProp(c3,(const xmlChar*)"id",(const xmlChar*)"eeC3");
//    //xmlNewProp(c3,(const xmlChar*)"name",(const xmlChar*)"C");
//    //xmlAddChild(cur,corrNode);

//    //map<std::string,BasisSetType*>::iterator it(myBasisSet.find("e"));
//    //return createCorrelation(corrNode,(*it).second);
// }
}
