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
    
    
#include "QMCWaveFunctions/Jastrow/PolyConstraints.h"
#include "QMCWaveFunctions/Jastrow/LinearCombinationFunctor.h"
#include "QMCWaveFunctions/Jastrow/CompositeFunctor.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
PolyConstraints::PolyConstraints(ParticleSet& p, TrialWaveFunction& psi, bool spinInd):
  OrbitalConstraintsBase(p,psi),IgnoreSpin(spinInd)
{
  //MULTIPLE is false
  JComponent.set(ONEBODY);
  JComponent.set(TWOBODY);
}

PolyConstraints::~PolyConstraints()
{
}

bool PolyConstraints::put(xmlNodePtr cur)
{
  //get generic parameters and grid information
  bool success=getVariables(cur);
  return success;
}

void PolyConstraints::resetParameters(OptimizableSetType& optVariables)
{
  std::map<std::string,BasisGroupType*>::iterator it(BasisGroups.begin()),it_end(BasisGroups.end());
  while(it != it_end)
  {
    (*it).second->resetParameters(optVariables);
    ++it;
  }
}

void PolyConstraints::addOptimizables(OptimizableSetType& outVars)
{
  std::map<std::string,BasisGroupType*>::iterator it(BasisGroups.begin()),it_end(BasisGroups.end());
  while(it != it_end)
  {
    (*it).second->addOptimizables(outVars);
    ++it;
  }
}

void PolyConstraints::addSingleBasisPerSpecies(xmlNodePtr cur)
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
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)(cur1->name));
        if(cname1 == "basisGroup")
        {
          createBasisGroup(cur1,elementType,rcut);
        }
        else
          if(cname1 == "grid")
          {
            gAttrib.put(cur1);
          }
        cur1=cur1->next;
      }
    }
    else
      if(cname == "basisGroup")
      {
        createBasisGroup(cur,elementType,rcut);
      }
      else
        if(cname == "grid")
          gAttrib.put(cur);
    cur=cur->next;
  }
}


void
PolyConstraints::createBasisGroup(xmlNodePtr cur, const std::string& elementType, RealType rcut)
{
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
  curBG->setL(rcut);
  curBG->put(cur);
}

OrbitalBase* PolyConstraints::createTwoBody()
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
    app_error() << "  PolyConstraints::createTwoBody fails to create a TwoBodyJastrow "
                << " due to missing <basisset/> " << std::endl;
    return 0;
  }
  //create a Jastrow function: spin-independent for now
  typedef TwoBodyJastrowOrbital<BasisGroupType> JeeType;
  JeeType *J2 = new JeeType(targetPtcl);
  BGContainerType::iterator it(BasisGroups.find(targetPtcl.getName()));
  if(it == BasisGroups.end())  //spin-dependent
  {
    delete J2;
    J2=0;
  }
  else  //spin-independent
  {
    BasisGroupType* curGroup=(*it).second;
    J2->insert("Jee",curGroup);
    for(int i=0; i<4; i++)
      J2->addFunc(curGroup);
  }
  return J2;
}

OrbitalBase* PolyConstraints::createOneBody(ParticleSet& source)
{
  std::map<std::string,BasisGroupType*> jnSet;
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
  typedef OneBodyJastrowOrbital<BasisGroupType> JneType;
  JneType* jne=new JneType(source,targetPtcl);
  BGContainerType::iterator jit(BasisGroups.begin()), jit_end(BasisGroups.end());
  while(jit != jit_end)
  {
    int ig=source.getSpeciesSet().findSpecies((*jit).first);
    if(ig < nSpecies) //should not add any species here
    {
      jne->addFunc(ig,(*jit).second);
    }
    ++jit;
  }
  return jne;
}

}
