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
    
    
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/Jastrow/NJAABuilder.h"
#include "QMCWaveFunctions/Jastrow/PadeJastrow.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCFactory/OneDimGridFactory.h"

namespace qmcplusplus
{


/** constructor
 * @param p target ParticleSet whose wave function is to be initialized
 *@param psi the wavefunction
 *@param psets a vector containing the ParticleSets
 *
 *Jastrow wavefunctions chose the needed distance tables and the
 *DistanceTableData objects are initialized based on the source
 *and target particle sets.
 */
NJAABuilder::NJAABuilder(ParticleSet& p, TrialWaveFunction& psi):
  OrbitalBuilderBase(p,psi),gridPtr(NULL),IgnoreSpin(true)
{ }

NJAABuilder::InFuncType*
NJAABuilder::createInFunc(const std::string& jastfunction, int speciesA, int speciesB)
{
  //return new PadeJastrow<RealType>;
  if(jastfunction == "rpa")
  {
    RPAJastrow<RealType>* newRPA=0;
    if(IgnoreSpin)
      newRPA= new RPAJastrow<RealType>(false);
    else
      newRPA = new RPAJastrow<RealType>(speciesA==speciesB);
    if(targetPtcl.Lattice.BoxBConds[0])
      newRPA->setDensity(targetPtcl.getTotalNum()/targetPtcl.Lattice.Volume);
    return newRPA;
  }
  return new PadeJastrow<RealType>;
}
bool NJAABuilder::putInFunc(xmlNodePtr cur)
{
  std::string corr_tag("correlation");
  std::string jastfunction("pade");
  const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
  if(ftype)
    jastfunction = (const char*) ftype;
  int	ng = targetPtcl.groups();
  ///happends only once
  if(InFunc.size() != ng*ng)
  {
    InFunc.resize(ng*ng,0);
  }
  int ia=0, ib=0, iab=0;
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "grid")
    {
      gridPtr=cur; //save the pointer
    }
    else
      if(cname ==corr_tag)
      {
        std::string spA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
        std::string spB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
        const xmlChar* refptr=xmlGetProp(cur,(const xmlChar *)"ref");
        const xmlChar* idptr=xmlGetProp(cur,(const xmlChar *)"id");
        if(!IgnoreSpin)
          //get the species
        {
          ia = targetPtcl.getSpeciesSet().findSpecies(spA);
          ib = targetPtcl.getSpeciesSet().findSpecies(spB);
          iab = ia*ng+ib;
        }
        if(!(InFunc[ia*ng+ib]))
        {
          InFuncType *j2=createInFunc(jastfunction,ia,ib);
          j2->put(cur);
          j2->addOptimizables(targetPsi.VarList);
          InFunc[ia*ng+ib]= j2;
          XMLReport("Added Jastrow Correlation between "<<spA<<" and "<<spB)
        }
        else
        {
          ERRORMSG("Using an existing Jastrow Correlation "<<spA<<" and "<<spB)
        }
      }
    cur = cur->next;
  } // while cur
  return true;
}

bool NJAABuilder::put(xmlNodePtr cur)
{
  const xmlChar* spin=xmlGetProp(cur,(const xmlChar*)"spin");
  if(spin != NULL)
  {
    std::string a((const char*)spin);
    if(a == "yes")
      IgnoreSpin=false;
  }
  //create analytic functions
  bool success = putInFunc(cur);
  bool insertGrid=false;
  //create a grid node if missing
  //use the lattice constant to set the cutoff but not GENERAL
  if(gridPtr == NULL)
  {
    gridPtr = xmlNewNode(NULL,(const xmlChar*)"grid");
    xmlNewProp(gridPtr,(const xmlChar*)"type",(const xmlChar*)"log");
    xmlNewProp(gridPtr,(const xmlChar*)"ri",(const xmlChar*)"1.0e-6");
    std::ostringstream rf;
    if(targetPtcl.Lattice.BoxBConds[0])
      rf << targetPtcl.Lattice.R(0,0)*0.48;
    else
      rf << 100.;
    xmlNewProp(gridPtr,(const xmlChar*)"rf",(const xmlChar*)rf.str().c_str());
    xmlNewProp(gridPtr,(const xmlChar*)"npts",(const xmlChar*)"101");
    insertGrid=true;
  }
  //create grid and initialize CubicSplineFunctions
  OneDimGridFactory::GridType* agrid = OneDimGridFactory::createGrid(gridPtr);
  //get the cutoff radius
  RealType rcut = OneDimGridFactory::setSmoothCutoff(agrid,gridPtr);
  app_log() << "  smoothing function starts at " << rcut << std::endl;
  DistanceTableData* d_table = DistanceTable::add(targetPtcl);
  TwoBodyJastrowOrbital<FuncType> *J2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl,d_table);
  int	ng = targetPtcl.groups();
  FuncType* dfunc= new FuncType;
  dfunc->setInFunc(InFunc[0]);
  dfunc->setOutFunc(new OutFuncType(agrid));
  dfunc->setCutoff(rcut,agrid->rmax());
  dfunc->reset();
  J2->F.resize(ng*ng,dfunc);
  J2->insert("j00",dfunc);
  if(!IgnoreSpin)
  {
    char oname[8];
    for(int i=0; i<ng-1; i++)
    {
      for(int j=i+1; j<ng; j++)
      {
        FuncType* ofunc= new FuncType;
        ofunc->setInFunc(InFunc[i*ng+j]);
        ofunc->setOutFunc(new OutFuncType(agrid));
        ofunc->setCutoff(rcut,agrid->rmax());
        ofunc->reset();
        J2->F[i*ng+j]=ofunc;
        J2->F[j*ng+i]=ofunc;
        sprintf(oname,"j%d%d",i,j);
        J2->insert(oname,ofunc);
      }
    }
  }
  J2->setOptimizable(true);
  targetPsi.addOrbital(J2);
  XMLReport("Added a Two-Body Jastrow Function")
  if(insertGrid)
    xmlAddChild(cur,gridPtr);
  return success;
}
}
