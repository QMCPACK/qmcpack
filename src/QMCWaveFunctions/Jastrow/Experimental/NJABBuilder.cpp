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
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/NJABBuilder.h"
#include "QMCWaveFunctions/Jastrow/PadeJastrow.h"
#include "QMCWaveFunctions/Jastrow/NoCuspJastrow.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
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
NJABBuilder::NJABBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(p,psi), ptclPool(psets), gridPtr(0), sourcePtcl(0)
{ }

NJABBuilder::InFuncType*
NJABBuilder::createInFunc(const std::string& jastfunction)
{
  if(jastfunction == "nocusp")
  {
    return new NoCuspJastrow<RealType>;
  }
  else
    if(jastfunction == "pade")
    {
      return new PadeJastrow<RealType>;
    }
    else
      if(jastfunction == "pade2")
      {
        return new PadeJastrow2<RealType>;
      }
  return 0;
}

/** create Input Analytic function
 *
 * \xmlonly
 * <jastrow name="Jne"
 *   type="Two-Body|One-Body|Polarization|Three-Body-Geminal"
 *   function="pade|pade2|no-cusp"
 *   transform="no|yes" spin="no|yes"
 *   source="ionic system">
 *   <grid/>
 *   <correlation speciesA="sourceSpecies" speciesB="targetSpecies" type="pade|pade2|no-cusp">
 *      <parameter name=" ">value<parameter>
 *   </correlation>
 * </jastrow>
 * \endxmlonly
 */
bool NJABBuilder::putInFunc(xmlNodePtr cur)
{
  std::string corr_tag("correlation");
  std::string jastfunction("pade");
  int	ng=1;
  const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
  if(ftype != NULL)
    jastfunction = (const char*) ftype;
  const xmlChar* s=xmlGetProp(cur,(const xmlChar*)"source");
  if(s != NULL)
  {
    std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find((const char*)s));
    if(pa_it == ptclPool.end())
      return false;
    sourcePtcl = (*pa_it).second;
    ng=sourcePtcl->getSpeciesSet().getTotalNum();
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
      if(cname == dtable_tag)
      {
        std::string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
        std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find(source_name));
        if(pa_it == ptclPool.end())
          return false;
        sourcePtcl=(*pa_it).second;
        ng = sourcePtcl->getSpeciesSet().getTotalNum();
        XMLReport("Number of sources " << ng)
        InFunc.resize(ng,0);
      }
      else
        if(cname ==corr_tag)
        {
          if(sourcePtcl==0)
            return false;
          std::string jfunctype(jastfunction);
          std::string spA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
          ftype = xmlGetProp(cur, (const xmlChar *)"type");
          if(ftype)
          {
            jfunctype=(const char*)ftype;
          }
          ia = sourcePtcl->getSpeciesSet().findSpecies(spA);
          if(!(InFunc[ia]))
          {
            InFuncType *j1=createInFunc(jfunctype);
            InFunc[ia]= j1;
            app_log() <<"   Added Jastrow Correlation ("<<jfunctype
                      << ") between " <<spA<<" and "<<targetPtcl.getName() << std::endl;
          }
          InFunc[ia]->put(cur);
          InFunc[ia]->addOptimizables(targetPsi.VarList);
        }
    cur = cur->next;
  } // while cur
  return true;
}

bool NJABBuilder::put(xmlNodePtr cur)
{
  //create analytic functions
  bool success = putInFunc(cur);
  //create grid and initialize CubicSplineFunctions
  OneDimGridFactory::GridType* agrid = OneDimGridFactory::createGrid(gridPtr);
  //get the cutoff value to smooth the function
  RealType rcut = OneDimGridFactory::setSmoothCutoff(agrid,gridPtr);
  OneBodyJastrow<FuncType> *J1 = new OneBodyJastrow<FuncType>(*sourcePtcl, targetPtcl);
  for(int i=0; i<InFunc.size(); i++)
  {
    if(InFunc[i])
    {
      FuncType* ofunc= new FuncType;
      ofunc->setInFunc(InFunc[i]);
      ofunc->setOutFunc(new OutFuncType(agrid));
      ofunc->setCutoff(rcut,agrid->rmax());
      ofunc->reset();
      J1->addFunc(i,ofunc);
    }
    else
    {
      J1->addFunc(i,0);
    }
  }
  J1->setOptimizable(true);
  targetPsi.addOrbital(J1);
  XMLReport("Added a One-Body Jastrow Function")
  return success;
}
}
