//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#if QMC_BUILD_LEVEL>2
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/ModPadeFunctor.h"
#include "QMCWaveFunctions/Jastrow/BesselZeroFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/OpenGaussianSlaterFunctor.h"
#endif
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"

namespace qmcplusplus
{

template<class FN>
bool JABBuilder::createJAB(xmlNodePtr cur, const std::string& jname)
{
  std::string corr_tag("correlation");
  std::vector<FN*> jastrow;
  int ng = 0;
  ParticleSet* sourcePtcl=0;
  const xmlChar* s=xmlGetProp(cur,(const xmlChar*)"source");
  if (s != NULL)
  {
    std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find((const char*)s));
    if (pa_it == ptclPool.end())
      return false;
    sourcePtcl = (*pa_it).second;
    ng=sourcePtcl->getSpeciesSet().getTotalNum();
    for (int i=0; i<ng; i++)
      jastrow.push_back(0);
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == dtable_tag)
    {
      std::string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
      std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find(source_name));
      if (pa_it == ptclPool.end())
        return false;
      sourcePtcl = (*pa_it).second;
      ng=sourcePtcl->getSpeciesSet().getTotalNum();
      XMLReport("Number of sources " << ng)
      for (int i=0; i<ng; i++)
        jastrow.push_back(0);
    }
    else
      if (cname == corr_tag)
      {
        if (sourcePtcl == 0)
          return false;
        std::string spA;
        OhmmsAttributeSet rAttrib;
        rAttrib.add(spA, "speciesA");
        rAttrib.add(spA, "elementType");
        rAttrib.put(cur);
        //string speciesB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
        int ia = sourcePtcl->getSpeciesSet().findSpecies(spA);
        if (!(jastrow[ia]))
        {
          jastrow[ia]= new FN;
          jastrow[ia]->put(cur);
          LOGMSG("  Added Jastrow Correlation between " << spA)
        }
      }
    cur = cur->next;
  } // while cur
  if (sourcePtcl == 0) //invalid input: cleanup
  {
    for (int ig=0; ig<ng; ig++)
      if (jastrow[ig])
        delete jastrow[ig];
    return false;
  }
  typedef OneBodyJastrowOrbital<FN> JneType;
  JneType* J1 = new JneType(*sourcePtcl,targetPtcl);
  typedef DiffOneBodyJastrowOrbital<FN> dJ1Type;
  dJ1Type *dJ1 = new dJ1Type(*sourcePtcl, targetPtcl);
  for (int ig=0; ig<ng; ig++)
  {
    J1->addFunc(ig,jastrow[ig]);
    dJ1->addFunc(ig,jastrow[ig]);
  }
  J1->dPsi=dJ1;
  std::string j1name="J1_"+jname;
  targetPsi.addOrbital(J1,j1name);
  XMLReport("Added a One-Body Jastrow Function")
  return true;
}

bool JABBuilder::put(xmlNodePtr cur)
{
  std::string jastfunction("pade");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(jastfunction,"function");
  aAttrib.put(cur);
  bool success=false;
  app_log() << "  One-Body Jastrow Function = " << jastfunction << std::endl;
  if (jastfunction == "pade")
  {
    success = createJAB<PadeFunctor<RealType> >(cur,jastfunction);
  }
  else
    if (jastfunction == "pade2")
    {
      success = createJAB<Pade2ndOrderFunctor<RealType> >(cur,jastfunction);
    }
#if QMC_BUILD_LEVEL>2
    else
      if (jastfunction == "short")
      {
        success = createJAB<ModPadeFunctor<RealType> >(cur,jastfunction);
      }
      else
        if (jastfunction == "gaussian")
        {
          success = createJAB<GaussianFunctor<RealType> >(cur,jastfunction);
        }
        else
          if (jastfunction == "opengaussianslater")
          {
            app_log()<<jastfunction<< std::endl;
            success = createJAB<OpenGaussianSlaterFunctor<RealType> >(cur,jastfunction);
          }
          else
            if (jastfunction == "besselzero")
            {
              success = createJAB<BesselZero<RealType> >(cur,jastfunction);
            }
            else
              if (jastfunction == "shiftedgaussian")
              {
                success = createJAB<TruncatedShiftedGaussianFunctor<RealType> >(cur,jastfunction);
              }
              else
                if (jastfunction == "bspline")
                {
                  success = createJAB<BsplineFunctor<RealType> >(cur,jastfunction);
                }
#endif
                else
                {
                  app_error() << "Unknown one body function: " << jastfunction << ".\n";
                }
  return success;
}
}
