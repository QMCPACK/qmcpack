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
    
    



#ifndef QMCPLUSPLUS_SIMPLE_UTILITIES_H
#define QMCPLUSPLUS_SIMPLE_UTILITIES_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{
/** A convenient function to evaluate the electron configuration
 *@param el the electron configuration
 *@param acontext xpath context
 */
inline
bool
determineNumOfElectrons(ParticleSet& el, xmlXPathContextPtr acontext)
{
  //initialize el with the wave function information
  //This is to determine the number of up and down electrons
  std::vector<int> N;
  std::string sdname("//");
  sdname.append(OrbitalBuilderBase::sd_tag);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)(sdname.c_str()),acontext);
  bool found_el=false;
  int nsd= result->nodesetval->nodeNr;
  XMLReport("Found " << nsd << " SlaterDeterminant.")
  if(nsd)
  {
    std::vector<xmlNodePtr> dset;
    xmlNodePtr cur=result->nodesetval->nodeTab[0]->children;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == OrbitalBuilderBase::det_tag)
        dset.push_back(cur);
      cur=cur->next;
    }
    if(dset.size())
    {
      XMLReport("Found " << dset.size() << OrbitalBuilderBase::det_tag)
      for(int i=0; i<dset.size(); i++)
      {
        int norb = 0;
        xmlChar* orb=xmlGetProp(dset[i],(const xmlChar*)"orbitals");
        if(orb)
        {
          norb = atoi((const char*)orb);
          XMLReport("Using attribute orbitals " << norb)
        }
        else
        {
          cur = dset[i]->children;
          while(cur != NULL)
          {
            std::string cname((const char*)(cur->name));
            if(cname == OrbitalBuilderBase::spo_tag)
              norb++;
            cur=cur->next;
          }
          XMLReport("Counting the number or ortbials " << norb)
        }
        N.push_back(norb);
      }
      el.create(N);
      found_el=true;
    }
  }
  xmlXPathFreeObject(result);
  if(N.empty())
  {
    result = xmlXPathEvalExpression((const xmlChar*)"//particleset",acontext);
    int n= result->nodesetval->nodeNr;
    int pset=0;
    while(!found_el && pset<n)
    {
      xmlChar* s= xmlGetProp(result->nodesetval->nodeTab[pset],(const xmlChar*)"name");
      if(s)
      {
        if(s[0] == 'e')
        {
          xmlNodePtr cur = result->nodesetval->nodeTab[pset]->children;
          found_el=true;
          while(cur != NULL)
          {
            std::string cname((const char*)(cur->name));
            if(cname == "group")
            {
              xmlChar* g= xmlGetProp(cur,(const xmlChar*)"size");
              if(s)
              {
                N.push_back(atoi((const char*)g));
              }
            }
            cur = cur->next;
          }
        }
      }
      pset++;
    }
    xmlXPathFreeObject(result);
  }
  if(N.empty())
  {
    ERRORMSG("Could not determine the number of electrons. Assuming He(1,1)")
    N.resize(2);
    N[0]=1;
    N[1]=1;
  }
  else
  {
    XMLReport("The electron configured by //slaterdeterminant or //particle/group/@name=\'e\'")
  }
  el.create(N);
  return true;
}
}


#endif


