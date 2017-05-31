//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Utilities/OhmmsInfo.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

bool MolecularOrbitalBuilder::put(xmlNodePtr cur)
{
  if(xmlHasProp(cur, (const xmlChar*)"href"))
  {
    std::string fname_in = (const char*)(xmlGetProp(cur, (const xmlChar *)"href"));
    LOGMSG("Opening external file " << fname_in)
    putOpen(fname_in);
  }
  else
  {
    putSpecial(cur);
  }
  return true;
}

/** process cur xml node
 *
 * Check attribute list
 * - transform if yes, using numerical grid
 * - source    the particle set providing the nuclei default is "i"
 *\xmlonly
 *  <determinantset type="MolecularOrbital" source="i" transform="yes"/>
 *\endxmlonly
 */
bool MolecularOrbitalBuilder::putSpecial(xmlNodePtr cur)
{
  bool usegrid = true;
  std::string transform("yes"), source("i"), radtype("sto");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(transform,"transform");
  aAttrib.add(transform,"usegrid");
  aAttrib.add(source,"source");
  aAttrib.add(radtype,"keyword");
  aAttrib.put(cur);
  if(transform == "no")
    usegrid=false;
  else
    usegrid=true;
  //ionic system for the molecualr orbitals
  ParticleSet* ions=0;
  //initialize with the source tag
  PtclPoolType::iterator pit(ptclPool.find(source));
  if(pit != ptclPool.end())
  {
    LOGMSG("Molecular orbital uses an ionic system " << source)
    ions=(*pit).second;
  }
  //check if distance table is used
  xmlNodePtr cur1=cur->children;
  while(cur1 != NULL)
  {
    std::string cname1((const char*)cur1->name);
    if(cname1 == "distancetable")
    {
      if(ions == 0)
      {
        const xmlChar*  a=xmlGetProp(cur1,(const xmlChar*)"source");
        if(a)
        {
          source = (const char*)a;
        }
        pit = ptclPool.find(source);
        if(pit != ptclPool.end())
          ions = (*pit).second;
      }
    }
    else
      if(cname1 == basisset_tag)
      {
        xmlNodePtr cur2=cur1->children;
        while(cur2 != NULL)
        {
          if(xmlStrEqual(cur2->name,(const xmlChar*)"atomicBasisSet"))
          {
            const xmlChar* a =xmlGetProp(cur2,(const xmlChar*)"type");
            if(xmlStrEqual(a,(const xmlChar*)"STO"))
            {
              radtype="sto";
            }
            else
            {
              radtype="gto";
            }
          }
          cur2=cur2->next;
        }//end-of-cur2
      }
    cur1=cur1->next;
  }
  if(ions == 0)
  {
    ERRORMSG("Any molecular orbital needs an ionic system. Missing " << source)
    return false;
  }
  if(usegrid)
  {
    LOGMSG("Using radial grids for molecular orbitals for " << targetPtcl.getName() << " with " << ions->getName())
    GridMolecularOrbitals a(targetPtcl,targetPsi,*ions);
    a.put(cur);
  }
  else
  {
    if(radtype == "gto")
    {
      LOGMSG("Using analytic GTO molecular orbitals for " << targetPtcl.getName() << " with " << ions->getName())
      GTOMolecularOrbitals a(targetPtcl,targetPsi,*ions);
      a.put(cur);
    }
    else
    {
      LOGMSG("Using analytic STO molecular orbitals for " << targetPtcl.getName() << " with " << ions->getName())
      STOMolecularOrbitals a(targetPtcl,targetPsi,*ions);
      a.put(cur);
    }
  }
  return true;
}

bool MolecularOrbitalBuilder::putOpen(const std::string& fname_in)
{
  xmlDocPtr doc=NULL;
  xmlNsPtr ns;
  xmlNodePtr cur;
  // build an XML tree from a the file;
  doc = xmlParseFile(fname_in.c_str());
  if (doc == NULL)
  {
    ERRORMSG(fname_in << " does not exist")
    return false;
  }
  ///using XPath instead of recursive search
  xmlXPathContextPtr context;
  // xmlXPathObjectPtr result;
  context = xmlXPathNewContext(doc);
  const char* name = OrbitalBuilderBase::detset_tag.c_str();
  std::string sdname("//");
  sdname.append(OrbitalBuilderBase::detset_tag);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)(sdname.c_str()),context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    ERRORMSG(fname_in << " does not contain any Data!")
  }
  else
  {
    putSpecial(result->nodesetval->nodeTab[0]);
  }
  //free local objects
  xmlXPathFreeObject(result);
  xmlXPathFreeContext(context);
  xmlFreeDoc(doc);
  return true;
}

}
