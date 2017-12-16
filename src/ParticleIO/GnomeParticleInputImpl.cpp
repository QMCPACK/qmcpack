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
    
    


#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include "OhmmsData/FileUtility.h"
#include "Utilities/OhmmsInfo.h"
#include "IO/ParticleIOUtility.h"
#include "Utilities/SpeciesCollection.h"
#include "ParticleBase/ParticleUtility.h"
#if defined(USE_HDF5)
#include "IO/HDFParticleIO.h"
#endif
#include "IO/XMLParticleIO.h"

namespace OHMMS
{

bool LatticeParser::put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur)
{
  double a0 = 1.0;
  cur = cur->xmlChildrenNode;
  std::vector<SingleParticleIndex_t> grid(3,SingleParticleIndex_t(1));
  TinyVector<std::string,OHMMS_DIM> bconds("p");
  while (cur != NULL)
  {
    if(!xmlStrcmp(cur->name, (const xmlChar *)"PARAMETER"))
    {
      std::string aname = (const char*)(xmlGetProp(cur, (const xmlChar *) "name"));
      if(aname == "scale")
      {
        putContent(a0,cur);
      }
      else
        if(aname == "lattice")
        {
          putContent(ref_.R,cur);
        }
        else
          if(aname == "grid")
          {
            putContent(grid[2],cur);
          }
          else
            if(aname == "omgrid")
            {
              putContent(grid[1],cur);
            }
            else
              if(aname == "mpigrid")
              {
                putContent(grid[0],cur);
              }
              else
                if(aname == "bconds")
                {
                  putContent(bconds,cur);
                }
    }
    cur = cur->next;
  }
  for(int idir=0; idir<OHMMS_DIM; idir++)
  {
    char b = bconds[idir][0];
    if(b == 'n' || b == 'N')
    {
      ref_.BConds[idir] = ParticleNoBCond;
      ref_.BConds.wrapper(idir) = ParticleNoBCond;
      ref_.BoxBConds[idir] = false;
    }
  }
  ref_.R *= a0;
  ref_.reset();
  PartitionGrid(ref_,grid);
  return true;
}


bool LatticeXMLWriter::get(std::ostream& os) const
{
  os<< "<UnitCell>" << std::endl;
  os<< "<PARAMETER name=\"lattice\">" << std::endl;
  os << ref_.R << "</PARAMETER>" << std::endl;
  os << "<PARAMETER name=\"bconds\">";
  for(int idir=0; idir<OHMMS_DIM; idir++)
  {
    if(ref_.BoxBConds[idir])
      os << "p ";
    else
      os << "n ";
  }
  os << "</PARAMETER>" << std::endl;
  os << "<PARAMETER name=\"grid\">";
  for(int idir=0; idir<OHMMS_DIM; idir++)
  {
    os <<ref_.Grid[idir] << " ";
  }
  os << "</PARAMETER>" << std::endl;
  os << "</UnitCell>" << std::endl;
  return true;
}


bool ParticleParser::put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur)
{
  if(xmlHasProp(cur, (const xmlChar *) "file"))
  {
    const char* fname
    =(const char*)(xmlGetProp(cur, (const xmlChar *) "file"));
    std::string pformat = getExtension(fname);
    if(pformat == "xml")
    {
      XMLParticleParser aHandle(ref_);
      return  aHandle.put(fname);
    }
    else
    {
      WARNMSG("Using old formats")
      std::ifstream fin(fname);
      if(fin)
      {
        ParticleInputFactory::createParticle(ref_,fin);
      }
      else
      {
        ERRORMSG("File " << fname << "is not found.");
      }
      return true;
    }
  }
  XMLParticleParser aHandle(ref_);
  aHandle.put(doc,ns,cur);
  return false;
}


}
