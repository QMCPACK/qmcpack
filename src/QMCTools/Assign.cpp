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
    
    



#include "Configuration.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/OneDimGridBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsApp/ProjectData.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "QMC/QMCUtilities.h"
#include <fstream>

int main(int argc, char **argv)
{
  using namespace qmcplusplus;
  xmlDocPtr m_doc;
  xmlNodePtr m_root;
  xmlXPathContextPtr m_context;
  enum {SourceIndex  = DistanceTableData::SourceIndex};
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  ///project description
  OHMMS::ProjectData myProject;
  ///random number controller
  OHMMS::RandomNumberControl myRandomControl;
  if(argc>1)
  {
    // build an XML tree from a the file;
    LOGMSG("Opening file " << argv[1])
    m_doc = xmlParseFile(argv[1]);
    if (m_doc == NULL)
    {
      ERRORMSG("File " << argv[1] << " is invalid")
      xmlFreeDoc(m_doc);
      return 1;
    }
    // Check the document is of the right kind
    m_root = xmlDocGetRootElement(m_doc);
    if (m_root == NULL)
    {
      ERRORMSG("Empty document");
      xmlFreeDoc(m_doc);
      return 1;
    }
  }
  else
  {
    WARNMSG("No argument is given. Assume that  does not need an input file")
  }
  m_context = xmlXPathNewContext(m_doc);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)"//project",m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    WARNMSG("Project is not defined")
    myProject.reset();
  }
  else
  {
    myProject.put(result->nodesetval->nodeTab[0]);
  }
  xmlXPathFreeObject(result);
  //initialize the random number generator
  xmlNodePtr rptr = myRandomControl.initialize(m_context);
  if(rptr)
  {
    xmlAddChild(m_root,rptr);
  }
  ///the ions
  ParticleSet ion;
  MCWalkerConfiguration el;
  el.setName("e");
  int iu = el.Species.addSpecies("u");
  int id = el.Species.addSpecies("d");
  int icharge = el.Species.addAttribute("charge");
  el.Species(icharge,iu) = -1;
  el.Species(icharge,id) = -1;
  bool init_els = determineNumOfElectrons(el,m_context);
  result
  = xmlXPathEvalExpression((const xmlChar*)"//particleset",m_context);
  xmlNodePtr el_ptr=NULL, ion_ptr=NULL;
  for(int i=0; i<result->nodesetval->nodeNr; i++)
  {
    xmlNodePtr cur=result->nodesetval->nodeTab[i];
    xmlChar* aname= xmlGetProp(cur,(const xmlChar*)"name");
    if(aname)
    {
      char fc = aname[0];
      if(fc == 'e')
      {
        el_ptr=cur;
      }
      else
        if(fc == 'i')
        {
          ion_ptr=cur;
        }
    }
  }
  bool donotresize = false;
  if(init_els)
  {
    el.setName("e");
    XMLReport("The configuration for electrons is already determined by the wave function")
    donotresize = true;
  }
  if(el_ptr)
  {
    XMLParticleParser pread(el,donotresize);
    pread.put(el_ptr);
  }
  if(ion_ptr)
  {
    XMLParticleParser pread(ion);
    pread.put(ion_ptr);
  }
  xmlXPathFreeObject(result);
  if(!ion.getTotalNum())
  {
    ion.setName("i");
    ion.create(1);
    ion.R[0] = 0.0;
  }
  //The ion-ion distance-table
  DistanceTableData* d_ii = DistanceTable::getTable(DistanceTable::add(ion));
  d_ii->create(1);
  d_ii->evaluate(ion);
  std::vector<double> Cut, Core;
  int Centers = ion.getTotalNum();
  //attribute id for cut
  int icut = ion.Species.addAttribute("cut");
  //store the max distance from atom
  Cut.resize(Centers);
  for(int iat=0; iat<Centers; iat++)
  {
    int id = ion.GroupID[iat];
    Cut[iat] = ion.Species(icut,id);
  }
  int icore = ion.Species.addAttribute("core");
  //store the max distance from atom
  Core.resize(Centers);
  for(int iat=0; iat<Centers; iat++)
  {
    Core[iat]=ion.Species(icore,ion.GroupID[iat]);
  }
  //3N-dimensional Gaussian
  ParticleSet::ParticlePos_t chi(el.getTotalNum());
  makeGaussRandom(chi);
  //determine if odd or even number of particles
  int irem = el.getTotalNum()%2;
  int ihalf = el.getTotalNum()/2;
  //assign the core
  int ncore(0);
  for(int iat=0; iat<Centers; iat++)
  {
    double sep=0.8*Cut[iat];
    for(int iel=0; iel<Core[iat]/2; iel++,ncore++)
    {
      el.R[ncore]=ion.R[iat]+sep*chi[ncore];
      el.R[ncore+ihalf]=ion.R[iat]+sep*chi[ncore+ihalf];
    }
  }
  int ipart = ncore;
  int isave_iat=0;
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++)
    {
      double bondlength = d_ii->r(nn);
      int jat = d_ii->J[nn];
      //only assign if the half bond-length < cutoff
      if(bondlength < Cut[iat]+Cut[jat])
      {
        if(ipart < ihalf)
        {
          XMLReport("Assigning particles = " << ipart << " and " << ipart+ihalf)
          /*place 2 electrons (an up and a down) at half
            the bond-length plus a random number multiplied
            by 10% of the bond-length*/
          el.R[ipart] = ion.R[iat]+0.5*d_ii->dr(nn)+0.1*bondlength*chi[ipart];
          el.R[ipart+ihalf] = ion.R[iat]+0.5*d_ii->dr(nn)+0.1*bondlength*chi[ipart+ihalf];
          ipart++;
          isave_iat = iat;
        }
      }
    }
  }
  //assign the last particle (if odd number of particles)
  int flag = 1;
  ipart = el.getTotalNum()-1;
  if(irem)
  {
    XMLReport("Assigning last particle.")
    for(int iat = isave_iat+1; iat<Centers; iat++)
    {
      for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++)
      {
        double bondlength = d_ii->r(nn);
        if((0.5*bondlength < Cut[iat]) && flag)
        {
          XMLReport("Assigning particle = " << ipart)
          el.R[ipart] = ion.R[iat]+0.5*d_ii->dr(nn)+0.1*bondlength*chi[ipart];
          flag = 0;
        }
      }
    }
  }
  std::cout << "Ionic configuration : " << ion.getName() << std::endl;
  ion.get(std::cout);
  std::cout << "Electronic configuration : " << el.getName() << std::endl;
  el.get(std::cout);
  std::string newxml(myProject.CurrentRoot());
  newxml.append(".ptcl.xml");
  std::ofstream ptcl_out(newxml.c_str());
  /*
    std::ofstream molmol("assign.xyz");

    molmol << Centers+el.getTotalNum() << std::endl;
    molmol << std::endl;

    for(int iat=0; iat<Centers; iat++)
    molmol << ion.Species.speciesName[ion.GroupID[iat]] << 0.5292*ion.R[iat] << std::endl;

    for(int ipart=0; ipart<el.getTotalNum(); ipart++)
    molmol << "He" << 0.5292*el.R[ipart] << std::endl;

    molmol.close();
  */
  xmlXPathFreeContext(m_context);
  xmlFreeDoc(m_doc);
  int nup = el.last(0);
  int ndown = el.last(1)-el.last(0);
  ptcl_out << "<?xml version=\"1.0\"?>" << std::endl;
  ptcl_out << "<particleset name=\"e\">" << std::endl;
  ptcl_out << "<group name=\"u\" size=\"" << nup << "\">" << std::endl;
  ptcl_out << "<parameter name=\"charge\">-1</parameter>" << std::endl;
  ptcl_out << "<attrib name=\"position\" datatype=\"posArray\">" << std::endl;
  for (int ipart=0; ipart<nup; ++ipart)
    ptcl_out << el.R[ipart] << std::endl;
  ptcl_out << "</attrib>" << std::endl;
  ptcl_out << "</group>" << std::endl;
  ptcl_out << "<group name=\"d\" size=\"" << ndown << "\">" << std::endl;
  ptcl_out << "<parameter name=\"charge\">-1</parameter>" << std::endl;
  ptcl_out << "<attrib name=\"position\" datatype=\"posArray\">" << std::endl;
  for (int ipart=nup; ipart<el.getTotalNum(); ++ipart)
    ptcl_out << el.R[ipart] << std::endl;
  ptcl_out << "</attrib>" << std::endl;
  ptcl_out << "</group>" << std::endl;
  ptcl_out << "</particleset>" << std::endl;
  OHMMS::Controller->finalize();
  return 0;
}

