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
#include "OhmmsData/ParameterSet.h"
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
  double ri = 0.0;
  double rf = 20.0;
  int npts = 100;
  ParameterSet m_param;
  m_param.add(ri,"ri","AU");
  m_param.add(rf,"rf","AU");
  m_param.add(npts,"npts","int");
  std::string type = "elel";
  result
  = xmlXPathEvalExpression((const xmlChar*)"//paircorrelation",m_context);
  xmlNodePtr cur=result->nodesetval->nodeTab[0];
  type = ((const char*)xmlGetProp(cur,(const xmlChar*)"type"));
  XMLReport("Type = " << type)
  xmlNodePtr tcur = cur->children;
  m_param.put(cur);
  std::vector<xmlNodePtr> wset;
  int pid=OHMMS::Controller->mycontext();
  while(tcur != NULL)
  {
    std::string cname((const char*)(tcur->name));
    if(cname == "mcwalkerset")
    {
      int pid_target=pid;
      xmlChar* anode=xmlGetProp(tcur,(const xmlChar*)"node");
      if(anode)
      {
        pid_target = atoi((const char*)anode);
      }
      if(pid_target == pid)
        wset.push_back(tcur);
    }
    tcur=tcur->next;
  }
  std::vector<std::string> ConfigFile;
  if(wset.size())
  {
    for(int i=0; i<wset.size(); i++)
    {
      xmlChar* att= xmlGetProp(wset[i],(const xmlChar*)"file");
      if(att)
        ConfigFile.push_back((const char*)att);
    }
  }
  xmlXPathFreeObject(result);
  int nup = el.last(0);
  int ndown = el.last(1)-el.last(0);
  int factor = 1;
  if(nup)
    factor*=nup;
  if(ndown)
    factor*=ndown;
  DistanceTableData* d_table;
  //create a grid
  LinearGrid<double> grid;
  grid.set(ri,rf,npts);
  double rcut = grid.rmax();
  XMLReport("Grid Values: ri = " << ri << " rf = " << rf << " npts = " << npts)
  int nsp = el.groups();
  int ng = nsp*nsp;
  Matrix<int> histogram;
  histogram.resize(grid.size(),ng);
  d_table = DistanceTable::getTable(DistanceTable::add(el));
  d_table->create(1);
  int count = 0;
  XMLReport("List of the configuration files:")
  for(int i=0; i<ConfigFile.size(); i++)
    XMLReport(ConfigFile[i]);
  for(int i=0; i<ConfigFile.size(); i++)
  {
    HDFWalkerInput wReader(ConfigFile[i]);
    while(wReader.put(el))
    {
      for(MCWalkerConfiguration::iterator it = el.begin();
          it != el.end(); ++it)
      {
        el.R = (*it)->R;
        DistanceTable::update(el);
        for(int i=0; i<d_table->size(SourceIndex); i++)
        {
          for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
          {
            double r = d_table->r(nn);
            if(r < rcut)
            {
              count++;
              int index = grid.index(r);
              int id = d_table->PairID[nn];
              histogram(index,id)++;
            }
          }
        }
      }
    }
  }
  std::string froot(myProject.CurrentRoot());
  froot.append(".pc");
  std::ofstream pc_out(froot.c_str());
  double vol = 4.0*4.0*atan(1.0)/3.0*pow(rcut,3);
  double norm = static_cast<double>(count*factor)/(vol);
  for(int ng=0; ng<histogram.rows()-1; ng++)
  {
    double r1 = grid[ng];
    //double r2 = (i<(grid.size()-1)) ? grid[i+1]:(2.0*grid[i]-grid[i-1]);
    double r2 = grid[ng+1];
    // double r = 0.5*(r1+r2);
    double r = 0.75*(pow(r2,4)-pow(r1,4))/(pow(r2,3)-pow(r1,3));
    pc_out << std::setw(20) << r;
    for(int j=0; j<histogram.cols(); j++)
    {
      //volume element dr = 4/3\pi[(r+dr)^3-r^3]
      //double binVol =  4.0*4.0*atan(1.0)/3.0*(pow(r2,3)-pow(r1,3));
      double binVol = 1.0;
      pc_out << std::setw(20) << static_cast<double>(histogram(ng,j))/(binVol*norm);
    }
    pc_out << std::endl;
  }
  std::cout << "Ionic configuration : " << ion.getName() << std::endl;
  ion.get(std::cout);
  std::cout << "Electronic configuration : " << el.getName() << std::endl;
  el.get(std::cout);
  xmlXPathFreeContext(m_context);
  xmlFreeDoc(m_doc);
  OHMMS::Controller->finalize();
  return 0;
}


