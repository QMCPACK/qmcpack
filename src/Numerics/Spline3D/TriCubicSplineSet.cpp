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
    
    



#include "Numerics/Spline3D/TriCubicSplineSet.h"
#include <string>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>



namespace ohmmsqmc
{


TriCubicSplineSet::TriCubicSplineSet(): DeviceGrid(NULL),
  Schr_Grid(NULL),
  m_set(NULL) { }

TriCubicSplineSet::~TriCubicSplineSet()
{
}

bool TriCubicSplineSet::put(xmlNodePtr cur,
                            Grid3D* agrid)
{
  ///  make these member functions (WHY???)
  std::vector<std::string> wfile;
  std::vector<int> norbs_t;
  DeviceGrid = agrid;  /// point to the whole Grid3D
  if(!m_set)
    m_set = new SetSplinePoint;   /// create r-point set
  xmlNodePtr node = cur->xmlChildrenNode;
  while( node != NULL)
  {
    std::string name((char*)node->name);
    if(name=="orbitals")
    {
      int n = atoi((char*)xmlGetProp(node,(xmlChar*)"norbs"));
      std::string w((char*)xmlGetProp(node,(xmlChar*)"file"));
      norbs_t.push_back(n);
      wfile.push_back(w);
    }
    if(name=="HeteroStructure")
    {
      xmlNodePtr node1 = node->xmlChildrenNode;
      while(node1 != NULL)
      {
        std::string cname((const char*)(node1->name));
        if(cname=="Grid3D")
        {
          xmlNodePtr node2 = node1->xmlChildrenNode;
          while(node2 != NULL)
          {
            std::string qname((const char*)(node2->name));
            if(qname=="qGridi")
            {
              ri[0] = atoi((char*)xmlGetProp(node2,(xmlChar*)"x"));
              ri[1] = atoi((char*)xmlGetProp(node2,(xmlChar*)"y"));
              ri[2] = atoi((char*)xmlGetProp(node2,(xmlChar*)"z"));
            }
            if(qname=="qGridf")
            {
              rf[0] = atoi((char*)xmlGetProp(node2,(xmlChar*)"x"));
              rf[1] = atoi((char*)xmlGetProp(node2,(xmlChar*)"y"));
              rf[2] = atoi((char*)xmlGetProp(node2,(xmlChar*)"z"));
            }
            node2 = node2->next;
          }
        }
        node1 =  node1->next;
      }
    }
    node = node->next;
  }
  /// allocate and initialise the Schr_Grid
  if(!Schr_Grid)
    Schr_Grid = new Grid3D;
  Schr_Grid->init(ri,rf,DeviceGrid);
  /// read in the input wave-function file
  norbs = 0;
  std::cout << "Reading Wafefunction file ..." ;
  for(int iset=0; iset<wfile.size(); iset++)
  {
    for(int iorb = 0; iorb < norbs_t[iset]; iorb++)
    {
      m_psi.push_back(new TriCubicSpline(m_set,Schr_Grid));
    }
    readwf(wfile[iset].c_str(),iset,norbs_t[iset]);
  }
  std::cout << "done!" << std::endl;
  if(wfile.size() ==1)
  {
    char aname[4];
    int num = orbital_map.size();
    for(int iorb = 0; iorb < num; iorb++)
    {
      sprintf(aname,"d%d",iorb);
      orbital_map[aname]=iorb;
    }
  }
  /// initialise the TriCubicSpline
  for(int iorb = 0; iorb < m_psi.size(); iorb++)
    m_psi[iorb]->Update(true);   /// Laplacian is required so true
  return true;
}

void TriCubicSplineSet::readwf(const char* wfile, int ispin, int num)
{
  char aname[4];
  int norbs_save = norbs;
  for(int iorb = 0; iorb < num; iorb++)
  {
    if(ispin == 0)
    {
      sprintf(aname,"u%d",iorb);
      orbital_map[aname]=norbs++;
    }
    else
    {
      sprintf(aname,"d%d",iorb);
      orbital_map[aname]=norbs++;
    }
  }
  std::ifstream infile(wfile, ios_base::in);
  double x;
  for(int i = 0; i < Schr_Grid->m_size; i++)
  {
    for(int iorb = norbs_save, j=0; j< num; iorb++,j++)
    {
      infile >>x;
      const gridvec_t& ir = Schr_Grid->ipt(i);
      (*m_psi[iorb])(ir[0],ir[1],ir[2]) = -x;
    }
  }
  return;
}

void TriCubicSplineSet::evaluate(const posvec_t& r,
                                 scalar_array_t& fval,
                                 posarray_t& gradf,scalar_array_t& lapf)
{
  for(int iorb = 0; iorb < norbs; iorb++)
    fval[iorb] = m_psi[iorb]->evaluate(r,gradf[iorb],lapf[iorb]);
  return;
}

}
