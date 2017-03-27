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
    
    



#include <fstream>
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCFactory/OneDimGridFactory.h"

using namespace qmcplusplus;

struct ECPTest
{

  typedef OneDimGridFactory::RealType RealType;
  typedef OneDimGridFactory::GridType GridType;
  typedef GaussianTimesRN<RealType> InFuncType;
  typedef OneDimCubicSpline<RealType> OutFuncType;

  std::string Name;
  RealType Zeff;

  ECPTest(const std::string& fname);
  void buildSemiLocal(xmlNodePtr cur);
  void buildLocal(xmlNodePtr cur);
};

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  int iargc=0;
  ECPTest ecp(argv[1]);
  return 0;
}

ECPTest::ECPTest(const std::string& fname)
{
  // build an XML tree from a the file;
  xmlDocPtr m_doc = xmlParseFile(fname.c_str());
  if (m_doc == NULL)
  {
    ERRORMSG("File " << fname << " is invalid")
    xmlFreeDoc(m_doc);
  }
  // Check the document is of the right kind
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if (cur == NULL)
  {
    ERRORMSG("Empty document");
    xmlFreeDoc(m_doc);
  }
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "header")
    {
      Zeff = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"zval"));
      Name = (const char*)xmlGetProp(cur,(const xmlChar*)"symbol");
    }
    else
      if(cname == "semilocal")
      {
        buildSemiLocal(cur);
      }
      else
        if(cname == "local")
        {
          buildLocal(cur);
        }
    cur=cur->next;
  }
  xmlFreeDoc(m_doc);
}

void ECPTest::buildSemiLocal(xmlNodePtr cur)
{
  GridType* grid_semilocal=0;
  std::vector<InFuncType*> semilocal;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "grid")
    {
      grid_semilocal=OneDimGridFactory::createGrid(cur);
    }
    else
      if(cname == "vps")
      {
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL)
        {
          if(xmlStrEqual(cur1->name,(const xmlChar*)"basisGroup"))
          {
            InFuncType* a=new InFuncType;
            a->putBasisGroup(cur1);
            semilocal.push_back(a);
          }
          cur1=cur1->next;
        }
      }
    cur=cur->next;
  }
  std::string fname(Name);
  fname.append(".semilocal.dat");
  std::ofstream fout(fname.c_str());
  fout.setf(std::ios::scientific, std::ios::floatfield);
  fout.precision(12);
  int ig=0;
  while(ig<grid_semilocal->size())
  {
    double r=(*grid_semilocal)[ig++];
    fout << std::setw(20) << setprecision(12) << r;
    for(int i=0; i<semilocal.size(); i++)
    {
      fout << std::setw(20) << semilocal[i]->f(r);
    }
    fout << std::endl;
  }
}

void ECPTest::buildLocal(xmlNodePtr cur)
{
  GridType* grid_semilocal=0;
  InFuncType* localpp;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "grid")
    {
      grid_semilocal=OneDimGridFactory::createGrid(cur);
    }
    else
      if(cname == "basisGroup")
      {
        localpp=new InFuncType;
        localpp->putBasisGroup(cur);
      }
    cur=cur->next;
  }
  std::cout << " Effective Z = " << Zeff << std::endl;
  std::string fname(Name);
  fname.append(".local.dat");
  std::ofstream fout(fname.c_str());
  fout.setf(std::ios::scientific, std::ios::floatfield);
  fout.precision(12);
  int ig=0;
  while(ig<grid_semilocal->size())
  {
    double r=(*grid_semilocal)[ig++];
    fout << std::setw(20) << setprecision(12) << r << std::setw(20) << localpp->f(r)-Zeff/r<< std::endl;
  }
  delete localpp;
}
