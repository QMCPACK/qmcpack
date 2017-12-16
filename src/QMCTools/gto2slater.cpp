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
    
    



#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/GaussianBasisSet.h"
#include "QMCTools/Any2Slater.h"

struct GTO2Slater
{

  typedef GaussianCombo<double>      GTOType;
  typedef LogGrid<double>            GridType;

  bool Normalized;
  GridType myGrid;
  xmlNodePtr gridPtr;
  std::map<std::string,xmlNodePtr> sPtr;

  GTO2Slater():Normalized(false),gridPtr(0) {}
  int parse(const char* fname);
  bool put(xmlNodePtr);
  void optimize();
};

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  GTO2Slater transformer;
  transformer.parse(argv[1]);
  return 0;
}

int GTO2Slater::parse(const char* fname)
{
  // build an XML tree from a the file;
  xmlDocPtr m_doc = xmlParseFile(fname);
  if (m_doc == NULL)
  {
    ERRORMSG("File " << fname << " is invalid")
    xmlFreeDoc(m_doc);
    return 1;
  }
  // Check the document is of the right kind
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if (cur == NULL)
  {
    ERRORMSG("Empty document");
    xmlFreeDoc(m_doc);
    return 1;
  }
  xmlXPathContextPtr m_context = xmlXPathNewContext(m_doc);
  xmlXPathObjectPtr result
  = xmlXPathEvalExpression((const xmlChar*)"//atomicBasisSet",m_context);
  if(!xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    for(int ic=0; ic<result->nodesetval->nodeNr; ic++)
    {
      std::cout << "Going to optimize" << std::endl;
      if(put(result->nodesetval->nodeTab[ic]))
      {
        optimize();
      }
    }
  }
  xmlXPathFreeObject(result);
  return 1;
}

bool GTO2Slater::put(xmlNodePtr cur)
{
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "grid")
      gridPtr = cur;
    else
      if(cname == "basisGroup")
      {
        std::string rid("invalid");
        std::string rtype("Gaussian");
        std::string norm("no");
        int l=0;
        OhmmsAttributeSet inAttrib;
        inAttrib.add(rid,"rid");
        inAttrib.add(l,"l");
        inAttrib.add(rtype,"type");
        inAttrib.add(norm,"normalized");
        inAttrib.put(cur);
        if(rtype == "Gaussian" && l == 0)
          //pick only S
        {
          //if Ngto==1, don't do it
          if(norm == "yes")
            Normalized=true;
          else
            Normalized=false;
          std::map<std::string,xmlNodePtr>::iterator it(sPtr.find(rid));
          if(it == sPtr.end())
          {
            sPtr[rid]=cur;
          }
        }
      }
    cur=cur->next;
  }
  if(sPtr.empty())
    return false;
  return
    true;
}

/** main functio to optimize multiple contracted S orbitals
 */
void GTO2Slater::optimize()
{
  //construct one-dim grid
  double ri = 1e-5;
  double rf = 10.0;
  int npts = 101;
  std::string gridType("log");
  if(gridPtr)
  {
    OhmmsAttributeSet radAttrib;
    radAttrib.add(gridType,"type");
    radAttrib.add(npts,"npts");
    radAttrib.add(ri,"ri");
    radAttrib.add(rf,"rf");
    radAttrib.put(gridPtr);
  }
  myGrid.set(ri,rf,npts);
  //create a numerical grid funtor
  typedef OneDimCubicSpline<double> RadialOrbitalType;
  RadialOrbitalType radorb(&myGrid);
  int L= 0;
  //Loop over all the constracted S orbitals
  std::map<std::string,xmlNodePtr>::iterator it(sPtr.begin()),it_end(sPtr.end());
  while(it != it_end)
  {
    //create contracted gaussian
    GTOType gset(L,Normalized);
    //read the radfunc's of basisGroup
    gset.putBasisGroup((*it).second);
    //convert to a radial functor
    Transform2GridFunctor<GTOType,RadialOrbitalType> transform(gset, radorb);
    transform.generate(myGrid.rmin(),myGrid.rmax(),myGrid.size());
    //optimize it with the radial functor
    Any2Slater gto2slater(radorb);
    gto2slater.optimize();
    ++it;
  }
  sPtr.clear();
}
