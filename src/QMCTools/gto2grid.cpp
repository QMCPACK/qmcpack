//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.h"
using namespace qmcplusplus;

void buildBasisSet(xmlNodePtr cur);

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  // build an XML tree from a the file;
  xmlDocPtr m_doc = xmlParseFile(argv[1]);
  if (m_doc == NULL)
  {
    ERRORMSG("File " << argv[1] << " is invalid")
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
      buildBasisSet(result->nodesetval->nodeTab[ic]);
    }
  }
  xmlXPathFreeObject(result);
  return 0;
}

void buildBasisSet(xmlNodePtr cur)
{
  xmlNodePtr anchor = cur;
  xmlNodePtr grid_ptr = 0;
  vector<xmlNodePtr> phi_ptr;
  vector<QuantumNumberType> nlms;
  int Lmax = 0;
  int current = 0;
  std::string acenter("none");
  const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"species");
  if(aptr)
    acenter = (const char*)aptr;
  cout << "Building basis set for " << acenter << endl;
  cur = anchor->children;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "grid")
      grid_ptr = cur;
    else
      if(cname == "basisSet")
      {
        phi_ptr.push_back(cur);
        nlms.push_back(QuantumNumberType());
        int n=1,l=0,m=0;
        const xmlChar* lptr = xmlGetProp(cur,(const xmlChar*)"l");
        if(lptr)
          l = atoi((const char*)lptr);
        const xmlChar* nptr = xmlGetProp(cur,(const xmlChar*)"n");
        if(nptr)
          n = atoi((const char*)nptr);
        const xmlChar* mptr = xmlGetProp(cur,(const xmlChar*)"m");
        if(mptr)
          m = atoi((const char*)mptr);
        Lmax = std::max(l,Lmax);
        nlms[current][0]=n;
        nlms[current][1]=l;
        nlms[current][2]=m;
        current++;
      }
    cur = cur->next;
  }
  RGFBuilderBase::CenteredOrbitalType aos(Lmax);
  RGFBuilderBase* rbuilder = new GTO2GridBuilder;
  rbuilder->setOrbitalSet(&aos,acenter);
  rbuilder->addGrid(grid_ptr);
  for(int i=0; i<nlms.size(); i++)
  {
    rbuilder->addRadialOrbital(phi_ptr[i],nlms[i]);
  }
  rbuilder->print(acenter,1);
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
