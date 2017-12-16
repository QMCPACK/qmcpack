//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <iostream>
#include "SQD/SQDFrame.h"
#include "OhmmsApp/ProjectData.h"

#ifdef HAVE_QT
#include <qapplication.h>

main(int argc, char **argv)
{
  OhmmsInfo Welcome("sqd");
  QApplication a( argc, argv );
  if(argc<2)
  {
    std::cerr << "Usage: sqd file-name [-nox]" << std::endl;
    std::cerr << "Use -nox if you like to skip plotting the results." << std::endl;
    return 1;
  }
  bool showplot = true;
  int i=1;
  while(i<argc)
  {
    std::string arg(argv[i]);
    if(arg == "-nox")
    {
      showplot = false;
    }
    i++;
  }
  SQDFrame *w = new SQDFrame;
  a.setMainWidget( w );
  if(w->solve(argv[1]))
  {
    w->show();
    return a.exec();
  }
  else
  {
    return 1;
  }
}
#else
main(int argc, char **argv)
{
  OhmmsInfo Welcome("sqd");
  if(argc<2)
  {
    std::cerr << "Usage: sqd file-name [-nox]" << std::endl;
    std::cerr << "Use -nox if you like to skip plotting the results." << std::endl;
    return 1;
  }
  SQDFrame sqd;
  if(sqd.solve(argv[1]))
  {
    sqd.show();
  }
}
#endif

bool
SQDFrame::solve(const char* fname)
{
  xmlNsPtr ns=NULL;
  // build an XML tree from a the file;
  m_doc = xmlParseFile(fname);
  if (m_doc == NULL)
    return false;
  //using XPath instead of recursive search
  xmlXPathContextPtr m_context = xmlXPathNewContext(m_doc);
  // Check the document is of the right kind
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if (cur == NULL)
  {
    fprintf(stderr,"empty document\n");
    xmlFreeDoc(m_doc);
    return false;
  }
  if (xmlStrcmp(cur->name, (const xmlChar *) "simulation"))
  {
    fprintf(stderr,"document of the wrong type, root node != simulation\n");
    xmlFreeDoc(m_doc);
    return false;
  }
  //project description, assign id and series
  qmcplusplus::ProjectData myProject;
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
  using namespace ohmmshf;
  HFSolver = new HartreeFock(Pot,Psi);
  bool success = HFSolver->put(cur);
  if(!success)
  {
    ERRORMSG("The input file does not conform. Exit")
    return false;
  }
  HFSolver->setRoot(myProject.CurrentRoot());
  success = HFSolver->solve();
  xmlFreeDoc(m_doc);
  xmlCleanupParser();
  return success;
}
