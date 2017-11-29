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
    
    



#include "Configuration.h"
#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus
{

QMCAppBase::QMCAppBase()
{
}

QMCAppBase::~QMCAppBase()
{
  while(!XmlDocStack.empty())
  {
    popDocument();
  }
}

bool QMCAppBase::pushDocument(const std::string& infile)
{
  Libxml2Document* adoc= new Libxml2Document();
  bool success = adoc->parse(infile);
  if(success)
  {
    XmlDocStack.push(adoc);
  }
  else
  {
    app_error() << "File " << infile << " is invalid" << std::endl;
    delete adoc;
  }
  return success;
}

void QMCAppBase::popDocument()
{
  if(!XmlDocStack.empty())
    //Check if the stack is empty
  {
    Libxml2Document* adoc=XmlDocStack.top();
    delete adoc;
    XmlDocStack.pop();
  }
}

/** parse an input file
 * @param infile name of an input file
 * @return true, if the document is a valid xml file.
 *
 * The data members m_doc and m_root point to the "current" document and
 * root element.
 */
bool QMCAppBase::parse(const std::string& infile)
{
  app_summary() << "  Input XML = " << infile << std::endl;
  return pushDocument(infile);
}

void QMCAppBase::saveXml()
{
  if(!XmlDocStack.empty())
  {
    std::string newxml(myProject.CurrentMainRoot());
    //string newxml(myProject.CurrentRoot());
    //myProject.PreviousRoot(newxml);
    //myProject.rewind();
    newxml.append(".cont.xml");
    app_log() << "\n========================================================="
              << "\n  A new xml input file : " << newxml << std::endl;
    XmlDocStack.top()->dump(newxml);
  }
}

std::string &QMCAppBase::getTitle()
{
  return myProject.m_title;
}
}
