//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
// -*- C++ -*-
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/Communicate.h"
#include "Utilities/Timer.h"
#include <fstream>

OhmmsXPathObject::OhmmsXPathObject(): NumObjects(0), result(NULL)
{
}

OhmmsXPathObject::OhmmsXPathObject(const char* expression,
                                   xmlXPathContextPtr context) : NumObjects(0), result(NULL)
{
  put(expression, context);
}

OhmmsXPathObject::~OhmmsXPathObject()
{
  if(result != NULL)
  {
    xmlXPathFreeObject(result);
  }
}

void OhmmsXPathObject::put(const char* expression, xmlXPathContextPtr context)
{
  //free the object if it has been used before
  if(result != NULL)
    xmlXPathFreeObject(result);
  result = xmlXPathEvalExpression((const xmlChar*)expression,context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    NumObjects=0;
  }
  else
  {
    NumObjects=result->nodesetval->nodeNr;
  }
}

Libxml2Document::Libxml2Document(): m_doc(NULL), m_root(NULL),
  m_context(NULL) { }

Libxml2Document::Libxml2Document(const std::string& xmlfile)
{
  parse(xmlfile);
}

Libxml2Document::~Libxml2Document()
{
  if(m_context != NULL)
    xmlXPathFreeContext(m_context);
  if(m_doc != NULL)
    xmlFreeDoc(m_doc);
}


xmlXPathContextPtr Libxml2Document::getXPathContext()
{
  if(m_context == NULL && m_doc != NULL)
  {
    m_context = xmlXPathNewContext(m_doc);
  }
  return m_context;
}

void Libxml2Document::dump(const std::string& newxml)
{
  xmlSaveFormatFile(newxml.c_str(),m_doc,1);
}

void Libxml2Document::addChild(xmlNodePtr newnode)
{
  xmlKeepBlanksDefault(1);
  xmlAddChild(m_root,newnode);
}

void Libxml2Document::addChild(const std::string& expression,
                               xmlNodePtr newnode)
{
  xmlXPathContextPtr acontext=getXPathContext();
  OhmmsXPathObject res(expression.c_str(),acontext);
  if(res.size())
  {
    xmlAddChild(res[0],newnode);
  }
}

#if defined(HAVE_MPI2)
/** this parse function uses MPI_Bcast but is not used for the moment */
bool Libxml2Document::parse(const std::string& xmlfile)
{
  if(m_doc != NULL)
    xmlFreeDoc(m_doc);
  qmcplusplus::Timer aClock;
  int length=0;
  char *buffer=0;
  aClock.restart();
  if(OHMMS::Controller->master())
  {
    std::ifstream is(xmlfile.c_str());
    // get length of file:
    is.seekg (0, std::ios::end);
    length = is.tellg();
    is.seekg (0, std::ios::beg);
    // allocate memory:
    buffer = new char [length+1];
    // read data as a block:
    is.read (buffer,length);
    is.close();
  }
  MPI_Bcast(&length,1, MPI_INT, 0, OHMMS::Controller->getID());
  OHMMS::Controller->barrier();
  if(!(OHMMS::Controller->master()))
  {
    buffer = new char [length+1];
  }
  MPI_Bcast(buffer,length, MPI_CHAR, 0, OHMMS::Controller->getID());
  //////BCast did not work so great
  //if(OHMMS::Controller->master()) {
  //  OOMPI_Request_array request;
  //  for(int ip=1, node=0; ip<OHMMS::Controller->ncontexts(); ip++,node++) {
  //    request[node]=OOMPI_COMM_WORLD[ip].Isend(buffer,length);
  //  }
  //  OOMPI_Status_array status = request.Waitall();
  //} else {
  //  buffer = new char[length+1];
  //  OOMPI_Message amsg(buffer,length);
  //  OOMPI_Request request = OOMPI_COMM_WORLD[0].Irecv(amsg);
  //  OOMPI_Status status = request.Wait();
  //}
  m_doc = xmlParseMemory(buffer,length);
  delete [] buffer;
  qmcplusplus::app_log() << " Parsing " << xmlfile << " : " << aClock.elapsed() << " seconds " << endl;
  if (m_doc == NULL)
  {
    return false;
  }
  m_root = xmlDocGetRootElement(m_doc);
  if (m_root == NULL)
  {
    return false;
  }
  InFileRoot = std::string(xmlfile,0,xmlfile.size()-4);
  return true;
}
#else
bool Libxml2Document::parse(const std::string& xmlfile)
{
  if(m_doc != NULL)
    xmlFreeDoc(m_doc);
  m_doc = xmlParseFile(xmlfile.c_str());
  if (m_doc == NULL)
  {
    return false;
  }
  m_root = xmlDocGetRootElement(m_doc);
  if (m_root == NULL)
  {
    return false;
  }
  InFileRoot = std::string(xmlfile,0,xmlfile.size()-4);
  return true;
}
#endif

//#if defined(HAVE_MPI)
//    xmlChar *xmlbuff;
//    char *charbuff;
//    int buffersize;
//
//    // build an XML tree from a the file;
//    if(OHMMS::Controler()->master()) {
//      m_doc = xmlParseFile(infile.c_str());
//      if (m_doc == NULL) {
//        ERRORMSG("File " << infile << " is invalid")
//        return false;
//      }
//      xmlDocDumpFormatMemory(m_doc,&xmlbuff,&buffersize,1);
//    }
//
//    OHMMS::Controller()->bcast(buffersize);
//
//    if(OHMMS::Controler()->master()) {
//      charbuff = (char*)xmlbuff;
//    } else {
//      charbuff = new char[buffersize];
//    }
//
//    OHMMS::Controller()->bcast(charbuff,buffersize);
//
//    if(OHMMS::Controller()->master()) {
//      xmlFreeDoc(xmlbuff);
//    } else {
//      m_doc = xmlReadMemory(charbuff, buffersize,
//          "noname.xml", NULL, 0);
//      delete [] charbuff;
//    }
//#else
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
