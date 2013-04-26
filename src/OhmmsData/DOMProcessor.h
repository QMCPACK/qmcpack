//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DOMPROCESSOR_H
#define QMCPLUSPLUS_DOMPROCESSOR_H
#if defined(ENABLE_LIBXML2)
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

namespace qmcplusplus
{

/**\class DOMProcessor
 *@brief a temporary class to save the Document Pointer
 */
struct DOMProcessor
{

  xmlDocPtr m_doc;
  xmlNsPtr m_ns;
  xmlNodePtr m_cur;
  xmlXPathContextPtr m_context;
  xmlXPathObjectPtr t_result;

  DOMProcessor(const char* fname)
  {
    // build an XML tree from a the file;
    m_doc = xmlParseFile(fname);
    if (m_doc == NULL)
      return;
    // Check the document is of the right kind
    m_cur = xmlDocGetRootElement(m_doc);
    if (m_cur == NULL)
    {
      fprintf(stderr,"empty document\n");
      xmlFreeDoc(m_doc);
      return;
    }
    m_ns = NULL;
    m_context = xmlXPathNewContext(m_doc);
    t_result = NULL;
  }

  ~DOMProcessor()
  {
    xmlFreeDoc(m_doc);
    if(t_result)
      xmlXPathFreeObject(t_result);
  }

  xmlXPathObjectPtr
  get(const char* apath)
  {
    //if(t_result) xmlXPathFreeObject(t_result);
    cout << "returning a xpathobject of " << apath << endl;
    return t_result = xmlXPathEvalExpression((const xmlChar*)apath,m_context);
  }
};

}
#endif
#endif
