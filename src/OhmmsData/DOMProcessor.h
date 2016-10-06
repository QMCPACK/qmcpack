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
    std::cout << "returning a xpathobject of " << apath << std::endl;
    return t_result = xmlXPathEvalExpression((const xmlChar*)apath,m_context);
  }
};

}
#endif
#endif
