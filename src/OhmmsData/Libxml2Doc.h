//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana Champain
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef LIBXML2_DOCUMENT_H
#define LIBXML2_DOCUMENT_H
#include "OhmmsData/libxmldefs.h"
#include <string>
/** class to handle xmlXPathObject
 */ 
struct OhmmsXPathObject
{

  //default constructor
  OhmmsXPathObject();

  /** constructor
   * @param expression  xpath expression
   * @param context xmlXPathContext with which expression is evaluated
   */
  OhmmsXPathObject(const char* expression, xmlXPathContextPtr context);

  /** constructor
   * @param expression  xpath expression
   * @param cur xmlNodePtr
   *
   * Create m_context 
   */
  OhmmsXPathObject(const char* expression, xmlNodePtr cur);

  ~OhmmsXPathObject();

  /** evaluate the expression and create the object
   * @param expression  xpath expression
   * @param context xmlXPathContext with which expression is evaluated
   */
  void put(const char* expression, xmlXPathContextPtr context);

  inline bool empty()
  {
    return NumObjects == 0;
  }

  inline int size()
  {
    return NumObjects;
  }

  inline xmlNodePtr operator[](int i)
  {
    if(result != NULL && i < NumObjects)
    {
      return result->nodesetval->nodeTab[i];
    }
    else
    {
      return NULL;
    }
  }

  int NumObjects;
  xmlXPathObjectPtr result;
  xmlXPathContextPtr m_context;
};

/** class that handles xmlDoc
 */
struct Libxml2Document
{

  Libxml2Document();
  Libxml2Document(const std::string& fname);
  ~Libxml2Document();

  bool parse(const std::string& fname);
  bool parseFromString(const std::string& data);

  inline xmlDocPtr getDocument()
  {
    return m_doc;
  }
  inline xmlNodePtr getRoot()
  {
    return m_root;
  }
  xmlXPathContextPtr getXPathContext();

  void dump(const std::string& newxml);
  void addChild(xmlNodePtr newnode);

  void addChild(const std::string& expression, xmlNodePtr newnode);

  xmlDocPtr m_doc;
  xmlNodePtr m_root;
  xmlXPathContextPtr m_context;
  std::string InFileRoot;
};

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
