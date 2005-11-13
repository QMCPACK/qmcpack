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
#ifndef LIBXML2_DOCUMENT_H
#define LIBXML2_DOCUMENT_H
#include "OhmmsData/libxmldefs.h"
#include <string>
/** class to handle xmlXPathObject 
 */
struct OhmmsXPathObject {

  //default constructor
  OhmmsXPathObject();

  /** constructor
   * @param expression  xpath expression
   * @param context xmlXPathContext with which expression is evaluated
   */
  OhmmsXPathObject(const char* expression, xmlXPathContextPtr context);

  ~OhmmsXPathObject();

  /** evaluate the expression and create the object
   * @param expression  xpath expression
   * @param context xmlXPathContext with which expression is evaluated
   */
  void put(const char* expression, xmlXPathContextPtr context);

  inline bool empty() { return NumObjects == 0;}

  inline int size() {return NumObjects;}

  inline xmlNodePtr operator[](int i) { 
    if(result != NULL && i < NumObjects) {
      return result->nodesetval->nodeTab[i];
    } else {
      return NULL;
    }
  }

  int NumObjects;
  xmlXPathObjectPtr result;
};

/** class that handles xmlDoc
 */
struct Libxml2Document {

  Libxml2Document();
  Libxml2Document(const std::string& fname);
  ~Libxml2Document();

  bool parse(const std::string& fname);

  inline xmlDocPtr getDocument() { return m_doc;}
  inline xmlNodePtr getRoot() { return m_root;}
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
