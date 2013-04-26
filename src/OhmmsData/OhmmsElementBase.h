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
// -*- C++ -*-
#ifndef OHMMS_XMLDATA_H
#define OHMMS_XMLDATA_H
/**@file OhmmsElementBase.h
 *@brief Declaration of OhmmsElementBase and define xml-related macros.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBXML2
#include "OhmmsData/libxmldefs.h"
#else /*HAVE_LIBXML2 */
#ifndef xmlDocPtr
#define xmlDocPtr void*
#define xmlNodePtr void*
#define xmlNsPtr void*
#endif
#endif  /*HAVE_LIBXML2 */

/**\class OhmmsElementBase
 *\brief Abstract class to provide xml-compatible I/O interfaces for the derived classes.
 *
 *Generic interfaces using std::iostream are much preferred. However,
 *there isn't any pure c++ xml parser that is based on std::iostream alone.
 *After evaluating several xml parsers, JK chose libxml
 *(The XML C parser and toolkit of gnome, http://www.xmlsoft.org)
 *based on its performance and availability on many platforms.
 *
 *The base class is written to be able to handle DTD or Schema in
 *future.  Current implementation assumes that each OhmmsElementBase
 *object handles a node and its child nodes. However, it does not
 *specify how the derived classes hanlde the child nodes.
 */
class OhmmsElementBase
{

public:

  ///enumeration to choose the xml parser
  enum {useLIBXML=0, /*!< using libxml2 library */
        useLIBXMLPP, /*!< using libxml++ library */
        usePLAIN     /*!< using ascii parser */
       };

  ///constructor with a name
  OhmmsElementBase(const char* aname = "none"):
    myIOMode(useLIBXML),  myName(aname)
  {
  }

  ///destructor
  virtual ~OhmmsElementBase() { }

  ///return the name
  inline const std::string& getName() const
  {
    return myName;
  }

  ///set name
  inline void setName(const std::string& aname)
  {
    myName = aname;
  }

  ///set iomode
  inline void setIOMode(int imode)
  {
    myIOMode = imode;
  }

  ///write to a ostream
  virtual bool get(std::ostream& ) const = 0;

  ///read from istream
  virtual bool put(std::istream& ) = 0;

  ///read from an xmlNode
  virtual bool put(xmlNodePtr cur) = 0;

  ///reset member data
  virtual void reset() = 0;

  ///add a xmlNode to the children list of parent
  virtual bool add(xmlNodePtr parent)
  {
    return true;
  }

  ///read from string
  void put(const std::string& s)
  {
    std::istringstream stream(s);
    put(stream);
  }

  ///write the start of a node
  virtual void begin_node(std::ostream& os) const { }

  ///write the end of a node
  virtual void end_node(std::ostream& os) const { }

protected:

  ///the type of IO mode: default is useLIBXML
  int myIOMode;

  ///the name of the node, corresponds to the xml tag
  std::string myName;
};

//add tolower function here

inline void tolower(std::string& s)
{
  for(int i=0; i<s.size(); ++i)
    s[i]=tolower(s[i]);
  //std::transform(s.begin(), s.end(), s.begin(), std::tolower);
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
