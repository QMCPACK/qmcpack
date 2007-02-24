//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
#ifndef OHMMS_LIBXML_NUMERICATTRIBIO_H
#define OHMMS_LIBXML_NUMERICATTRIBIO_H
#include "OhmmsData/libxmldefs.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
/** assign vector<T> from a node. Create a temporary vector and make assignment.
 *\param a reference vector<T> 
 *\param cur current node to which a content is copied
 *\return ture if successful
 */
template<class T>
inline bool 
putContent(APPNAMESPACE::Vector<T>& a, xmlNodePtr cur){
  std::istringstream 
    stream((const char*)
	   (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  int i=0;
  int n(a.size());
  while(!stream.eof() && i<n ){ stream >> a(i++);}
  return true;
}

template<typename T>
inline bool
putContent(APPNAMESPACE::Matrix<T>& a, xmlNodePtr cur){
  std::istringstream
    stream((const char*)
        (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  int i=0, ntot=a.size();
  while(!stream.eof() && i<ntot){ stream >> a(i++);}
  return true;
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
