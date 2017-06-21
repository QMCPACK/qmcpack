//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_LIBXML_NUMERICATTRIBIO_H
#define OHMMS_LIBXML_NUMERICATTRIBIO_H
#include "OhmmsData/libxmldefs.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
/** assign std::vector<T> from a node. Create a temporary vector and make assignment.
 *\param a reference std::vector<T>
 *\param cur current node to which a content is copied
 *\return ture if successful
 */
template<class T>
inline bool
putContent(qmcplusplus::Vector<T>& a, xmlNodePtr cur)
{
  std::istringstream
  stream((const char*)
         (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  int i=0;
  int n(a.size());
  while(!stream.eof() && i<n )
  {
    stream >> a(i++);
  }
  return true;
}

template<typename T>
inline bool
putContent(qmcplusplus::Matrix<T>& a, xmlNodePtr cur)
{
  std::istringstream
  stream((const char*)
         (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  int i=0, ntot=a.size();
  while(!stream.eof() && i<ntot)
  {
    stream >> a(i++);
  }
  return true;
}

#endif
