//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_OPTIMIZE_VARREGISTRY_H
#define QMCPLUSPLUS_OPTIMIZE_VARREGISTRY_H
#include <map>

/** Set of variables that can be modified during optimizations
 */
template<class T>
struct VarRegistry: public std::map<std::string,T>
{

  typedef typename std::map<std::string,T>::iterator iterator;

  void print(std::ostream& os)
  {
    iterator it(this->begin());
    os << "Optimizable variable list " << std::endl;
    while(it != this->end())
    {
      os << (*it).first << " " << (*it).second << std::endl;
      ++it;
    }
  }

  //add an function not
  inline int addVariable(const std::string& vname, T val)
  {
    this->operator[](vname)=val;
    return this->size();
  }

};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
