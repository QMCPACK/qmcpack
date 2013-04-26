//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file VarList.h
 * @brief Define VarRegistry
 */
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
