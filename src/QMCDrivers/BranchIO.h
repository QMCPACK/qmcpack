//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#include <boost/archive/text_oarchive.hpp>

#ifndef QMCPLUSPLUS_BRANCHIO_H
#define QMCPLUSPLUS_BRANCHIO_H
namespace qmcplusplus 
{
  struct BranchIO
  {
    typedef SimpleFixedNodeBranch::RealType RealType;
    typedef SimpleFixedNodeBranch::BranchModeType BranchModeType;
    typedef SimpleFixedNodeBranch::IParamType IParamType;
    typedef SimpleFixedNodeBranch::VParamType VParamType;

    SimpleFixedNodeBranch& ref;
    Communicate* myComm;
    BranchIO(SimpleFixedNodeBranch& source, Communicate* c): ref(source),myComm(c){}
    bool write(const string& fname);
    bool read(const string& fname);
  };
}
#endif

/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.cpp,v $   $Author: jnkim $
 * $Revision: 2756 $   $Date: 2008-06-23 14:09:25 -0500 (Mon, 23 Jun 2008) $
 * $Id: SimpleFixedNodeBranch.cpp 2756 2008-06-23 19:09:25Z jnkim $ 
 ***************************************************************************/

