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

#ifdef HAVE_ADIOS
#include "adios_read.h"
#include "adios_error.h"
#endif

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
  BranchIO(SimpleFixedNodeBranch& source, Communicate* c): ref(source),myComm(c) {}
#ifdef HAVE_ADIOS
  int64_t get_Checkpoint_size();
  void adios_checkpoint(int64_t adios_handle);
#ifdef ADIOS_VERIFY
  void adios_checkpoint_verify(ADIOS_FILE* fp);
#endif
#endif

  bool write(const string& fname);
  bool read(const string& fname);
  bool read_adios(const string& fname);
  void bcast_state();

  static vector<string> vParamName;
  static vector<string> iParamName;

  static void initAttributes();
};
}
#endif

/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.cpp,v $   $Author: jnkim $
 * $Revision: 2756 $   $Date: 2008-06-23 14:09:25 -0500 (Mon, 23 Jun 2008) $
 * $Id: SimpleFixedNodeBranch.cpp 2756 2008-06-23 19:09:25Z jnkim $
 ***************************************************************************/

