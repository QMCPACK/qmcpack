//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_RANDOMGENERATOR_IO_H
#define QMCPLUSPLUS_RANDOMGENERATOR_IO_H
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/HDFStringAttrib.h"
namespace qmcplusplus
{
/** Specialization for Random Generator */
template<>
struct HDFAttribIO<RandomGenerator_t>: public HDFAttribIOBase
{

  RandomGenerator_t& ref;

  HDFAttribIO<RandomGenerator_t>(RandomGenerator_t& a): ref(a) { }

  /** write the state to a hdf5 group
   * @param gid group id
   * @param overwrite if true, open the dataset
   */
  inline void write(hid_t grp, const char* name)
  {
    std::ostringstream a;
    ref.write(a);
    qmcplusplus::HDFAttribIO<std::ostringstream> o(a);
    o.write(grp,ref.EngineName.c_str());
  }

  /** read the state from a hdf5 group
   * @param gid group id
   */
  inline void read(hid_t grp, const char* name)
  {
    std::string s;
    HDFAttribIO<string> r(s);
    r.read(grp,ref.EngineName.c_str());
    std::istringstream is(s);
    ref.read(is);
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: BoostRandom.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/
