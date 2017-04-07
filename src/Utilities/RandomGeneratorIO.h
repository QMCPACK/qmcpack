//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
    HDFAttribIO<std::string> r(s);
    r.read(grp,ref.EngineName.c_str());
    std::istringstream is(s);
    ref.read(is);
  }
};
}
#endif
