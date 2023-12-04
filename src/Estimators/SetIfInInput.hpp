//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SETIFININPUT_H
#define QMCPLUSPLUS_SETIFININPUT_H
namespace qmcplusplus {

template<>
inline bool InputSection::setIfInInput<float>(float& var, const std::string& tag)
{
  if (has(tag))
  {
    var = get<FullPrecReal>(tag);
    return true;
  }
  else
    return false;
}


template<typename T>
inline bool InputSection::setIfInInput(T& var, const std::string& tag)
{
  if (has(tag))
  {
    var = get<T>(tag);
    return true;
  }
  else
    return false;
}

}
#endif
