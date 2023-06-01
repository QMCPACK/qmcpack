//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_FILEUTILITY_H
#define OHMMS_FILEUTILITY_H

#include <string_view>

inline std::string_view getExtension(const std::string_view str)
{
  size_t pos = str.find_last_of('.');
  if (pos == std::string_view::npos)
    return std::string_view();
  return str.substr(pos + 1);
}
#endif
