//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MEMORYUSAGE_H
#define QMCPLUSPLUS_MEMORYUSAGE_H

#include <iostream>
#include <string>

namespace qmcplusplus
{

void print_mem(const std::string& title, std::ostream& log);

}
#endif
