//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H
#define QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H

#include "QMCWaveFunctions/EinsplineSetBuilder.h"

namespace qmcplusplus
{
  
class MockEinsplineSetBuilder : public EinsplineSetBuilder
{
public:
  MockEinsplineSetBuilder() {};
  MockEinsplineSetBuilder(const MockEinsplineSetBuilder&) = default;
};

}
#endif //QMCPLUSPLUS_MOCKEINSPLINESETBUILDER_H
