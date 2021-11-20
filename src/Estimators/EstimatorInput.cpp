//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include "EstimatorInput.h"

/** \file
 *  collected input checks common across estimators
 */
namespace qmcplusplus
{
namespace estimatorinput
{

void checkCenterCorner(InputSection& input_section, const std::string& error_tag)
{
  if (input_section.has("center") && input_section.has("corner"))
    throw UniformCommunicateError(error_tag + " cannot defined both center and corner.");
}

}
} // namespace qmcplusplus
