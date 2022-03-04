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

#ifndef QMCPLUSPLUS_ESIMATORINPUT_H
#define QMCPLUSPLUS_ESIMATORINPUT_H

#include <string>
#include "Configuration.h"
#include "InputSection.h"

/** If tag is present in input_secution set its variable.
 *  Many variables have default values we don't want overwritten and that we want expressed in native c++
 *  Define macro to avoid repeating code for this lambda
 */
#define LAMBDA_setIfInInput [&](auto& var, const std::string& tag) -> bool { return input_section_.setIfInInput(var, tag); }

namespace qmcplusplus
{

namespace estimatorinput
{

void checkCenterCorner(InputSection& input_section, const std::string& error_tag);


} // namespace estimatorinput
} // namespace qmcplusplus
#endif
