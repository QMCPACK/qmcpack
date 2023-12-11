////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Communicate.h
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_APPABORT_H
#define QMCPLUSPLUS_APPABORT_H

#include <string>
#include <sstream>

/** break on this function to catch any APP_ABORT call in debugger
 */
void breakableAppAbort(const std::string& str_msg);

/** Widely used but deprecated fatal error macros from legacy code
 *
 *  they allow rather odd use of the << operator for the msg argument
 *  so removing them is non trivial.
 */
#define APP_ABORT(msg)                                   \
  {                                                      \
    std::ostringstream error_message;                    \
    error_message << "Fatal Error. Aborting at " << msg; \
    breakableAppAbort(error_message.str());              \
  }

#define APP_ABORT_TRACE(f, l, msg)                                                  \
  {                                                                                 \
    std::ostringstream error_message;                                               \
    error_message << "Fatal Error. Aborting at " << l << "::" << f << "\n " << msg; \
    breakableAppAbort(error_message.str());                                         \
  }

#endif
