//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_UNIFORMCOMMUNICATEERROR_H
#define QMCPLUSPLUS_UNIFORMCOMMUNICATEERROR_H
#include <stdexcept>

namespace qmcplusplus
{
/** This a subclass for runtime errors that will occur on all ranks.
 *
 *  This is intended to be thrown and caught within a single Communicate context.
 *  
 *  Intended use is with a catch block like this:
 *
 *      catch (const UniformMPIerror& ue)
 *      {
 *        my_local_comm->barrier_and_abort(ue.what());
 *      }
 *
 *  Which insures that the error message actually makes it out on the head
 *  rank before the MPI processes are aborted.
 *  
 *  If even one member of the comm i.e. my_local_comm does not experience the issue
 *  you will HANG the job until the walltime expires. When in doubt it is better
 *  to lose the error message until the crash is examined in a debugger.
 *  
 *  Because of this you should catch the exception as soon as you have access to a Communicate
 *  instance. However unlike the unit test hostile `comm->barrier_and_abort` you can and
 *  should check that an exception is emitted so don't catch it in the same method/function.
 */
class UniformCommunicateError : public std::runtime_error
{
public:
  UniformCommunicateError(const std::string& what_arg) : std::runtime_error(what_arg) {}
  UniformCommunicateError(const char* what_arg) : std::runtime_error(what_arg) {}
};

} // namespace qmcplusplus

#endif
