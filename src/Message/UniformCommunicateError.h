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

#include <exception>
namespace qmcplusplus
{
/** This a subclass for runtime errors that will occur on all ranks.
 *
 *  This is intended to be thrown and caught within a single Communicate context.
 *  
 *  The intended use is with a catch block like this:
 *
 *  catch (const UniformMPIerror& ue)
 *  {
 *    my_local_comm->barrier_and_abort(ue.what());
 *  }
 *
 *  Which insures that the error message actually makes it out on the head
 *  rank before the MPI processes are aborted.
 *  
 *  If even one member of the comm i.e. my_local_comm does not experience the issue
 *  you will HANG the job until the walltime expires.
 *  
 *  If you split a communicator it is your responsibility to add a try catch block in your
 *  communicator scope to prevent the classes using your new communicator from propogating
 *  UniformCommunicateError into wider Comunicate scopes.
 */
class UniformCommunicateError : public std::runtime_error
{
public:
  UniformCommunicateError(const std::string& what_arg) : std::runtime_error(what_arg) {}
  UniformCommunicateError(const char* what_arg) : std::runtime_error(what_arg) {}
};

} // namespace qmcplusplus

#endif
