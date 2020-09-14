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
 *  This combined with a block like which must call barrier_and_abort for the 
 *  MPI comm the error is thrown in.
 * 
 *  catch (const UniformMPIerror& ue)
 *  {
 *    my_local_comm->barrier_and_abort(ue.what());
 *  }
 *
 *  insures that the error message actually makes it out on the head
 *  rank before the MPI processes are aborted.
 *  
 *  If even one rank does not experience the issue you may hang the
 *  job until the walltime expires.
 */
class UniformCommunicateError : public std::runtime_error
{
public:
  UniformCommunicateError(const std::string& what_arg) : std::runtime_error(what_arg) {}
  UniformCommunicateError(const char* what_arg) : std::runtime_error(what_arg) {}
};

} // namespace qmcplusplus

#endif
