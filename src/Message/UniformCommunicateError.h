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
/** Subclass of std::runtime_error for errors occurring uniformly across all ranks of a communicator.
 *
 *  This exception is intended for collective, deterministic runtime errors that are guaranteed
 *  to occur on all ranks within a specific Communicate context (e.g., parsing invalid XML input,
 *  incompatible simulation parameters, or missing required shared resources).
 *
 *  Intended usage pattern:
 *  \code
 *      catch (const UniformCommunicateError& ue)
 *      {
 *        my_local_comm->barrier_and_abort(formatted_msg);
 *      }
 *  \endcode
 *
 *  Calling `barrier_and_abort()` synchronizes the ranks of `my_local_comm` at a barrier before
 *  aborting, ensuring that the head rank (rank 0 of `my_local_comm`) can flush the error message
 *  cleanly before the MPI processes are terminated.
 *
 *  Key considerations:
 *  - Communicator Scoping:
 *    The exception MUST be caught at the appropriate communicator level (`my_local_comm`).
 *    For example, in ensemble runs, errors occurring within an individual calculation must be
 *    caught and aborted on that job's sub-communicator (e.g., `qmcComm`), and not allowed to escape
 *    to the global communicator (`OHMMS::Controller`). Invoking `barrier_and_abort()` on a broader
 *    communicator while other groups are still running will cause an MPI deadlock.
 *  - Non-Uniform Errors:
 *    If an issue occurs on only a subset of ranks (e.g., node-local memory allocation failure or
 *    unexpected runtime failure), do NOT use `UniformCommunicateError` or `barrier_and_abort()`,
 *    as waiting at a barrier will hang the job. Instead, throw a standard exception such as
 *    `std::runtime_error`. Avoid calling `APP_ABORT(...)` directly in library or algorithm code,
 *    as this terminates the process immediately and breaks unit testability. Unhandled exceptions
 *    will unwind the stack to `main()`, where the outer exception handler safely aborts the global
 *    communicator with the failing rank ID prepended.
 *  - Contextual Error Reporting:
 *    When constructing the message for `barrier_and_abort()`, prepend relevant context such as
 *    the parallel scope, ensemble job index, and input file name to `ue.what()`. Prepending avoids
 *    needing to manipulate punctuation or trailing newlines in the original error string.
 *  - Catch Placement & Testability:
 *    Catch this exception as close to the bottom of the call stack (outermost driver layer)
 *    as possible, right before leaving the scope of the consistent communicator.
 *    Intermediate functions, builders, and factories should let the exception propagate so
 *    that unit tests can verify failure modes without terminating the process via
 *    `barrier_and_abort()`. Furthermore, the outermost caller possesses the high-level
 *    context (such as input file names and job indices) needed for informative diagnostics.
 */
class UniformCommunicateError : public std::runtime_error
{
public:
  UniformCommunicateError(const std::string& what_arg) : std::runtime_error(what_arg) {}
  UniformCommunicateError(const char* what_arg) : std::runtime_error(what_arg) {}
};

} // namespace qmcplusplus

#endif
