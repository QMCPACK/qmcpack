//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Alfredo A. Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_COMMUNICATE_H
#define OHMMS_COMMUNICATE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MPI
#include "mpi3/environment.hpp"
namespace mpi3 = boost::mpi3;
#endif

#ifdef HAVE_MPI
struct CommunicatorTraits
{
  using mpi_comm_type = MPI_Comm;
};

#else
struct CommunicatorTraits
{
  using mpi_comm_type            = int;
  static const int MPI_COMM_NULL = 0;
};
#endif

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <unistd.h>

#include "Message/AppAbort.h"

/**@class Communicate
 * @ingroup Message
 * @brief
 *  Wrapping information on parallelism.
 *  Very limited in functions. Currently, only single-mode or mpi-mode
 *  is available (mutually exclusive).
 * @todo Possibly, make it a general manager class for mpi+openmp, mpi+mpi
 */
class Communicate : public CommunicatorTraits
{
public:
  ///constructor
  Communicate();

#ifdef HAVE_MPI
  /** Constructor that takes ownership of an existing or temporary communicator.
   *  @param in_comm Communicator to move from (e.g., returned by env.world(), split(), or split_shared()).
   *  Transfers ownership via move semantics directly without the collective overhead of `MPI_Comm_dup`.
   *  `in_comm` is left in an empty/null state.
   */
  Communicate(mpi3::communicator&& in_comm);
#endif

  /**destructor
   * Call proper finalization of Communication library
   */
  ~Communicate();

  ///disable copy constructor
  Communicate(const Communicate&) = delete;

  Communicate(Communicate&&);

  /** Create a new communicator by partitioning the current communicator into `nparts` subgroups.
   *
   *  The ranks are partitioned fairly using FairDivideLow. When `multiple` > 1 is specified,
   *  ranks are grouped in contiguous chunks of size `multiple` so that the size of each
   *  resulting sub-communicator is guaranteed to be a multiple of `multiple`.
   *
   *  @param nparts   The number of sub-communicator groups to divide into.
   *  @param multiple The granularity of each sub-communicator size. `size()` must be divisible
   *                  by `multiple`, and `multiple` must be >= 1. Defaults to 1.
   *  @return The sub-communicator corresponding to the calling rank's subgroup.
   *
   *  Example 1: size() == 12, nparts = 2, multiple = 1
   *    | group 0      | group 1       |
   *    |  0 1 2 3 4 5 | 6 7 8 9 10 11 |
   *
   *  Example 2: size() == 10, nparts = 3, multiple = 2
   *    | group 0 | group 1 | group 2 |
   *    |  0 1    |  2 3 4 5 | 6 7 8 9 |
   */
  Communicate split(int nparts, int multiple = 1) const;

  /** create a new communicator by splitting the current communicator with stride
   * size() == 12 ; nparts = 2; stride = 3;
   *   | group 0 | group 1 |
   *   |  0 1 2  |  3 4 5  |
   *   |  6 7 8  | 9 10 11 |
   */
  Communicate splitByStride(int nparts, int stride) const;

  /** create a node/shared-memory communicator from current (parent) communicator
   *
   *  Note: the size of the node communicator may differ on different nodes.
   */
  Communicate createNodeComm() const;


  /** return the rank count of a uniform node
   *
   *  A "uniform node" is defined as a logical grouping of ranks that resides entirely within
   *  a single physical node, where all uniform nodes across the communicator have an identical
   *  number of ranks.
   *
   *  Because different physical nodes may host varying numbers of ranks, the rank count of a
   *  uniform node is determined by the greatest common divisor (GCD) of all local node communicator
   *  sizes across the communicator. This ensures every physical node can be partitioned into an
   *  integer number of identically sized sub-node groups.
   *
   *  In non-MPI (serial) builds, this returns 1.
   */
  int getUniformNodeRankCount() const;

  void barrier() const;
  void abort() const;
  void barrier_and_abort(const std::string& msg) const;
  void set_world();

#if defined(HAVE_MPI)
  ///operator for implicit conversion to MPI_Comm
  operator MPI_Comm() const { return myMPI; }
  ///return reference to mpi3::communicator
  mpi3::communicator& getCommMPI3() { return comm; }
  const mpi3::communicator& getCommMPI3() const { return comm; }
#endif

  ///return the Communicator ID (typically MPI_WORLD_COMM)
  mpi_comm_type getMPI() const { return myMPI; }

  ///return the rank
  int rank() const { return d_mycontext; }
  ///return the number of tasks
  int size() const { return d_ncontexts; }

  ///return the group id
  int getGroupID() const { return d_groupid; }
  ///return the number of intra_comms which belong to the same group
  int getNumGroups() const { return d_ngroups; }
  // Avoid public access to unique_ptr.
  Communicate& getInterGroupComm() const;

  ///return true if the current MPI rank is the group lead
  bool isGroupLeader() { return d_mycontext == 0; }

  // MMORALES: leaving this here temprarily, but it doesn;t belong here.
  // MMORALES: FIX FIX FIX
#ifdef HAVE_MPI

// For Mac OS X
#ifndef HOST_NAME_MAX
#ifdef _POSIX_HOST_NAME_MAX
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#endif
#endif

#endif

#ifdef HAVE_MPI
  /** A hack to get around Communicate not supporting flexible processor subgroups
   *
   *  MMORALES:
   *  right now there is no easy way to use Communicate
   *  for generic processor subgroups, so calling split on myMPI
   *  and managing the communicator directly
   *  \todo THIS MUST BE FIXED!!!
   */
  void split_comm(int key, MPI_Comm& comm)
  {
    int myrank = rank();
    MPI_Comm_split(myMPI, key, myrank, &comm);
  }
#endif

  template<typename T>
  void allreduce(T&);
  template<typename T>
  void reduce(T&, int dest = 0);
  template<typename T>
  void reduce(const T* sb, T* rb, int n, int dest = 0);
  template<typename T>
  void reduce_in_place(T*, int n, int dest = 0);
  template<typename T>
  void bcast(T&, int root = 0);
  template<typename T>
  void bcast(T*, int n, int root = 0);
  template<typename T>
  void gather(T& sb, T& rb, int dest = 0);
  template<typename T, typename IT>
  void gatherv(T& sb, T& rb, IT& counts, IT& displ, int dest = 0);
  template<typename T>
  void scatter(T& sb, T& rb, int dest = 0);
  template<typename T, typename IT>
  void scatterv(T& sb, T& rb, IT& counts, IT& displ, int source = 0);
  template<typename T, typename IT>
  void gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest = 0);
  template<typename T, typename TMPI, typename IT>
  void gatherv_in_place(T* buf, const TMPI& datatype, IT& counts, IT& displ, int dest = 0);
  template<typename T>
  void allgather(T& sb, T& rb);
  template<typename T>
  void allgather(const T* sb, T* rb, int count);


private:
  /** Raw communicator
   *
   *  Currently it is only owned by Communicate which manages its creation and destruction
   *  After switching to mpi3::communicator, myMPI is only a reference to the raw communicator owned by mpi3::communicator
   */
  mpi_comm_type myMPI;
  /// Rank
  int d_mycontext;
  /// Size
  int d_ncontexts;
  /// Group ID of the current communicator in the parent communicator
  int d_groupid;
  /// Total number of groups in the parent communicator
  int d_ngroups;
  /// inter group communicator
  std::unique_ptr<Communicate> inter_group_comm_;

#ifdef HAVE_MPI
  /// mpi3 communicator wrapper
  mutable mpi3::communicator comm;
#endif
};


namespace OHMMS
{
/** Global Communicator for a process
 */
#include <memory>
extern std::unique_ptr<Communicate> Controller;
} // namespace OHMMS


#endif // OHMMS_COMMUNICATE_H
