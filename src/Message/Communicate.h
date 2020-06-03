//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
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
  typedef MPI_Comm mpi_comm_type;
  typedef MPI_Status status;
  typedef MPI_Request request;
};
#define APP_ABORT(msg)                                            \
  {                                                               \
    std::cerr << "Fatal Error. Aborting at " << msg << std::endl; \
    MPI_Abort(MPI_COMM_WORLD, 1);                                 \
  }

#define APP_ABORT_TRACE(f, l, msg)                                                           \
  {                                                                                          \
    std::cerr << "Fatal Error. Aborting at " << l << "::" << f << "\n " << msg << std::endl; \
    MPI_Abort(MPI_COMM_WORLD, 1);                                                            \
  }

#else
struct CommunicatorTraits
{
  typedef int mpi_comm_type;
  typedef int status;
  typedef int request;
  static const int MPI_COMM_NULL    = 0;
  static const int MPI_REQUEST_NULL = 1;
};

#define APP_ABORT(msg)                                            \
  {                                                               \
    std::cerr << "Fatal Error. Aborting at " << msg << std::endl; \
    std::cerr.flush();                                            \
    exit(1);                                                      \
  }

#define APP_ABORT_TRACE(f, l, msg)                                                           \
  {                                                                                          \
    std::cerr << "Fatal Error. Aborting at " << l << "::" << f << "\n " << msg << std::endl; \
    exit(1);                                                                                 \
  }

#endif

#include <string>
#include <vector>
#include <utility>
#include <unistd.h>
#include <cstring>

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

  ///constructor from mpi3 environment
#ifdef HAVE_MPI
  Communicate(const mpi3::environment& env);

  ///constructor with communicator
  Communicate(const mpi3::communicator& in_comm);
#endif

  ///constructor with communicator
  Communicate(const Communicate& in_comm, int nparts);

  /**destructor
   * Call proper finalization of Communication library
   */
  virtual ~Communicate();

  ///disable constructor
  Communicate(const Communicate&) = delete;

  // Only for unit tests
  void initialize(int argc, char** argv);

#ifdef HAVE_MPI
  void initialize(const mpi3::environment& env);
#endif
  /// initialize this as a node/shared-memory communicator
  void initializeAsNodeComm(const Communicate& parent);
  void finalize();
  void barrier() const;
  void abort() const;
  void barrier_and_abort(const std::string& msg) const;
  void set_world();

#if defined(HAVE_MPI)
  ///operator for implicit conversion to MPI_Comm
  inline operator MPI_Comm() const { return myMPI; }
#endif

  ///return the Communicator ID (typically MPI_WORLD_COMM)
  inline mpi_comm_type getMPI() const { return myMPI; }

  ///return the rank
  inline int rank() const { return d_mycontext; }
  ///return the number of tasks
  inline int size() const { return d_ncontexts; }

  ///return the group id
  inline int getGroupID() const { return d_groupid; }
  ///return the number of intra_comms which belong to the same group
  inline int getNumGroups() const { return d_ngroups; }
  //inline bool master() const { return (d_mycontext == 0);}
  //intra_comm_type split(int n);
  void cleanupMessage(void*);
  inline void setNodeID(int i) { d_mycontext = i; }
  inline void setNumNodes(int n) { d_ncontexts = n; }

  inline void setName(const std::string& aname) { myName = aname; }
  inline const std::string& getName() const { return myName; }

  ///return true if the current MPI rank is the group lead
  inline bool isGroupLeader() { return d_mycontext == 0; }

  // MMORALES: leaving this here temprarily, but it doesn;t belong here.
  // MMORALES: FIX FIX FIX
#ifdef HAVE_MPI

// For Mac OS X
#ifndef HOST_NAME_MAX
#ifdef _POSIX_HOST_NAME_MAX
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#endif
#endif

  inline bool head_nodes(MPI_Comm& MPI_COMM_HEAD_OF_NODES)
  {
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    int myrank = rank(), nprocs = size();
    char* dummy = new char[nprocs * HOST_NAME_MAX];
    MPI_Allgather(hostname, HOST_NAME_MAX, MPI_CHAR, dummy, HOST_NAME_MAX, MPI_CHAR, myMPI);
    bool head_of_node = true;
    for (int i = 0; i < myrank; i++)
      if (strcmp(hostname, dummy + i * HOST_NAME_MAX) == 0)
      {
        head_of_node = false;
        break;
      }
    int key = head_of_node ? 0 : 10;
    MPI_Comm_split(myMPI, key, myrank, &MPI_COMM_HEAD_OF_NODES);
    delete[] dummy;
    return head_of_node;
  }
#endif

  // MMORALES:
  // right now there is no easy way to use Communicate
  // for generic processor subgroups, so calling split on myMPI
  // and managing the communicator directly
  // THIS MUST BE FIXED!!!
#ifdef HAVE_MPI
  inline void split_comm(int key, MPI_Comm& comm)
  {
    int myrank = rank();
    MPI_Comm_split(myMPI, key, myrank, &comm);
  }
#endif

  template<typename T>
  void allreduce(T&);
  template<typename T>
  void reduce(T&);
  template<typename T>
  void reduce(T* restrict, T* restrict, int n);
  template<typename T>
  void reduce_in_place(T* restrict, int n);
  template<typename T>
  void bcast(T&);
  template<typename T>
  void bcast(T* restrict, int n);
  template<typename T>
  void send(int dest, int tag, T&);
  template<typename T>
  void gather(T& sb, T& rb, int dest = 0);
  template<typename T, typename IT>
  void gatherv(T& sb, T& rb, IT& counts, IT& displ, int dest = 0);
  template<typename T>
  void allgather(T& sb, T& rb, int count);
  template<typename T, typename IT>
  void allgatherv(T& sb, T& rb, IT& counts, IT& displ);
  template<typename T>
  void scatter(T& sb, T& rb, int dest = 0);
  template<typename T, typename IT>
  void scatterv(T& sb, T& rb, IT& counts, IT& displ, int source = 0);
  template<typename T>
  request irecv(int source, int tag, T&);
  template<typename T>
  request isend(int dest, int tag, T&);
  template<typename T>
  request irecv(int source, int tag, T*, int n);
  template<typename T>
  request isend(int dest, int tag, T*, int n);
  template<typename T, typename IT>
  void gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest = 0);
  template<typename T, typename TMPI, typename IT>
  void gatherv_in_place(T* buf, TMPI& datatype, IT& counts, IT& displ, int dest = 0);
  template<typename T>
  void allgather(T* sb, T* rb, int count);
  template<typename T>
  void gsum(T&);

protected:
  /** Raw communicator
   *
   *  Currently it is only owned by Communicate which manages its creation and destruction
   *  After switching to mpi3::communicator, myMPI is only a reference to the raw communicator owned by mpi3::communicator
   */
  mpi_comm_type myMPI;
  /// Communicator name
  std::string myName;
  /// Rank
  int d_mycontext;
  /// Size
  int d_ncontexts;
  /// Group ID of the current communicator in the parent communicator
  int d_groupid;
  /// Total number of groups in the parent communicator
  int d_ngroups;

public:
  /// Group Lead Communicator
  Communicate* GroupLeaderComm;

#ifdef HAVE_MPI
  /// mpi3 communicator wrapper
  mpi3::communicator comm;
#endif
};


namespace OHMMS
{
/** Global Communicator for a process
 */
extern Communicate* Controller;
} // namespace OHMMS


#endif // OHMMS_COMMUNICATE_H
