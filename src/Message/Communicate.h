//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef OHMMS_COMMUNICATE_H
#define OHMMS_COMMUNICATE_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_OOMPI
#include "oompi.h"
struct CommunicatorTraits
{
  typedef MPI_Comm         mpi_comm_type;
  typedef MPI_Status       status;
  typedef MPI_Request      request;
  typedef OOMPI_Intra_comm intra_comm_type;
};
#define APP_ABORT(msg) \
{std::cerr << "Fatal Error. Aborting at " << msg << endl; OOMPI_COMM_WORLD.Abort();}

#else
struct CommunicatorTraits
{
  typedef int  mpi_comm_type;
  typedef int  status;
  typedef int  request;
  typedef int  intra_comm_type;
};

#define APP_ABORT(msg) \
{std::cerr << "Fatal Error. Aborting at " << msg << std::endl; std::cerr.flush(); abort();}

#endif

#include <string>
#include <vector>
#include <utility>


/**@class Communicate
 * @ingroup Message
 * @brief
 *  Wrapping information on parallelism.
 *  Very limited in functions. Currently, only single-mode or mpi-mode
 *  is available (mutually exclusive).
 * @todo Possibly, make it a general manager class for mpi+openmp, mpi+mpi
 */
class Communicate: public CommunicatorTraits
{

public:

  ///constructor
  Communicate();

  ///constructor with arguments
  Communicate(int argc, char **argv);

  ///constructor with communicator
  //Communicate(const intra_comm_type& c);
  Communicate(const Communicate& comm, int nparts);

  /** constructor
   * @param comm a communicator which will be split into groups
   * @param jobs number of tasks per group
   */
  Communicate(const Communicate& comm, const std::vector<int>& jobs);

  /**destructor
   * Call proper finalization of Communication library
   */
  virtual ~Communicate();

  void initialize(int argc, char **argv);
  void finalize();
  void abort();
  void abort(const char* msg);
  void set_world();

#if defined(HAVE_MPI)
  ///operator for implicit conversion to MPI_Comm
  inline operator MPI_Comm() const
  {
    return myMPI;
  }
#endif

  ///return the Communicator ID (typically MPI_WORLD_COMM)
  inline mpi_comm_type getMPI() const
  {
    return myMPI;
  }

  inline intra_comm_type& getComm()
  {
    return myComm;
  }
  inline const intra_comm_type& getComm() const
  {
    return myComm;
  }

  ///return the rank
  inline int rank() const
  {
    return d_mycontext;
  }
  ///return the number of tasks
  inline int size() const
  {
    return d_ncontexts;
  }

  ///return the group id
  inline int getGroupID() const
  {
    return d_groupid;
  }
  ///return the number of intra_comms which belong to the same group
  inline int getGroupSize() const
  {
    return d_ngroups;
  }
  //inline bool master() const { return (d_mycontext == 0);}
  //intra_comm_type split(int n);
  void cleanupMessage(void*);
  inline void setNodeID(int i)
  {
    d_mycontext = i;
  }
  inline void setNumNodes(int n)
  {
    d_ncontexts = n;
  }

  void barrier();

  inline void setName(const std::string& aname)
  {
    myName=aname;
  }
  inline const std::string& getName() const
  {
    return myName;
  }

  template<typename T> void allreduce(T&);
  template<typename T> void reduce(T&);
  template<typename T> void reduce(T* restrict, T* restrict, int n);
  template<typename T> void bcast(T&);
  template<typename T> void bcast(T* restrict, int n);
  template<typename T> void send(int dest, int tag, T&);
  template<typename T> void gather(T& sb, T& rb, int dest=0);
  template<typename T, typename IT> void gatherv(T& sb, T& rb, IT& counts, IT& displ, int dest=0);
  template<typename T> void allgather(T& sb, T& rb, int count);
  template<typename T, typename IT> void allgatherv(T& sb, T& rb, IT& counts, IT& displ);
  template<typename T> void scatter(T& sb, T& rb, int dest=0);
  template<typename T, typename IT> void scatterv(T& sb, T& rb, IT& counts, IT& displ, int source=0);
  template<typename T> request irecv(int source, int tag, T&);
  template<typename T> request isend(int dest, int tag, T&);
  template<typename T> request irecv(int source, int tag, T*, int n);
  template<typename T> request isend(int dest, int tag, T*, int n);

protected:

  mpi_comm_type myMPI;
  intra_comm_type myComm;
  std::string myName;
  int d_mycontext;
  int d_ncontexts;
  int d_groupid;
  int d_ngroups;
};


namespace OHMMS
{
/** Global Communicator for a process
 */
extern Communicate* Controller;
}


#endif // OHMMS_COMMUNICATE_H
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
