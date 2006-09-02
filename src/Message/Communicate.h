//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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
#include "ohmms-config.h"
#endif

#ifdef HAVE_OOMPI
#include "oompi.h"
#else
  #ifdef HAVE_MPI
  #include <mpi.h>
  #else
  typedef int MPI_Status;
  typedef int MPI_Request;
  #endif
#endif

/**@class Communicate
 * @ingroup Message
 * @brief 
 *  Wrapping information on parallelism.
 *  Very limited in functions. Currently, only single-mode or mpi-mode 
 *  is available (mutually exclusive).
 * @todo Possibly, make it a general manager class for mpi+openmp, mpi+mpi
 */
class Communicate {
public:

#ifdef HAVE_MPI
   typedef MPI_Comm mpi_comm_type;
#else
   typedef int mpi_comm_type;
#endif

  ///constructor
  Communicate();

  ///constructor with arguments
  Communicate(int argc, char **argv);

  /**destructor
   * Call proper finalization of Communication library
   */
  virtual ~Communicate();

  void initialize(int argc, char **argv);
  void finalize();
  void abort();

  ///return the Communicator ID (typically MPI_WORLD_COMM)
  inline mpi_comm_type getID() const { return CommID;}


  ///return the rank of this node
  inline int getNodeID() const { return d_mycontext;}
  inline int mycontext() const { return d_mycontext;}

  ///return the number of nodes
  inline int getNumNodes() const { return d_ncontexts;}
  inline int ncontexts() const { return d_ncontexts;}

  inline bool master() const { return (d_mycontext == 0);}

  void cleanupMessage(void*);
  inline void setNodeID(int i) { d_mycontext = i;}
  inline void setNumNodes(int n) { d_ncontexts = n;}
  inline void setCommID(mpi_comm_type i) { CommID = i;}
  void barrier();

protected:

  mpi_comm_type CommID; 
  int d_mycontext; 
  int d_ncontexts;

};


namespace OHMMS {
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
