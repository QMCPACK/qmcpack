//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/** Communicate class that wraps information on parallelism.
 * @ingroup Message
 *
 *  Very limited in functions. Currently, only single-mode or mpi-mode 
 *  is available (mutually exclusive).
 * @todo Possibly, make it a general manager class for mpi+openmp, mpi+mpi
 */
class Communicate {
public:

  ///constructor
  Communicate();

  ///constructor with arguments
  Communicate(int argc, char **argv);

  /**destructor
   *
   * Call proper finalization of Communication library
   */
  virtual ~Communicate();

  void initialize(int argc, char **argv);
  void finalize();

  ///return the Communicator ID (typically MPI_WORLD_COMM)
  inline int getID() const { return CommID;}


  ///return the rank of this node
  inline int getNodeID() const { return d_mycontext;}
  inline int mycontext() const { return d_mycontext;}

  ///return the number of nodes
  inline int getNumNodes() const { return d_ncontexts;}
  inline int ncontexts() const { return d_ncontexts;}

  inline bool master() const { return (d_mycontext == 0);}

  void cleanupMessage(void*);
  inline void setNodeID(int i) { d_mycontext = i;}
  inline void setCommID(int i) { CommID = i;}

  void barrier();

protected:

  int CommID; 
  int d_mycontext; 
  int d_ncontexts;

};


namespace OHMMS {
  /** Global Communicator for a process */
  extern Communicate* Controller;
}
#endif // OHMMS_COMMUNICATE_H
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
