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

#ifndef OHMMS_COMMCREATE_H
#define OHMMS_COMMCREATE_H

#include "Message/TagMaker.h"
#include "Message/Communicate.h"
#ifdef USE_MPI
#include <mpi.h>
#endif 

/*!\class CommCreate
 * \brief Singleton Pattern to create one instance of Communicate
 * during a run. 
 */
class CommCreate {
public:

  static Communicate* get();
  static Communicate* get(int, char**);
  static void remove();

private:
  // factory
  CommCreate() { }
  ~CommCreate() { }
  static Communicate* Comm;
};
#endif // OHMMS_COMMCREATE_H
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
