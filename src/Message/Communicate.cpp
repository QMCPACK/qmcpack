//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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

#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif
#include "Message/Communicate.h"
#include "Message/TagMaker.h"

//static data of TagMaker::CurrentTag is initialized.
int TagMaker::CurrentTag = 1000;

//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;

//default constructor: ready for a serial execution
Communicate::Communicate():d_mycontext(0), d_ncontexts(1), CommID(0){}

Communicate::Communicate(int argc, char **argv){
  initialize(argc,argv);
}

//exclusive:  OOMPI, MPI or Serial
#ifdef HAVE_OOMPI

//================================================================
// Implements Communicate with OOMPI library
//================================================================
#include "oompi.h"

Communicate::~Communicate(){ }

void Communicate::initialize(int argc, char **argv){
  OOMPI_COMM_WORLD.Init(argc, argv);
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();
  CommID = OOMPI_COMM_WORLD.Get_mpi();
  //LOGMSG("Initializating MPI/OOMPI communicator")
}

void Communicate::finalize() {
  OOMPI_COMM_WORLD.Finalize();
}

void Communicate::cleanupMessage(void*) { }

#else

#ifdef HAVE_MPI

//================================================================
// Implements Communicate with standard MPI library
//================================================================
#include <mpi.h>
Communicate::~Communicate(){ }

void Communicate::initialize(int argc, char **argv){
  int flag;
  MPI_Initialized(&flag);
  if(!flag) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&d_ncontexts);
    MPI_Comm_rank(MPI_COMM_WORLD,&d_mycontext);
    CommID = MPI_COMM_WORLD;
  }
}

void Communicate::finalize(){
  MPI_Finalize();
}

void Communicate::cleanupMessage(void*) { }

#else //HAVE_MPI

Communicate::~Communicate(){}
void Communicate::initialize(int argc, char **argv){ }
void Communicate::finalize(){ }
void Communicate::cleanupMessage(void*) { }

#endif // !HAVE_MPI

#endif // !HAVE_OOMPI
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
