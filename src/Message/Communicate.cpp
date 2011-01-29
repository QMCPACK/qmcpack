/////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003- by Jeongnim Kim
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

#include "Message/Communicate.h"
#include "Message/TagMaker.h"
#include <iostream>
#include <cstdio>
#include <Platforms/sysutil.h>
#include <tau/profiler.h>

//static data of TagMaker::CurrentTag is initialized.
int TagMaker::CurrentTag = 1000;

//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;

//default constructor: ready for a serial execution
Communicate::Communicate():
myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0)
{
}


Communicate::Communicate(int argc, char **argv)
{
  initialize(argc,argv);
}

//exclusive:  OOMPI, MPI or Serial
#ifdef HAVE_OOMPI

void
Communicate::set_world()
{
  myComm = OOMPI_COMM_WORLD;
  myMPI = myComm.Get_mpi();
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();
}


Communicate::Communicate(const Communicate& comm, int nparts)
{
  //this is a workaround due to the OOMPI bug with split
  if(nparts>1) {
    MPI_Comm row;
    int n=comm.size()/nparts;
    int p=comm.rank()/n;
    int q=comm.rank()%n;
    MPI_Comm_split(comm.getMPI(),p,q,&row);
    myComm=OOMPI_Intra_comm(row);
    d_groupid=p;
  } else {
    myComm=OOMPI_Intra_comm(comm.getComm());
    d_groupid=0;
  }

  myMPI = myComm.Get_mpi();
  d_mycontext=myComm.Rank();
  d_ncontexts=myComm.Size();
}


//================================================================
// Implements Communicate with OOMPI library
//================================================================
Communicate::~Communicate()
{ 

}

void Communicate::initialize(int argc, char **argv)
{

  OOMPI_COMM_WORLD.Init(argc, argv);
  myComm = OOMPI_COMM_WORLD;
  myMPI = myComm.Get_mpi();
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();

#ifdef __linux__
  for (int proc=0; proc<OHMMS::Controller->size(); proc++) {
    if (OHMMS::Controller->rank() == proc) {
      fprintf (stderr, "Rank = %4d  Free Memory = %5ld MB\n", proc, freemem());
    }
    barrier();
  }
  barrier();
#endif

  std::string when="qmc."+getDateAndTime("%Y%m%d_%H%M");
  hpmInit(QMC_MAIN_EVENT,when.c_str());
}

void Communicate::finalize() 
{
  OOMPI_COMM_WORLD.Finalize();
}

void Communicate::cleanupMessage(void*) 
{ 
}

void Communicate::abort()
{
  OOMPI_COMM_WORLD.Abort();
}

void Communicate::barrier()
{
  myComm.Barrier();
}

void Communicate::abort(const char* msg)
{ 
  std::cerr << msg << std::endl;
  OOMPI_COMM_WORLD.Abort();
}


#else

Communicate::~Communicate(){}

void Communicate::initialize(int argc, char **argv)
{ 
  std::string when="qmc."+getDateAndTime("%Y%m%d_%H%M");
  hpmInit(QMC_MAIN_EVENT,when.c_str());
}

void Communicate::set_world()
{ 
}

void Communicate::finalize()
{ 
  hpmTerminate(QMC_MAIN_EVENT);
}

void Communicate::abort()
{ 
  abort();
}

void Communicate::abort(const char* msg){ 
  std::cerr << msg << std::endl;
  abort();
}

void Communicate::barrier()
{
}

void Communicate::cleanupMessage(void*) 
{ 
}

Communicate::Communicate(const Communicate& comm, int nparts):
myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0)
{
}

#endif // !HAVE_OOMPI
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
