//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003- by Jeongnim Kim
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

#include "Message/Communicate.h"
#include "Message/TagMaker.h"
#include <iostream>

//static data of TagMaker::CurrentTag is initialized.
int TagMaker::CurrentTag = 1000;

//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;

//default constructor: ready for a serial execution
Communicate::Communicate():
d_mycontext(0), d_ncontexts(1), myCommID(0)
{
}

Communicate::Communicate(int argc, char **argv){
  initialize(argc,argv);
}

//exclusive:  OOMPI, MPI or Serial
#ifdef HAVE_OOMPI

Communicate::Communicate(const intra_comm_type& c):myComm(c) {
  myCommID = myComm.Get_mpi();
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
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();
  myComm = OOMPI_COMM_WORLD;
  myCommID = OOMPI_COMM_WORLD.Get_mpi();
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
  OOMPI_COMM_WORLD.Barrier();
}

void Communicate::abort(const char* msg)
{ 
  std::cerr << msg << std::endl;
  OOMPI_COMM_WORLD.Abort();
}

/** split a communicator into ng groups
 * @param number of groups
 * @return a communicator
 */
Communicate::intra_comm_type
Communicate::split(int ng) 
{
  //this is a workaround due to the OOMPI bug with split
  if(ng>1) {
    MPI_Comm row;
    int n=d_ncontexts/ng;
    int p=d_mycontext/n;
    int q=d_mycontext%n;
    MPI_Comm_split(myCommID,p,q,&row);
    return OOMPI_Intra_comm(row);
  } else {
    return OOMPI_Intra_comm(myCommID);
  }
}

#else

Communicate::~Communicate(){}
void Communicate::initialize(int argc, char **argv){ }
void Communicate::finalize(){ }
void Communicate::abort(){ 
  abort();
}

void Communicate::abort(const char* msg){ 
  std::cerr << msg << std::endl;
  abort();
}
void Communicate::barrier(){
}
void Communicate::cleanupMessage(void*) { }

Communicate::mpi_comm_type Communicate::split(int n) {
  return myCommID;
}

Communicate::Communicate(const intra_comm_type& c):
d_mycontext(0), d_ncontexts(1), myCommID(0)
{
}
#endif // !HAVE_OOMPI
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
