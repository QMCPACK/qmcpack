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

#include "Message/CommCreate.h"

// initialize static data members

int TagMaker::CurrentTag = 1000;

Communicate* CommCreate::Comm =0;

Communicate* CommCreate::get() {
  if(!Comm)  {
    Comm = new Communicate;
  }
  return Comm;
}

Communicate* CommCreate::get(int argc, char **argv) {
  if(!Comm)  {
    Comm = new Communicate(argc,argv);
  }
  return Comm;
}

void CommCreate::remove() {
  if(Comm) delete Comm;
  Comm = 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
