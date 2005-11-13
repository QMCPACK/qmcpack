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
/** @file OhmmsInform.cpp
 * @brief Definition of OhmmsInform class.
 */
#include "Utilities/OhmmsInform.h"

#include <iomanip>
#include <fstream>
using namespace std;

OhmmsInform::OhmmsInform(bool allcanwrite, bool writenode) 
:thisStream(0), OwnStream(false) 
{ 
  thisPrompt = string("ohmms>");
  if(allcanwrite || writenode) { 
    thisStream = &cout;
    CanWrite = true;
  } else {
    thisStream = new ostringstream();
    CanWrite = false;
  }
}

OhmmsInform::OhmmsInform(const char* prompt, bool allcanwrite, bool writenode)
:thisStream(0), thisPrompt(prompt), OwnStream(false)
{ 
  if(allcanwrite || writenode) { 
    thisStream = &cout;
    CanWrite = true;
  } else {
    thisStream = new ostringstream();
    CanWrite = false;
  }
}

OhmmsInform::OhmmsInform(const char* prompt, const char* fname, int appmode)
  :thisPrompt(prompt), OwnStream(true), CanWrite(true) 
{ 
  // file mode
  if(appmode == OVERWRITE) 
    thisStream = new ofstream(fname);
  else
    thisStream = new ofstream(fname,ios::app);
}


OhmmsInform::OhmmsInform(const char* prompt, ostream& o)
  :thisPrompt(prompt), thisStream(&o), CanWrite(true)
{
  OwnStream = false;
}


void OhmmsInform::set(const char* fname, int appmode)
{
  if(OwnStream && thisStream) delete thisStream;
  OwnStream = true;
  if(appmode == OVERWRITE)
    thisStream = new ofstream(fname);
  else
    thisStream = new ofstream(fname,ios::app);
}

void OhmmsInform::set(OhmmsInform& o)
{
  if(OwnStream) {
    if(thisStream) delete thisStream;
  }
  thisStream = o.thisStream;
  OwnStream = false;
}

OhmmsInform::~OhmmsInform() {
  if(OwnStream) delete thisStream;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
