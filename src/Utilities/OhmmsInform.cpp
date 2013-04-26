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

OhmmsInform::OhmmsInform(bool allcanwrite, bool writenode):bgStream(0),myPrompt("qmc>")
{
  if(allcanwrite || writenode)
  {
    myStream = &cout;
    OwnStream =false;
  }
  else
  {
    myStream = new ostringstream();
    OwnStream=true;
  }
  Blanks=0;
}

OhmmsInform::OhmmsInform(const char* prompt, bool allcanwrite, bool writenode):myStream(0),bgStream(0),myPrompt(prompt)
{
  if(allcanwrite || writenode)
  {
    myStream = &cout;
    OwnStream=false;
  }
  else
  {
    myStream = new ostringstream();
    OwnStream=true;
  }
  Blanks=0;
}

OhmmsInform::OhmmsInform(const char* prompt, const char* fname, int appmode):OwnStream(true), bgStream(0),myPrompt(prompt)
{
  // file mode
  if(appmode == OVERWRITE)
    myStream = new ofstream(fname);
  else
    myStream = new ofstream(fname,ios::app);
  Blanks=0;
}


OhmmsInform::OhmmsInform(const char* prompt, ostream& o)
  :OwnStream(false),myStream(&o),myPrompt(prompt)
{
  Blanks=0;
}


void OhmmsInform::set(const char* fname, int appmode)
{
  if(OwnStream && myStream)
    delete myStream;
  OwnStream = true;
  if(appmode == OVERWRITE)
    myStream = new ofstream(fname);
  else
    myStream = new ofstream(fname,ios::app);
}


void OhmmsInform::set(OhmmsInform& o)
{
  if(OwnStream)
  {
    if(myStream)
      delete myStream;
  }
  myStream =o.myStream;
  myPrompt=o.myPrompt;
  OwnStream = false;
}

void OhmmsInform::set(OhmmsInform& o, const string& s)
{
  if(OwnStream)
  {
    if(myStream)
      delete myStream;
  }
  myStream =o.myStream;
  myPrompt=s;
  OwnStream = false;
}
void OhmmsInform::setPrompt(const string& s)
{
  myPrompt=s;
}

void OhmmsInform::setStdError()
{
  if(OwnStream)
  {
    if(myStream)
      delete myStream;
  }
  myStream=&cerr;
  OwnStream=false;
}

void OhmmsInform::turnoff()
{
  bgStream=myStream;
  myStream = new ostringstream();
}

void OhmmsInform::reset()
{
  if(bgStream)
  {
    delete myStream;
    myStream=bgStream;
    bgStream=0;
  }
}

OhmmsInform::~OhmmsInform()
{
  if(OwnStream)
    delete myStream;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
