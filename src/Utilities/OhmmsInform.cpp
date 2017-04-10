//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file OhmmsInform.cpp
 * @brief Definition of OhmmsInform class.
 */
#include "Utilities/OhmmsInform.h"

#include <iomanip>
#include <fstream>

OhmmsInform::OhmmsInform(bool allcanwrite, bool writenode):bgStream(0),myPrompt("qmc>")
{
  if(allcanwrite || writenode)
  {
    myStream = &std::cout;
    OwnStream =false;
  }
  else
  {
    myStream = new std::ostringstream();
    OwnStream=true;
  }
  Blanks=0;
}

OhmmsInform::OhmmsInform(const char* prompt, bool allcanwrite, bool writenode):myStream(0),bgStream(0),myPrompt(prompt)
{
  if(allcanwrite || writenode)
  {
    myStream = &std::cout;
    OwnStream=false;
  }
  else
  {
    myStream = new std::ostringstream();
    OwnStream=true;
  }
  Blanks=0;
}

OhmmsInform::OhmmsInform(const char* prompt, const char* fname, int appmode):OwnStream(true), bgStream(0),myPrompt(prompt)
{
  // file mode
  if(appmode == OVERWRITE)
    myStream = new std::ofstream(fname);
  else
    myStream = new std::ofstream(fname,std::ios::app);
  Blanks=0;
}


OhmmsInform::OhmmsInform(const char* prompt, std::ostream& o)
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
    myStream = new std::ofstream(fname);
  else
    myStream = new std::ofstream(fname,std::ios::app);
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

void OhmmsInform::set(OhmmsInform& o, const std::string& s)
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
void OhmmsInform::setPrompt(const std::string& s)
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
  myStream=&std::cerr;
  OwnStream=false;
}

void OhmmsInform::turnoff()
{
  bgStream=myStream;
  myStream = new std::ostream(0);
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

