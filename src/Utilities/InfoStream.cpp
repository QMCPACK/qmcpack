//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <Utilities/InfoStream.h>
#include <fstream>

InfoStream::~InfoStream()
{
  if (currStream != nullStream) {
    delete nullStream;
  }

  if (ownStream && currStream) {
    delete(currStream);
  }
}

void InfoStream::pause() {
  if (currStream != nullStream)
  {
    prevStream = currStream;
    currStream = nullStream;
  }
}

void InfoStream::resume() {
  if (prevStream) {
    currStream = prevStream;
    prevStream = NULL;
   }
}

void InfoStream::shutOff() {
  prevStream = NULL;
  currStream = nullStream;
}

void InfoStream::redirectToFile(const std::string &fname) {
  currStream = new std::ofstream(fname);
  ownStream = true;
}

void InfoStream::redirectToSameStream(InfoStream &info)
{
  currStream = &info.getStream();
  ownStream = false;
}

