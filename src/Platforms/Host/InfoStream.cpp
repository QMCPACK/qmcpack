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


#include <Platforms/Host/InfoStream.h>
#include <fstream>

void InfoStream::pause()
{
  if (currStream != &nullStream)
  {
    prevStream = currStream;
    currStream = &nullStream;
  }
}

void InfoStream::resume()
{
  if (prevStream)
  {
    currStream = prevStream;
    prevStream = nullptr;
  }
}

void InfoStream::shutOff()
{
  prevStream = nullptr;
  currStream = &nullStream;
}

void InfoStream::redirectToFile(const std::string& fname)
{
  fileStream = std::unique_ptr<std::ofstream>(new std::ofstream(fname));
  currStream = fileStream.get();
}

void InfoStream::redirectToSameStream(InfoStream& info)
{
  currStream = &info.getStream();
  fileStream.reset();
}
