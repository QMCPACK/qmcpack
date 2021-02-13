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


#include "InfoStream.h"
#include <fstream>
#include <cassert>

InfoStream::InfoStream(std::ostream* output_stream)
    : currStream(output_stream), outputStream(output_stream), nullStream(nullptr)
{
  if (output_stream == nullptr)
    currStream = outputStream = &nullStream;
}

void InfoStream::setStream(std::ostream* output_stream)
{
  if (outputStream == output_stream)
    return;
  assert(checkCurr());
  if (output_stream == nullptr)
    outputStream = &nullStream;
  else
    outputStream = output_stream;
  if (currStream != &nullStream)
    currStream = outputStream;
  fileStream.reset();
  assert(checkCurr());
}

void InfoStream::flush() { outputStream->flush(); }

void InfoStream::pause()
{
  assert(checkCurr());
  currStream = &nullStream;
}

void InfoStream::resume()
{
  assert(checkCurr());
  currStream = outputStream;
}

void InfoStream::shutOff()
{
  assert(checkCurr());
  currStream = outputStream = &nullStream;
  fileStream.reset();
}

void InfoStream::redirectToFile(const std::string& fname)
{
  assert(checkCurr());
  fileStream   = std::make_unique<std::ofstream>(fname);
  outputStream = fileStream.get();
  if (currStream != &nullStream)
    currStream = outputStream;
  assert(checkCurr());
}

void InfoStream::redirectToSameStream(InfoStream& info)
{
  if (this == &info)
    return;
  assert(checkCurr());
  outputStream = info.outputStream;
  if (currStream != &nullStream)
    currStream = outputStream;
  fileStream.reset();
  assert(checkCurr());
}
