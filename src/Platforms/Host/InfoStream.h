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


/** @file InfoStream.h
 * @brief Declaration of InfoStream class.
 */

#ifndef INFOSTREAM_H
#define INFOSTREAM_H

#include <fstream>
#include <memory>
#include <iostream>

/**
 *  Interface to output streams.  Can redirect output to stdout/stderr, a file, or a null stream.
 */

class InfoStream
{
public:
  InfoStream(std::ostream* output_stream);

  InfoStream(const InfoStream& in) = delete;
  InfoStream& operator=(const InfoStream& in) = delete;

  /// returns current stream
  std::ostream& getStream(const std::string& tag = "") { return *currStream; }

  void setStream(std::ostream* output_stream);

  /// flush stream buffer
  void flush();

  /// Stop output (redirect to a null stream)
  void pause();

  ///  Continue output on the stream used before pausing
  void resume();

  /// check if the stream is active
  bool isActive() const { return currStream != &nullStream; }

  /// Open a file and output to that file
  void redirectToFile(const std::string& fname);

  /// Copy a stream
  void redirectToSameStream(InfoStream& info);

  ///  Permanently turn off the stream
  void shutOff();

private:
  // check currStream, it can only be outputStream or &nullStream
  bool checkCurr() const { return currStream == outputStream || currStream == &nullStream; }

  // current stream that InfoStream represents
  std::ostream* currStream;

  // the output stream when InfoStream IsActive. (1. from external, 2. nullStream (off), 3. fileStream)
  std::ostream* outputStream;

  // Created at construction. Used during pause and shutoff
  std::ostream nullStream;

  // ofstream pointer allocated when file output is used
  std::unique_ptr<std::ofstream> fileStream;
};

template<class T>
inline InfoStream& operator<<(InfoStream& o, const T& val)
{
  o.getStream() << val;
  return o;
}


#endif
