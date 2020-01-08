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
  InfoStream(std::ostream* output_stream) : currStream(output_stream), prevStream(nullptr), nullStream(nullptr)
  { }

  InfoStream(InfoStream& in) : currStream(&in.getStream()), prevStream(nullptr), nullStream(nullptr)
  { }

  std::ostream& getStream(const std::string& tag = "") { return *currStream; }

  void setStream(std::ostream* output_stream) { currStream = output_stream; }


  void flush() { currStream->flush(); }

  /// Stop output (redirect to a null stream)
  void pause();

  ///  Continue output on the stream used before pausing
  void resume();

  /// Open a file and output to that file
  void redirectToFile(const std::string& fname);

  /// Copy a stream
  void redirectToSameStream(InfoStream& info);

  ///  Permanently turn off the stream
  void shutOff();

private:
  // ofstream pointer allocated when file output is used
  std::unique_ptr<std::ofstream> fileStream;

  // current stream that InfoStream represents
  std::ostream* currStream;

  // save stream during pause
  std::ostream* prevStream;

  // Created at construction. Used during pause
  std::ostream nullStream;
};

template<class T>
inline InfoStream& operator<<(InfoStream& o, const T& val)
{
  o.getStream() << val;
  return o;
}


#endif
