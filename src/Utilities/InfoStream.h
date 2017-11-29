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

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

// Bit bucket to discard output
// From https://stackoverflow.com/questions/7818371/printing-to-nowhere-with-ostream
class NullStreamBuf : public std::streambuf
{
protected:
    virtual int overflow(int c) { return c; }
};

class NullStream : public NullStreamBuf, public std::ostream
{
public:
    NullStream() : std::ostream(this) {}
};

/**
 *  Interface to output streams.  Can redirect output to stdout/stderr, a file, or a null stream.
 */

class InfoStream
{
public:
  InfoStream(std::ostream *output_stream): prevStream(NULL), nullStream(new NullStream),
      ownStream(false) {
    currStream = output_stream;
  }

  ~InfoStream();

  std::ostream &getStream(const std::string &tag = "") {
    return *currStream;
  }

  void setStream(std::ostream *output_stream) {
    currStream = output_stream;
  }


  void flush() {
    currStream->flush();
  }

  /// Stop output (redirect to a null stream)
  void pause();

  ///  Continue output on the stream used before pausing
  void resume();

  /// Open a file and output to that file
  void redirectToFile(const std::string &fname);

  /// Copy a stream
  void redirectToSameStream(InfoStream &info);

private:

  // Keep track of whether we should delete the stream or not
  bool ownStream;

  std::ostream *currStream;

  // save stream during pause
  std::ostream *prevStream;

  // Created at construction. Used during pause
  std::ostream *nullStream;
};


#endif
