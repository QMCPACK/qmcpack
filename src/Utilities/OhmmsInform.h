//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file OhmmsInform.h
 * @brief Declaration of OhmmsInform class.
 */
#ifndef OHMMSINFORM_H
#define OHMMSINFORM_H

#include <string>
#include <iostream>
#include <iomanip>
#if (__GNUC__ == 2)
#include <strstream.h>
#else
#include <sstream>
#endif
#include <strstream>

/**  Wrapper of std::ostream to provide uniform environments for run-time report
 *
 * Each object can choose a prompt, a mode (only master processor or all the processors),
 * and an std::ostream.
 */
class OhmmsInform
{

public:

  enum {OVERWRITE, APPEND};
  OhmmsInform(bool allcanwrite=true, bool writenode=true);
  OhmmsInform(const char* prompt, bool allcanwrite=true, bool writenode=true);
  OhmmsInform(const char* prompt, const char* fname, int appmode = OVERWRITE);
  OhmmsInform(const char* prompt, std::ostream& );

  ~OhmmsInform();

  void set(const char* fname, int appmode = OVERWRITE);
  void set(OhmmsInform&);
  void set(OhmmsInform&,const std::string& s);
  void setPrompt(const std::string& s);
  inline void flush()
  {
    myStream->flush();
  }

  inline std::ostream& getStream()
  {
    return *myStream;
  }

  ///switch to std::cerr
  void setStdError();
  ///incremenet indentation
  inline void pushd()
  {
    Blanks+=2;
  }
  ///decrease indentation
  inline void popd()
  {
    if(Blanks)
      Blanks-=2;
  }

  //inline bool open() const { return CanWrite;}

  /** temporarily turn off the stream
   */
  void turnoff();

  /** reset the stream to the original
   */
  void reset();

private:
  bool OwnStream;
  int Blanks;
  std::ostream* myStream;
  std::ostream* bgStream;
  std::string   myPrompt;
};

// templated version of operator<< for Inform objects
template<class T>
inline
OhmmsInform& operator<<(OhmmsInform& o, const T& val)
{
  o.getStream() << val;
  return o;
}

// specialized function for sending strings to Inform object
inline OhmmsInform& operator<<(OhmmsInform& o, const std::string& s)
{
  o.getStream() << s.c_str();
  return o;
}


#endif//OHMMSINFO_H

