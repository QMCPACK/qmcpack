// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called POOMA (Parallel Object-Oriented Methods and Applications) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-98-65.
// 
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No. W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without 
// charge, provided that this Notice and any statement of authorship are 
// reproduced on all copies.  Neither the Government nor the University 
// makes any warranty, express or implied, or assumes any liability or 
// responsibility for the use of this SOFTWARE.
// 
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
// 
// For more information about POOMA, send e-mail to pooma@acl.lanl.gov,
// or visit the POOMA web page at http://www.acl.lanl.gov/pooma/.
// ----------------------------------------------------------------------
// ACL:license

//
// Fri Apr 26 18:26:12 CDT 2002 Modification by Jeongnim Kim jnkim@uiuc.edu
//
#ifndef POOMA_UTILITIES_CLOCK_H
#define POOMA_UTILITIES_CLOCK_H

#include <time.h>

namespace Pooma {

//-----------------------------------------------------------------------
// Clock provides a running timer, utilizing high-speed SGI timers if 
// available.
//-----------------------------------------------------------------------

class Clock {
public:

  //---------------------------------------------------------------------
  // Set a static const that tells whether or not this class is utilizing
  // high-speed timers.
  
#if defined(CLOCK_SGI_CYCLE)
  static const bool highSpeed = true;
#else
  static const bool highSpeed = false;
#endif  

  //---------------------------------------------------------------------
  // Return the current value of the timer [sec].
  
  inline static double value()
  {
#if defined(CLOCK_SGI_CYCLE)
  timespec ts;
  clock_gettime(CLOCK_SGI_CYCLE, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
#else
  return double(clock()) / CLOCKS_PER_SEC;
#endif
  }

  //////////////////////////////////////////////////
  //functions to enable multiple clock objects
  // added by Jeongnim Kim
  //////////////////////////////////////////////////
  inline void start() { start_time = value();}
  inline void stop() { stop_time = value();}
  inline double cpu_time() const { return stop_time - start_time;}
  //////////////////////////////////////////////////

private:
  //////////////////////////////////////////////////
  //data members to enable multiple clock objects
  // added by Jeongnim Kim
  //////////////////////////////////////////////////
  double start_time, stop_time;

};

} // namespace Pooma


//////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_CLOCK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile$   $Author$
// $Revision$   $Date$
// ----------------------------------------------------------------------
// ACL:rcsinfo
