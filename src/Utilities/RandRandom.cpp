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

#include "Utilities/RandRandom.h"
extern "C" {
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
}

#ifdef SP2
#define struct_t
#else
#define struct_t struct
#endif

const double RandRandom::rand_max_inv = 1.0/(double(RAND_MAX) + 1);

void RandRandom::init(int i, int nstr, int iseed) {

  thisStreamID = i;
  nStreams = nstr;
  if(iseed < 0) {
    struct_t timeval tvbuf;//Values from call to gettimeofday
    struct_t timezone tzbuf; //Timezone
    gettimeofday(&tvbuf, &tzbuf);
    unsigned long last_secs = tvbuf.tv_sec;
    thisSeed = last_secs%16081+thisStreamID*nStreams;
  } else if(iseed == 0) {
    thisSeed = (thisStreamID+1)*nStreams;
  } else {
    thisSeed = iseed;//(thisStreamID+1)*nStreams;
  }
  srand(thisSeed);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
