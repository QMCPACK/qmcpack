//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


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

void RandRandom::init(int i, int nstr, int iseed)
{
  thisStreamID = i;
  nStreams = nstr;
  if(iseed < 0)
  {
    struct_t timeval tvbuf;//Values from call to gettimeofday
    struct_t timezone tzbuf; //Timezone
    gettimeofday(&tvbuf, &tzbuf);
    unsigned long last_secs = tvbuf.tv_sec;
    thisSeed = last_secs%16081+thisStreamID*nStreams;
  }
  else
    if(iseed == 0)
    {
      thisSeed = (thisStreamID+1)*nStreams;
    }
    else
    {
      thisSeed = iseed;//(thisStreamID+1)*nStreams;
    }
  srand(thisSeed);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
