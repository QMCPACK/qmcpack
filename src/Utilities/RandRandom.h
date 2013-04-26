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
#ifndef OHMMS_RANDRANDOM_H
#define OHMMS_RANDRANDOM_H
#include <math.h>
#include <stdlib.h>

/*! \class RandRandom
 *  A wrapper class to generate a random number using C srand and rand.
 */
class RandRandom
{

public:
  typedef double Return_t;


  RandRandom(): thisStreamID(0), nStreams(1) { }

  RandRandom(int i, int nstr, int iseed)
  {
    init(i,nstr,iseed);
  }

  ~RandRandom() { }

  void init(int i = 0, int nstr =1, int iseed=-1);

  inline Return_t getRandom()
  {
    return Return_t(rand())*rand_max_inv;
  }

  inline Return_t operator()()
  {
    return getRandom();
  }

  inline int irand()
  {
    return rand();
  }

  inline void bivariate(Return_t& g1, Return_t& g2)
  {
    Return_t v1, v2, r;
    do
    {
      v1 = 2.0e0*getRandom() - 1.0e0;
      v2 = 2.0e0*getRandom() - 1.0e0;
      r = v1*v1+v2*v2;
    }
    while(r > 1.0e0);
    Return_t fac = sqrt(-2.0e0*log(r)/r);
    g1 = v1*fac;
    g2 = v2*fac;
  }

private:
  static const Return_t rand_max_inv;
  int thisStreamID;
  int nStreams;
  int thisSeed;
};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
