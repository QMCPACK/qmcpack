//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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

