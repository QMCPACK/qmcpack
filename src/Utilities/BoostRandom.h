//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_BOOSTRANDOM_H
#define OHMMS_BOOSTRANDOM_H

#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif

#include <ctime>        
#include <boost/random.hpp>
#include "Message/Communicate.h"

template<class T>
class BoostRandom {

public:

  typedef T Return_t;
  typedef boost::minstd_rand base_generator_type;
  typedef boost::mt19937     generator_type;
  typedef boost::variate_generator<generator_type,boost::uniform_real<T> > uniform_generator_type;

  BoostRandom(): baseSeed(0), unit_dist(0,1), uni(generator,unit_dist) { 
    init(OHMMS::Controller->mycontext(),OHMMS::Controller->ncontexts(),-1);
  }

  BoostRandom(int iseed): 
    baseSeed(iseed), unit_dist(0,1), uni(generator,unit_dist)
  {
    init(OHMMS::Controller->mycontext(),OHMMS::Controller->ncontexts(),iseed);
  }
 
  void init(int i, int nstr, int iseed) {
    if(iseed<0) {
      iseed=static_cast<unsigned int>(std::time(0))%16081+(i+1)*nstr+i;
    } 
    baseSeed=iseed;
    base_generator.seed(iseed);
    generator.seed(base_generator);
  }

  inline uniform_generator_type& getGenerator() { return uni;}

  //!< return [0,1)
  inline T getRandom() { return uni(); }

  inline Return_t operator()() { return uni();} 
  inline int irand() { return 1;}

  inline void bivariate(Return_t& g1, Return_t &g2) {
    Return_t v1, v2, r;
    do {
    v1 = 2.0e0*uni() - 1.0e0;
    v2 = 2.0e0*uni() - 1.0e0;
    r = v1*v1+v2*v2;
    } while(r > 1.0e0);
    Return_t fac = sqrt(-2.0e0*log(r)/r);
    g1 = v1*fac;
    g2 = v2*fac;
  }

private:
  unsigned int baseSeed;
  base_generator_type    base_generator;
  generator_type         generator;
  boost::uniform_real<T> unit_dist;
  uniform_generator_type uni;
};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
