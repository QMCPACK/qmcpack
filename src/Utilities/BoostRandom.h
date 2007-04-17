//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_BOOSTRANDOM_H
#define OHMMS_BOOSTRANDOM_H

#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif

#include <ctime>        
#include <sstream>        
#include <boost/random.hpp>
#include "Message/Communicate.h"

template<class T>
class BoostRandom {

public:

  typedef T Return_t;
  typedef boost::minstd_rand      base_generator_type;
  typedef boost::mt19937          generator_type;
  typedef boost::variate_generator<generator_type,boost::uniform_real<T> > uniform_generator_type;

  std::string ClassName;
  std::string EngineName;

  explicit BoostRandom(int iseed=0): ClassName("boost"), EngineName("mt19937"), 
  baseSeed(iseed),  unit_dist(0,1), generator(0), uni(0) 
  { }

  ~BoostRandom() {
    if(uni) delete uni;
    if(generator) delete generator;
  }

  /** initialize the generator
   * @param i thread index
   * @param nstr number of threads
   * @param iseed_in input seed
   *
   * Initialize generator with the seed. 
   */
  void init(int i, int nstr, int iseed_in) 
  {
    baseSeed=iseed_in;
    myContext=i;
    nContexts=nstr;
    if(iseed_in<=0) 
      baseSeed=static_cast<uint32_t>(std::time(0))%16081+(i+1)*nstr+i;

    //use the constructor with the seed
    //generator->seed(baseSeed) is not working!!!
    if(generator == 0) 
    {
      generator = new generator_type(baseSeed);
      uni = new uniform_generator_type(*generator,unit_dist);
      std::cout << "  BoostRandom::init " << myContext << " " << baseSeed << std::endl;
    }
  }

  //randomize again
  void reset() 
  {
    delete uni;
    delete generator;
    baseSeed=static_cast<uint32_t>(std::time(0))%16081+(myContext+1)*nContexts+myContext;
    generator = new generator_type(baseSeed);
    uni = new uniform_generator_type(*generator,unit_dist);
    std::cout << "  BoostRandom::reset " << myContext << " " << baseSeed << std::endl;
  }

  inline uniform_generator_type& getGenerator() { return *uni;}

  //!< return [0,1)
  inline T getRandom() { return (*uni)(); }

  inline Return_t operator()() { return (*uni)();} 
  inline int irand() { return 1;}

  inline void bivariate(Return_t& g1, Return_t &g2) {
    Return_t v1, v2, r;
    do {
    v1 = 2.0e0*((*uni)()) - 1.0e0;
    v2 = 2.0e0*((*uni)()) - 1.0e0;
    r = v1*v1+v2*v2;
    } while(r > 1.0e0);
    Return_t fac = sqrt(-2.0e0*log(r)/r);
    g1 = v1*fac;
    g2 = v2*fac;
  }

  inline void read(std::istream& rin) {
    rin >> uni->engine();
  }

  inline void write(std::ostream& rout) const {
    rout << uni->engine();
  }

private:
  uint32_t baseSeed;
  int myContext;
  int nContexts;
  base_generator_type    base_generator;
  boost::uniform_real<T> unit_dist;
  generator_type         *generator;
  uniform_generator_type *uni;
};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
