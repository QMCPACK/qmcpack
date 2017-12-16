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
    
    

#ifndef OHMMS_BOOSTRANDOM_H
#define OHMMS_BOOSTRANDOM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <ctime>
#include <sstream>
#include <limits>
#include <boost/random.hpp>

/** random number generator using boost::random
 *
 * A wrapper of boost::random class to work with applicatoins.
 */
template<typename T, typename RNG=boost::mt19937>
class BoostRandom
{

public:
  /// real result type
  typedef T result_type;
  /// randmon number generator [0,max) where max depends on the generator type
  typedef RNG generator_type;
  /// unsigned integer type
  typedef typename generator_type::result_type uint_type;
  /// real random generator [0,1)
  typedef boost::variate_generator<generator_type, boost::uniform_real<T> > uniform_generator_type;

  std::string ClassName;
  std::string EngineName;

  ///default constructor
  explicit BoostRandom(uint_type iseed=911, const std::string& aname="mt19937")
    : ClassName("boost"), EngineName(aname),
      myContext(0), nContexts(1), baseOffset(0),
      uni(generator_type(iseed),boost::uniform_real<T>(0,1))
  {
  }

  ///copy constructor
  BoostRandom(const BoostRandom& rng): ClassName(rng.ClassName), EngineName(rng.EngineName),
    myContext(rng.myContext), nContexts(rng.nContexts), baseOffset(rng.baseOffset),
    uni(rng.uni)
  {
  }

  ///copy operator (unnecessary but why not)
  BoostRandom<T,RNG>& operator=(const BoostRandom& r)
  {
    ClassName=r.ClassName;
    EngineName=r.EngineName;
    myContext=r.myContext;
    nContexts=r.nContexts;
    baseOffset=r.baseOffset,
    uni=r.uni;
    return *this;
  }

  ~BoostRandom() { }

  /** initialize the generator
   * @param i thread index
   * @param nstr number of threads
   * @param iseed_in input seed
   *
   * Initialize generator with the seed.
   */
  void init(int i, int nstr, int iseed_in, uint_type offset=1)
  {
    uint_type baseSeed=iseed_in;
    myContext=i;
    nContexts=nstr;
    if(iseed_in<=0)
      baseSeed=make_seed(i,nstr);
    baseOffset=offset;
    uni.engine().seed(baseSeed);
  }

  ///get baseOffset
  inline int offset() const
  {
    return baseOffset;
  }
  ///assign baseOffset
  inline int& offset()
  {
    return baseOffset;
  }

  ///assign seed
  inline void seed(uint_type aseed)
  {
    uni.engine().seed(aseed);
  }

  uniform_generator_type& engine()
  {
    return uni;
  }
  /////reset the seed
  //inline void reset()
  //{
  //  uni.engine().seed(make_seed(myContext,nContexts,baseOffset));
  //}

  /** return a random number [0,1)
   */
  inline result_type rand()
  {
    return uni();
  }

  /** return a random number [0,1)
   */
  inline result_type operator()()
  {
    return uni();
  }

  /** return a random integer
   */
  inline uint_type irand()
  {
    return uni.engine()()%std::numeric_limits<uint_type>::max();
  }

  /** generate a series of random numbers */
  template<typename T1>
  inline void generate_uniform(T1* restrict d, int n)
  {
    for(int i=0; i<n; ++i) d[i]=uni();
  }

  inline void generate_normal(T* restrict d, int n)
  {
    BoxMuller2::generate(*this,d,n);
  }

  //inline void bivariate(resul_type& g1, resul_type &g2) {
  //  resul_type v1, v2, r;
  //  do {
  //  v1 = 2.0e0*uni() - 1.0e0;
  //  v2 = 2.0e0*uni() - 1.0e0;
  //  r = v1*v1+v2*v2;
  //  } while(r > 1.0e0);
  //  resul_type fac = sqrt(-2.0e0*log(r)/r);
  //  g1 = v1*fac;
  //  g2 = v2*fac;
  //}

  inline int state_size() const
  {
    return uni.engine().state_size;
  }

  inline void read(std::istream& rin)
  {
    rin >> uni.engine();
  }

  inline void write(std::ostream& rout) const
  {
    rout << uni.engine();
  }

  inline void save(std::vector<uint_type>& curstate) const
  {
    curstate.clear();
    std::stringstream otemp;
    otemp << uni.engine();
    copy(std::istream_iterator<uint_type>(otemp)
              , std::istream_iterator<uint_type>(),back_inserter(curstate));
  }

  inline void load(const std::vector<uint_type>& newstate)
  {
    std::stringstream otemp;
    copy(newstate.begin(),newstate.end(),std::ostream_iterator<uint_type>(otemp," "));
    otemp >> uni.engine();
  }

private:
  ///context number
  int myContext;
  ///number of contexts
  int nContexts;
  ///offset of the random seed
  int baseOffset;
  ///random number generator [0,1)
  uniform_generator_type uni;
};
#endif

