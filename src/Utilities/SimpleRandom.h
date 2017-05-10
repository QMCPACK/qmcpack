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
    
    


#ifndef QMCPLUSPLUS_SLOWMT_RANDOM_H
#define QMCPLUSPLUS_SLOWMT_RANDOM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <ctime>
#include <vector>
#include <Utilities/MersenneTwister.h>

/** random number generator using boost::random
 *
 * A wrapper of boost::random class to work with applicatoins.
 */
template<typename ENG>
class SimpleRandom
{

public:
  /// randmon number generator [0,max) where max depends on the generator type
  typedef ENG generator_type;
  /// unsigned integer type
  typedef typename generator_type::uint_type uint_type;
  /// real result type
  typedef double result_type;

  std::string ClassName;
  std::string EngineName;

  ///default constructor
  explicit SimpleRandom(uint_type iseed=911, const std::string& aname="mtrand")
    : ClassName("mtrand"), EngineName(aname),
      myContext(0), nContexts(1), baseOffset(0), uni(iseed)
  {
  }


  inline ~SimpleRandom() { }

  /** initialize the generator
   * @param i thread index
   * @param nstr number of threads
   * @param iseed_in input seed
   *
   * Initialize generator with the seed.
   */
  inline void init(int i, int nstr, int iseed_in, uint_type offset=1)
  {
    uint_type baseSeed=iseed_in;
    myContext=i;
    nContexts=nstr;
    if(iseed_in<=0)
      baseSeed=make_seed(i,nstr);
    baseOffset=offset;
    uni.seed(baseSeed);
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

  /** return a random number [0,1)
   */
  inline result_type rand()
  {
    return uni.rand();
  }

  /** return a random number [0,1)
   */
  inline result_type operator()()
  {
    return uni.rand();
  }

  /** return a random integer [0,2^32-1]
   */
  inline uint_type irand()
  {
    return uni.randInt();
  }

  inline generator_type& engine()
  {
    return uni;
  }

  inline int state_size() const
  {
    return uni.SAVE;
  }

  inline void read(std::istream& rin)
  {
    rin >> uni;
  }

  inline void write(std::ostream& rout)
  {
    rout << uni;
  }

  inline void save(std::vector<uint_type>& curstate) const
  {
    curstate.resize(625);
    uni.save(&curstate[0]);
  }

  inline void load(std::vector<uint_type>& newstate)
  {
    uni.load(&newstate[0]);
  }

private:
  ///context number
  int myContext;
  ///number of contexts
  int nContexts;
  ///offset of the random seed
  int baseOffset;
  ///random number generator [0,1)
  generator_type uni;
};
#endif

