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
#include <sstream>
#include <limits>
#include <boost/config.hpp>
#ifdef BOOST_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
#endif
#include <boost/random.hpp>
#ifdef BOOST_CLANG
#pragma clang diagnostic pop
#endif

/** random number generator using boost::random
 *
 * A wrapper of boost::random class to work with applicatoins.
 */
template<typename T, typename RNG = boost::mt19937>
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
  typedef boost::variate_generator<generator_type, boost::uniform_real<T>> uniform_generator_type;

  std::string ClassName;
  std::string EngineName;

  ///default constructor
  explicit BoostRandom(uint_type iseed = 911, const std::string& aname = "mt19937")
      : ClassName("boost"),
        EngineName(aname),
        uni(generator_type(iseed), boost::uniform_real<T>(0, 1))
  {}

  /** initialize the generator
   * @param iseed_in input seed
   *
   * Initialize generator with the seed.
   */
  void init(int iseed_in);

  ///assign seed
  inline void seed(uint_type aseed) { uni.engine().seed(aseed); }

  /** return a random number [0,1)
   */
  inline result_type rand() { return uni(); }

  /** return a random number [0,1)
   */
  inline result_type operator()() { return uni(); }

  inline size_t state_size() const { return uni.engine().state_size; }

  inline void read(std::istream& rin) { rin >> uni.engine(); }

  inline void write(std::ostream& rout) const { rout << uni.engine(); }

  inline void save(std::vector<uint_type>& curstate) const
  {
    curstate.clear();
    std::stringstream otemp;
    otemp << uni.engine();
    copy(std::istream_iterator<uint_type>(otemp), std::istream_iterator<uint_type>(), back_inserter(curstate));
  }

  inline void load(const std::vector<uint_type>& newstate)
  {
    std::stringstream otemp;
    copy(newstate.begin(), newstate.end(), std::ostream_iterator<uint_type>(otemp, " "));
    otemp >> uni.engine();
  }

private:
  ///random number generator [0,1)
  uniform_generator_type uni;
};

extern template class BoostRandom<float>;
extern template class BoostRandom<double>;
#endif
