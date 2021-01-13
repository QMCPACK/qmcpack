//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_STDRAND_H
#define QMCPLUSPLUS_STDRAND_H

#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <iterator>

namespace qmcplusplus
{
template<typename T>
class StdRand
{
public:
  using result_type = T;
  using uint_type   = std::uint_fast64_t;
  void init(int i, int nstr, int iseed_in, uint_type offset = 1)
  {
    uint_type baseSeed = iseed_in;
    myContext          = i;
    nContexts          = nstr;
    // if (iseed_in <= 0)
    //   baseSeed = makeSeed(i, nstr);
    engine.seed(baseSeed);
  }

  int offset() const { return baseOffset; }
  void seed(uint_type aseed) { engine.seed(aseed); }

  result_type operator()() { return distribution(engine); }
  void write(std::ostream& rout) const { rout << engine; }
  void read(std::istream& rin) { rin >> engine; }
  size_t state_size() { return engine.state_size; }

  void load(const std::vector<uint_type>& newstate)
  {
    std::stringstream otemp;
    copy(newstate.begin(), newstate.end(), std::ostream_iterator<uint_type>(otemp, " "));
    otemp >> engine;
  }

  void save(std::vector<uint_type>& curstate) const
  {
    curstate.clear();
    std::stringstream otemp;
    otemp << engine;
    copy(std::istream_iterator<uint_type>(otemp), std::istream_iterator<uint_type>(), back_inserter(curstate));
  }

public:
  // Non const allows use of default copy constructor
  // We copy construct from children a boggling amount, is the aim to repeat sequences?
  // But then why does MomentumEstimator make a copy too...
  // Shouldn't the chain of custody for the prng's be practically
  // the most important thing in a QMC code?
  std::string ClassName{"StdRand"};
  std::string EngineName{"std::mersenne_twister_engine"};
  
private:
  static constexpr double min = 0.0;
  static constexpr double max = 1.0;
  std::uniform_real_distribution<T> distribution{min, max};
  std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
                             0x9908b0df, 11,
                             0xffffffff, 7,
                             0x9d2c5680, 15,
                             0xefc60000, 18, 1812433253>
      engine;

  ///context number
  int myContext;
  ///number of contexts
  int nContexts;
  ///offset of the random seed
  int baseOffset;
  ///random number generator [0,1)
};

} // namespace qmcplusplus

#endif
