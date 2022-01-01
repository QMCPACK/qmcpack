//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

/** @file
 *  A minimally functional wrapper for the since c++11 <random>
 *
 *  This supports what I saw as the minimal functionality a RandomGenerator type
 *  needed to abstract. Patterned on by radically cut down from BoostRandom.h
 *
 *  only used by test code
 */
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
class StdRandom
{
public:
  using result_type = T;
  using Engine = std::mt19937;
  using uint_type   = Engine::result_type;

  StdRandom(uint_type iseed = 911) : engine(iseed) { }

  void init(int iseed_in)
  {
    uint_type baseSeed = iseed_in;
    engine.seed(baseSeed);
  }

  void seed(uint_type aseed) { engine.seed(aseed); }

  result_type rand() { return distribution(engine); }
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
  std::string ClassName{"StdRand"};
  std::string EngineName{"std::mt19937"};

private:
  ///random number generator [0,1)
  static constexpr double min = 0.0;
  static constexpr double max = 1.0;
  std::uniform_real_distribution<T> distribution{min, max};
  Engine engine;
};

} // namespace qmcplusplus

#endif
