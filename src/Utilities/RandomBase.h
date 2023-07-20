//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Steven Hahn, hahnse@ornl.gov,
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RANDOMBASE_H
#define QMCPLUSPLUS_RANDOMBASE_H

#include <istream>
#include <memory>
#include <ostream>
#include <vector>

namespace qmcplusplus
{
template<typename T>
class RandomBase
{
public:
  using result_type                                         = T;
  using uint_type                                           = uint_fast32_t;
  virtual ~RandomBase()                                     = default;
  virtual void init(int iseed_in)                           = 0;
  virtual void seed(uint_type aseed)                        = 0;
  virtual T operator()()                                    = 0;
  virtual void write(std::ostream& rout) const              = 0;
  virtual void read(std::istream& rin)                      = 0;
  virtual void load(const std::vector<uint_type>& newstate) = 0;
  virtual void save(std::vector<uint_type>& curstate) const = 0;
  virtual size_t state_size() const                         = 0;
  virtual std::unique_ptr<RandomBase<T>> makeClone() const  = 0;
};

} // namespace qmcplusplus
#endif
