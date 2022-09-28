//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PREFETCHED_RANGE_H
#define QMCPLUSPLUS_PREFETCHED_RANGE_H

namespace qmcplusplus
{
/// helper class for the prefetched range of a vector
class PrefetchedRange
{
  // [first, last) rows of Ainv
  int first, last;

public:
  PrefetchedRange() : first(0), last(0){};
  void setRange(int first_in, int last_in)
  {
    first = first_in;
    last  = last_in;
  }
  inline void clear() { first = last = 0; };
  inline int getOffset(int index) const
  {
    if (!checkRange(index))
      throw std::runtime_error("index not in range \n");
    return index - first;
  }
  inline bool checkRange(int index) const { return (index >= first) && (index < last); };
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_PREFETCHED_RANGE_H
