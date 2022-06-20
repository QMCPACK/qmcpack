///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_STRIDED_RANGE_KERNELS_HPP
#define AFQMC_STRIDED_RANGE_KERNELS_HPP

namespace kernels
{
/* Note: Taken from thrust examples.
 */

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

template<typename Iterator>
class strided_range
{
public:
  using difference_type = typename thrust::iterator_difference<Iterator>::type;

  struct stride_functor : public thrust::unary_function<difference_type, difference_type>
  {
    difference_type stride;

    stride_functor(difference_type stride) : stride(stride) {}

    __host__ __device__ difference_type operator()(const difference_type& i) const { return stride * i; }
  };

  using CountingIterator    = typename thrust::counting_iterator<difference_type>;
  using TransformIterator   = typename thrust::transform_iterator<stride_functor, CountingIterator>;
  using PermutationIterator = typename thrust::permutation_iterator<Iterator, TransformIterator>;

  // type of the strided_range iterator
  using iterator = PermutationIterator;

  // construct strided_range for the range [first,last)
  strided_range(Iterator first, Iterator last, difference_type stride) : first(first), last(last), stride(stride) {}

  iterator begin(void) const
  {
    return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
  }

  iterator end(void) const { return begin() + ((last - first) + (stride - 1)) / stride; }

protected:
  Iterator first;
  Iterator last;
  difference_type stride;
};

} // namespace kernels

#endif
