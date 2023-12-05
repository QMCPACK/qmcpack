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

#ifndef AFQMC_STRIDED_2DRANGE_KERNELS_HPP
#define AFQMC_STRIDED_2DRANGE_KERNELS_HPP

namespace kernels
{
/* Note: Taken from thrust examples.
 */

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

template<typename Iterator>
class strided_2Drange
{
public:
  using difference_type = typename thrust::iterator_difference<Iterator>::type;

  struct stride_functor : public thrust::unary_function<difference_type, difference_type>
  {
    difference_type stride;
    difference_type nrows;

    stride_functor(difference_type stride, difference_type nr) : stride(stride), nrows(nr) {}

    __host__ __device__ difference_type operator()(const difference_type& n) const
    {
      int j = n / nrows;
      int i = n % nrows;
      return j * stride + i;
    }
  };

  using CountingIterator    = typename thrust::counting_iterator<difference_type>;
  using TransformIterator   = typename thrust::transform_iterator<stride_functor, CountingIterator>;
  using PermutationIterator = typename thrust::permutation_iterator<Iterator, TransformIterator>;

  // type of the strided_2Drange iterator
  using iterator = PermutationIterator;

  // construct strided_2Drange for the 2Drange [first,last)
  strided_2Drange(Iterator first, Iterator last, difference_type stride, difference_type nr)
      : first(first), last(last), stride(stride), nrows(nr)
  {}

  iterator begin(void) const
  {
    return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride, nrows)));
  }

  iterator end(void) const
  {
    int j = (last - first) / stride;
    int i = (last - first) % stride;
    return begin() + j * nrows + i;
  }

protected:
  Iterator first;
  Iterator last;
  difference_type stride;
  difference_type nrows;
};

} // namespace kernels

#endif
