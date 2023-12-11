////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////


#ifndef AFQMC MATRIX_TRAITS_HPP
#define AFQMC MATRIX_TRAITS_HPP

namespace qmcplusplus
{
namespace afqmc
{
namespace detail
{
// taken from Jean Guegant github/blog
template<typename UnnamedType>
struct container
{
  // Let's put the test in private.
private:
  // We use std::declval to 'recreate' an object of 'UnnamedType'.
  // We use std::declval to also 'recreate' an object of type 'Param'.
  // We can use both of these recreated objects to test the validity!
  template<typename... Params>
  constexpr auto test_validity(int /* unused */)
      -> decltype(std::declval<UnnamedType>()(std::declval<Params>()...), std::true_type())
  {
    // If substitution didn't fail, we can return a true_type.
    return std::true_type();
  }

  template<typename... Params>
  constexpr std::false_type test_validity(...)
  {
    // Our sink-hole returns a false_type.
    return std::false_type();
  }

public:
  // A public operator() that accept the argument we wish to test onto the UnnamedType.
  // Notice that the return type is automatic!
  template<typename... Params>
  constexpr auto operator()(Params&&...)
  {
    // The argument is forwarded to one of the two overloads.
    // The SFINAE on the 'true_type' will come into play to dispatch.
    return test_validity<Params...>(int());
  }
};

template<typename UnnamedType>
constexpr auto is_valid(UnnamedType&& t)
{
  // We used auto for the return type: it will be deduced here.
  return container<UnnamedType>();
}

//  template<typename T, typename
//  is_sparse_matrix


} // namespace detail

} // namespace afqmc

} // namespace qmcplusplus

#endif
