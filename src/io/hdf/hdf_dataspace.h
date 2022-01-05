//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H
#define QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H

/**@file hdf_dataspace.h
 * @brief define h5_space_type to handle basic datatype for hdf5
 *
 * h5_space_type is a helper class for h5data_proxy and used internally
 * - h5_space_type<T,RANK>
 * - h5_space_type<std::complex<T>,RANK>
 * - h5_space_type<std::array<T,D>,RANK>
 * - h5_space_type<TinyVector<T,D>,RANK>
 * - h5_space_type<TinyVector<std::complex<T>,D>,RANK> // removed, picked up by template recursion
 * - h5_space_type<Tensor<T,D>,RANK>
 * - h5_space_type<Tensor<std::complex<T>,D>,RANK> // removed, picked up by template recursion
 */

#include <complex>
#include <array>
#include "hdf_datatype.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{
/** default struct to define a h5 dataspace, any intrinsic type T
 *
 * @tparam T intrinsic datatype
 * @tparam RANK rank of the multidimensional h5dataspace
 */
template<typename T, hsize_t RANK>
struct h5_space_type
{
  ///shape of the dataspace, protected for zero size array, hdf5 support scalar as rank = 0
  hsize_t dims[RANK > 0 ? RANK : 1];
  ///rank of the multidimensional dataspace
  static constexpr hsize_t rank = RANK;
  ///new rank added due to T
  static constexpr int added_rank() { return 0; }
  ///return the address
  inline static auto get_address(T* a) { return a; }
  inline static auto get_address(const T* a) { return a; }
};

/** specialization of h5_space_type for std::complex<T>
 *
 * Raize the dimension of the space by 1 and set the last dimension=2
 */
template<typename T, hsize_t RANK>
struct h5_space_type<std::complex<T>, RANK> : public h5_space_type<T, RANK + 1>
{
  using Base = h5_space_type<T, RANK + 1>;
  using Base::dims;
  using Base::rank;
  static constexpr int added_rank() { return Base::added_rank() + 1; }
  inline h5_space_type() { dims[RANK] = 2; }
  inline static auto get_address(std::complex<T>* a) { return Base::get_address(reinterpret_cast<T*>(a)); }
  inline static auto get_address(const std::complex<T>* a) { return Base::get_address(reinterpret_cast<const T*>(a)); }
};

/** specialization of h5_space_type for std::array<T,D> for any intrinsic type T
 */
template<typename T, std::size_t D, hsize_t RANK>
struct h5_space_type<std::array<T, D>, RANK> : public h5_space_type<T, RANK + 1>
{
  using Base = h5_space_type<T, RANK + 1>;
  using Base::dims;
  using Base::rank;
  inline h5_space_type() { dims[RANK] = D; }
  static constexpr int added_rank() { return Base::added_rank() + 1; }
  inline static auto get_address(std::array<T, D>* a) { return Base::get_address(a->data()); }
  inline static auto get_address(const std::array<T, D>* a) { return Base::get_address(a->data()); }
};

/** specialization of h5_space_type for TinyVector<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t RANK>
struct h5_space_type<TinyVector<T, D>, RANK> : public h5_space_type<T, RANK + 1>
{
  using Base = h5_space_type<T, RANK + 1>;
  using Base::dims;
  using Base::rank;
  inline h5_space_type() { dims[RANK] = D; }
  static constexpr int added_rank() { return Base::added_rank() + 1; }
  inline static auto get_address(TinyVector<T, D>* a) { return Base::get_address(a->data()); }
  inline static auto get_address(const TinyVector<T, D>* a) { return Base::get_address(a->data()); }
};

/** specialization of h5_space_type for Tensor<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t RANK>
struct h5_space_type<Tensor<T, D>, RANK> : public h5_space_type<T, RANK + 2>
{
  using Base = h5_space_type<T, RANK + 2>;
  using Base::dims;
  using Base::rank;
  inline h5_space_type()
  {
    dims[RANK]     = D;
    dims[RANK + 1] = D;
  }
  static constexpr int added_rank() { return Base::added_rank() + 2; }
  inline static auto get_address(Tensor<T, D>* a) { return Base::get_address(a->data()); }
  inline static auto get_address(const Tensor<T, D>* a) { return Base::get_address(a->data()); }
};

} // namespace qmcplusplus
#endif
