//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H
#define QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H
/**@file hdf_dataspace.h
 * @brief define h5_space_type to handle basic datatype for hdf5
 *
 * h5_space_type is a helper class for h5data_proxy and used internally
 * - h5_space_type<T,DS>
 * - h5_space_type<std::complex<T>,DS>
 * - h5_space_type<TinyVector<T,D>,DS>
 * - h5_space_type<Tensor<T,D>,DS>
 */


#include <io/hdf_datatype.h>
#include <complex>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>

namespace qmcplusplus
{
/** default struct to define a h5 dataspace, any intrinsic type T
 *
 * \tparm T intrinsic datatype
 * \tparm D dimension of the h5dataspace
 */
template<typename T, hsize_t D, typename ENABLE = void>
struct h5_space_type
{
  ///shape of the dataspace
  hsize_t dims[D];
  ///dimension of the dataspace
  constexpr int size() const { return D; }
  ///new dimension added due to T
  constexpr int added_size() const { return 0; }
  ///return the address
  inline static auto get_address(T* a) { return a; }
};

/** specialization of h5_space_type for scalar intrinsic type T
 */
template<typename T>
struct h5_space_type<T, 0, std::enable_if_t<std::is_floating_point<T>::value || std::is_integral<T>::value>>
{
  hsize_t dims[1];
  inline h5_space_type() { dims[0] = 1; }
  constexpr int size() const { return 1; }
  constexpr int added_size() const { return 1; }
  inline static auto get_address(T* a) { return a; }
};

/** specialization of h5_space_type for std::complex<T>
 *
 * Raize the dimension of the space by 1 and set the last dimension=2
 */
template<typename T, hsize_t D>
struct h5_space_type<std::complex<T>, D> : public h5_space_type<T, D + 1>
{
  using base = h5_space_type<T, D + 1>;
  using base::dims;
  using base::size;
  constexpr int added_size() const { return base::added_size() + 1; }
  inline h5_space_type() { dims[D] = 2; }
  inline static auto get_address(std::complex<T>* a) { return base::get_address(reinterpret_cast<T*>(a)); }
};

/** specialization of h5_space_type for TinyVector<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<TinyVector<T, D>, DS> : public h5_space_type<T, DS + 1>
{
  using base = h5_space_type<T, DS + 1>;
  using base::dims;
  using base::size;
  inline h5_space_type() { dims[DS] = D; }
  constexpr int added_size() const { return base::added_size() + 1; }
  inline static auto get_address(TinyVector<T, D>* a) { return base::get_address(a->data()); }
};

/** specialization of h5_space_type for TinyVector<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<Tensor<T, D>, DS> : public h5_space_type<T, DS + 2>
{
  using base = h5_space_type<T, DS + 2>;
  using base::dims;
  using base::size;
  inline h5_space_type()
  {
    dims[DS]     = D;
    dims[DS + 1] = D;
  }
  constexpr int added_size() const { return base::added_size() + 2; }
  inline static auto get_address(Tensor<T, D>* a) { return base::get_address(a->data()); }
};

} // namespace qmcplusplus
#endif
