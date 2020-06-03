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


#ifndef QMCPLUSPLUS_HDF_PETE_TRAITS_H
#define QMCPLUSPLUS_HDF_PETE_TRAITS_H

#include <io/hdf_dataproxy.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>

namespace qmcplusplus
{
/** specialization for Vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T>
struct h5data_proxy<Vector<T>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  typedef Vector<T> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a) : ref_(a) { dims[0] = ref_.size(); }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref_.resize(dims[0]);
    return h5d_read(grp, aname, get_address(ref_.data()), xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(ref_.data()), xfer_plist);
  }
};


template<typename T>
struct h5data_proxy<Matrix<T>> : public h5_space_type<T, 2>
{
  using FileSpace = h5_space_type<T, 2>;
  using FileSpace::dims;
  using FileSpace::get_address;
  typedef Matrix<T> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a) : ref_(a)
  {
    dims[0] = ref_.rows();
    dims[1] = ref_.cols();
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref_.resize(dims[0], dims[1]);
    return h5d_read(grp, aname, get_address(ref_.data()), xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(ref_.data()), xfer_plist);
  }
};


template<typename T, unsigned D>
struct h5data_proxy<Array<T, D>> : public h5_space_type<T, D>
{
  using FileSpace = h5_space_type<T, D>;
  using FileSpace::dims;
  using FileSpace::get_address;
  typedef Array<T, D> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a) : ref_(a)
  {
    for (int i = 0; i < D; ++i)
      dims[i] = ref_.size(i);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref_.resize(dims);
    return h5d_read(grp, aname, get_address(ref_.data()), xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(ref_.data()), xfer_plist);
  }
};

} // namespace qmcplusplus
#endif
