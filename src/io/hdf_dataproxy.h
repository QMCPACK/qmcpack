//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_H5_DATAPROXY_H
#define QMCPLUSPLUS_HDF_H5_DATAPROXY_H

#include <io/hdf_wrapper_functions.h>
#include <io/hdf_dataspace.h>

namespace qmcplusplus
{

/** generic h5data_proxy<T> for scalar basic datatypes defined in hdf_dataspace.h
 */
template<typename T>
struct h5data_proxy : public h5_space_type<T, 0>
{
  using data_type = T;
  using FileSpace = h5_space_type<T, 0>;
  using FileSpace::dims;
  using FileSpace::get_address;
  data_type& ref_;

  inline h5data_proxy(data_type& a) : ref_(a) { }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_read(grp, aname, get_address(&ref_), xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(&ref_), xfer_plist);
  }

};

/** specialization for bool, convert to int
 */
template<>
struct h5data_proxy<bool> : public h5_space_type<int, 0>
{
  using data_type = bool;
  using FileSpace = h5_space_type<int, 0>;
  using FileSpace::dims;
  using FileSpace::get_address;
  data_type& ref_;

  inline h5data_proxy(data_type& a) : ref_(a) { }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    int copy = static_cast<int>(ref_);
    bool okay = h5d_read(grp, aname, get_address(&copy), xfer_plist);
    ref_ = static_cast<bool>(copy);
    return okay;
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    int copy = static_cast<int>(ref_);
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(&copy), xfer_plist);
  }

};

} // namespace qmcplusplus
#endif
