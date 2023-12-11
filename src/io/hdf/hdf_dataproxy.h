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

#include "hdf_wrapper_functions.h"
#include "hdf_dataspace.h"

namespace qmcplusplus
{

/** generic h5data_proxy<T> for scalar basic datatypes defined in hdf_dataspace.h
 * Note if the dataset to be written has const specifier, T should not carry const.
 */
template<typename T>
struct h5data_proxy : public h5_space_type<T, 0>
{
  using data_type = T;
  using FileSpace = h5_space_type<T, 0>;
  using FileSpace::dims;
  using FileSpace::get_address;

  inline h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_read(grp, aname, get_address(&ref), xfer_plist);
  }

  inline bool check_existence(hid_t grp, const std::string& aname) { return h5d_check_existence(grp, aname); }

  inline bool check_type(hid_t grp, const std::string& aname) { return h5d_check_type<data_type>(grp, aname); }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(&ref), xfer_plist);
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

  inline h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    int copy  = static_cast<int>(ref);
    bool okay = h5d_read(grp, aname, get_address(&copy), xfer_plist);
    ref       = static_cast<bool>(copy);
    return okay;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    int copy = static_cast<int>(ref);
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(&copy), xfer_plist);
  }
};

} // namespace qmcplusplus
#endif
