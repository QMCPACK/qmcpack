//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF5_ERROR_SUPPRESSION_H
#define QMCPLUSPLUS_HDF5_ERROR_SUPPRESSION_H

#include "hdf5.h"

namespace qmcplusplus
{
/// class suppressing warnings from the HDF5 library
class hdf_error_suppression
{
public:
  /// constructor
  hdf_error_suppression()
  {
    H5Eget_auto2(H5E_DEFAULT, &err_func, &client_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    enabled = true;
  }
  /// We don't want more than one instance at a time.
  hdf_error_suppression(const hdf_error_suppression&) = delete;
  hdf_error_suppression& operator=(const hdf_error_suppression&) = delete;
  /// destructor reset HDF5 error handling
  ~hdf_error_suppression()
  {
    H5Eset_auto2(H5E_DEFAULT, err_func, client_data);
    enabled = false;
  }

  /// status of hdf_error_suppression. An instance of this class changes enabled to true.
  static inline bool enabled{false};

private:
  ///error type
  H5E_auto2_t err_func;
  ///error handling
  void* client_data{nullptr};
};

} // namespace qmcplusplus
#endif //QMCPLUSPLUS_HDF5_ERROR_SUPPRESSION_H

