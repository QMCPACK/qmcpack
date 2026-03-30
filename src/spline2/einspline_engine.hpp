//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EINSPLINE_ENGINE_HPP
#define QMCPLUSPLUS_EINSPLINE_ENGINE_HPP

#include "bspline_traits.hpp"
#include "hdf/hdf_archive.h"

namespace qmcplusplus
{
/** einspline_engine
   *
   * for handling hdf5 I/O
   */
template<typename T, unsigned D>
class einspline_engine
{
public:
  using real_type  = typename bspline_traits<T, D>::real_type;
  using value_type = typename bspline_traits<T, D>::value_type;
  using SplineType = typename bspline_traits<T, D>::SplineType;
  ///spline engine
  SplineType& spliner;

  einspline_engine(SplineType& s) : spliner(s) {}
};

template<unsigned DIM>
struct dim_traits
{};

// for 3D multi
template<>
struct dim_traits<4>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner.x_grid.num + 3;
    dims[1] = a.spliner.y_grid.num + 3;
    dims[2] = a.spliner.z_grid.num + 3;
    dims[3] = a.spliner.z_stride;
  }
};

// for 1D multi
template<>
struct dim_traits<2>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner.x_grid.num + 2;
    dims[1] = a.spliner.x_stride;
  }
};

/** specialization of h5data_proxy for einspline_engine
   */
template<typename T, unsigned D>
struct h5data_proxy<einspline_engine<T, D>> : public h5_space_type<T, D + 1>
{
  using value_type = T;
  using Base       = h5_space_type<T, D + 1>;
  using Base::dims;
  using Base::get_address;
  using data_type = einspline_engine<T, D>;

  inline h5data_proxy(const data_type& a) { dim_traits<D + 1>::setdim(a, dims); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  { return h5d_read(grp, aname, get_address(ref.spliner.coefs), xfer_plist); }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  { return h5d_write(grp, aname.c_str(), Base::rank, dims, get_address(ref.spliner.coefs), xfer_plist); }
};

} // namespace qmcplusplus

#endif
