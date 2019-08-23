/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_HYPERSLAB_IO_H
#define QMCPLUSPLUS_HDF_HYPERSLAB_IO_H
#include <type_traits/container_traits.h>
#include <io/hdf_datatype.h>
#include <io/hdf_dataspace.h>
#include <io/hdf_dataproxy.h>

#ifdef BUILD_AFQMC
#ifdef QMC_CUDA
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#endif
#endif

namespace qmcplusplus
{
/** class to use file space hyperslab with a serialized container
 *
 * container_proxy<CT> handles the size and datatype
 */
template<typename CT, unsigned DIM>
struct hyperslab_proxy
{
  ///rank of hyperslab
  static const unsigned int slab_rank = DIM;
  using SpaceType = h5_space_type<typename CT::value_type, DIM>;
  ///global dimension of the hyperslab
  SpaceType global_space;
  ///local dimension of the hyperslab
  SpaceType local_space;
  ///offset of the hyperslab
  std::array<hsize_t, SpaceType::size()> slab_offset;
  ///container reference
  CT& ref_;

  template<typename IC>
  inline hyperslab_proxy(CT& a, const IC& dims_in) : ref_(a)
  {
/*
    use_slab  = false;
    slab_rank = dims_in.size();
    for (int i = 0; i < dims_in.size(); ++i)
      slab_dims[i] = static_cast<hsize_t>(dims_in[i]);
    slab_dims_local = slab_dims;
    if (element_size > 1)
    {
      slab_dims[slab_rank] = element_size;
      slab_rank += 1;
    }
*/
  }

  inline hyperslab_proxy(CT& a, const std::array<int, DIM>& dims_in, const std::array<int, DIM>& dims_loc, const std::array<int, DIM>& offsets_in) : ref_(a)
  {
    for (int i = 0; i < slab_rank; ++i)
      global_space.dims[i] = static_cast<hsize_t>(dims_in[i]);
    for (int i = 0; i < slab_rank; ++i)
      local_space.dims[i] = static_cast<hsize_t>(dims_loc[i]);
    for (int i = 0; i < slab_rank; ++i)
      slab_offset[i] = static_cast<hsize_t>(offsets_in[i]);
    /// value_type related dimensions always have offset 0
    for (int i = slab_rank; i < SpaceType::size(); ++i)
      slab_offset[i] = 0;
  }

  /** return the size of the i-th dimension of global space
   * @param i dimension
   */
  inline hsize_t size(int i) const { return (i > SpaceType::size()) ? 0 : global_space.dims[i]; }

  inline void change_shape() { 
//this->resize(slab_dims.data(), MAXDIM + 1);
 }
};

template<typename CT, unsigned MAXDIM>
struct h5data_proxy<hyperslab_proxy<CT, MAXDIM>>
{
  hyperslab_proxy<CT, MAXDIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT, MAXDIM>& a) : ref_(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_read(grp, aname.c_str(), ref_.slab_rank, ref_.global_space.dims, ref_.local_space.dims,
                    ref_.slab_offset.data(), hyperslab_proxy<CT, MAXDIM>::SpaceType::get_address(ref_.ref_.data()), xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_write(grp, aname.c_str(), ref_.slab_rank, ref_.global_space.dims, ref_.local_space.dims,
                     ref_.slab_offset.data(), hyperslab_proxy<CT, MAXDIM>::SpaceType::get_address(ref_.ref_.data()), xfer_plist);
  }
};

#ifdef BUILD_AFQMC
#ifdef QMC_CUDA
template<typename T, unsigned MAXDIM>
struct h5data_proxy<hyperslab_proxy<boost::multi::array<T, 2, qmc_cuda::cuda_gpu_allocator<T>>, MAXDIM>>
{
  typedef boost::multi::array<T, 2, qmc_cuda::cuda_gpu_allocator<T>> CT;
  hyperslab_proxy<CT, MAXDIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT, MAXDIM>& a) : ref_(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref_.use_slab)
    {
      // later on specialize h5d_read for fancy pointers
      std::size_t sz = ref_.ref.num_elements();
      boost::multi::array<T, 1> buf(typename boost::multi::layout_t<1u>::extensions_type{sz});
      auto ret = h5d_read(grp, aname.c_str(), ref_.slab_rank, ref_.slab_dims.data(), ref_.slab_dims_local.data(),
                          ref_.slab_offset.data(), buf.origin(), xfer_plist);
      qmc_cuda::copy_n(buf.data(), sz, ref_.ref.origin());
      return ret;
    }
    else
    {
      int rank = ref_.slab_rank;
      if (!get_space(grp, aname, rank, ref_.slab_dims.data(), true))
      {
        std::cerr << " Disabled hyperslab resize with boost::multi::array<gpu_allocator>.\n";
        return false;
      }
      return h5d_read(grp, aname, ref_.data(), xfer_plist);
    }
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    std::cerr << " Disabled hyperslab write with boost::multi::array<gpu_allocator>.\n";
    return false;
  }
};

template<typename T, unsigned MAXDIM>
struct h5data_proxy<hyperslab_proxy<boost::multi::array_ref<T, 2, qmc_cuda::cuda_gpu_ptr<T>>, MAXDIM>>
{
  typedef boost::multi::array_ref<T, 2, qmc_cuda::cuda_gpu_ptr<T>> CT;
  hyperslab_proxy<CT, MAXDIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT, MAXDIM>& a) : ref_(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref_.use_slab)
    {
      // later on specialize h5d_read for fancy pointers
      std::size_t sz = ref_.ref.num_elements();
      boost::multi::array<T, 1> buf(typename boost::multi::layout_t<1u>::extensions_type{sz});
      auto ret = h5d_read(grp, aname.c_str(), ref_.slab_rank, ref_.slab_dims.data(), ref_.slab_dims_local.data(),
                          ref_.slab_offset.data(), buf.origin(), xfer_plist);
      qmc_cuda::copy_n(buf.data(), sz, ref_.ref.origin());
      return ret;
    }
    else
    {
      int rank = ref_.slab_rank;
      if (!get_space(grp, aname, rank, ref_.slab_dims.data(), true))
      {
        std::cerr << " Disabled hyperslab resize with boost::multi::array_ref<gpu_ptr>.\n";
        return false;
      }
      return h5d_read(grp, aname, ref_.data(), xfer_plist);
    }
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    std::cerr << " Disabled hyperslab write with boost::multi::array_ref<gpu_ptr>.\n";
    return false;
  }
};
#endif
#endif
} // namespace qmcplusplus
#endif
