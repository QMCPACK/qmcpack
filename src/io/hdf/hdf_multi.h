//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
////
//// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////
//// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_HDF_MULTI_INTERFACE_H
#define QMCPLUSPLUS_HDF_MULTI_INTERFACE_H

#include <multi/array.hpp>
#include <multi/array_ref.hpp>
#include "hdf_dataproxy.h"
#include "hdf_hyperslab.h"

#ifdef BUILD_AFQMC
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/device_pointers.hpp"
#endif
#endif

namespace qmcplusplus
{
/** specialization for vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T, class Alloc>
struct h5data_proxy<boost::multi::array<T, 1, Alloc>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array<T, 1, Alloc>;

  inline h5data_proxy(const data_type& a) { dims[0] = a.num_elements(); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    using iextensions = typename boost::multi::iextensions<1u>;
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref.reextent(iextensions{static_cast<boost::multi::size_t>(dims[0])});
    return h5d_read(grp, aname, get_address(std::addressof(*ref.origin())), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(std::addressof(*ref.origin())), xfer_plist);
  }
};

template<typename T, class Alloc>
struct h5data_proxy<boost::multi::array<T, 2, Alloc>> : public h5_space_type<T, 2>
{
  using FileSpace = h5_space_type<T, 2>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array<T, 2, Alloc>;

  inline h5data_proxy(const data_type& a)
  {
    dims[0] = a.size(0);
    dims[1] = a.size(1);
  }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref.reextent({static_cast<boost::multi::size_t>(dims[0]), static_cast<boost::multi::size_t>(dims[1])});
    return h5d_read(grp, aname, get_address(std::addressof(*ref.origin())), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(std::addressof(*ref.origin())), xfer_plist);
  }
};

template<typename T, class Ptr>
struct h5data_proxy<boost::multi::array_ref<T, 1, Ptr>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array_ref<T, 1, Ptr>;

  inline h5data_proxy(const data_type& a) { dims[0] = a.num_elements(); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
    {
      if (dims[0] > 0)
      {
        std::cerr << " Error: multi::array_ref can't be resized in h5data_proxy<>::read." << std::endl;
        std::cerr << dims[0] << " " << ref.size(0) << std::endl;
      }
      return false;
    }
    return h5d_read(grp, aname, get_address(std::addressof(*ref.origin())), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(std::addressof(*ref.origin())), xfer_plist);
  }
};

template<typename T, class Ptr>
struct h5data_proxy<boost::multi::array_ref<T, 2, Ptr>> : public h5_space_type<T, 2>
{
  using FileSpace = h5_space_type<T, 2>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array_ref<T, 2, Ptr>;

  inline h5data_proxy(const data_type& a)
  {
    dims[0] = a.size(0);
    dims[1] = a.size(1);
  }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
    {
      if (dims[0] * dims[1] > 0)
      {
        std::cerr << " Error: multi::array_ref can't be resized in h5data_proxy<>::read." << std::endl;
        std::cerr << dims[0] << " " << dims[1] << " " << ref.size(0) << " " << ref.size(1) << std::endl;
      }
      return false;
    }
    return h5d_read(grp, aname, get_address(std::addressof(*ref.origin())), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(std::addressof(*ref.origin())), xfer_plist);
  }
};


#ifdef BUILD_AFQMC
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
// Specializations for cuda_gpu_allocator
// Need buffered I/O and copies to gpu
template<typename T>
struct h5data_proxy<boost::multi::array<T, 1, device::device_allocator<T>>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array<T, 1, device::device_allocator<T>>;

  inline h5data_proxy(const data_type& a) { dims[0] = a.num_elements(); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref.reextent({dims[0]});
    auto sz           = ref.num_elements();
    using iextensions = typename boost::multi::iextensions<1u>;
    boost::multi::array<T, 1> buf(iextensions{sz});
    auto ret = h5d_read(grp, aname, get_address(buf.data()), xfer_plist);
    device::copy_n(buf.data(), sz, ref.origin());
    return ret;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(to_address(ref.origin())), xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array<T, 2, device::device_allocator<T>>> : public h5_space_type<T, 2>
{
  using FileSpace = h5_space_type<T, 2>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array<T, 2, device::device_allocator<T>>;

  inline h5data_proxy(const data_type& a)
  {
    dims[0] = a.size(0);
    dims[1] = a.size(1);
  }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref.reextent({dims[0], dims[1]});
    auto sz           = ref.num_elements();
    using iextensions = typename boost::multi::iextensions<1u>;
    boost::multi::array<T, 1> buf(iextensions{sz});
    auto ret = h5d_read(grp, aname, get_address(buf.data()), xfer_plist);
    device::copy_n(buf.data(), sz, ref.origin());
    return ret;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(to_address(ref.origin())), xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array_ref<T, 1, device::device_pointer<T>>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array_ref<T, 1, device::device_pointer<T>>;

  inline h5data_proxy(const data_type& a) { dims[0] = a.num_elements(); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
    {
      if (dims[0] > 0)
      {
        std::cerr << " Error: multi::array_ref can't be resized in h5data_proxy<>::read." << std::endl;
        std::cerr << dims[0] << " " << ref.size(0) << std::endl;
      }
      return false;
    }
    auto sz           = ref.num_elements();
    using iextensions = typename boost::multi::iextensions<1u>;
    boost::multi::array<T, 1> buf(iextensions{sz});
    auto ret = h5d_read(grp, aname, get_address(buf.data()), xfer_plist);
    device::copy_n(buf.data(), sz, ref.origin());
    return ret;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(to_address(ref.origin())), xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array_ref<T, 2, device::device_pointer<T>>> : public h5_space_type<T, 2>
{
  using FileSpace = h5_space_type<T, 2>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = boost::multi::array_ref<T, 2, device::device_pointer<T>>;

  inline h5data_proxy(const data_type& a)
  {
    dims[0] = a.size(0);
    dims[1] = a.size(1);
  }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
    {
      if (dims[0] * dims[1] > 0)
      {
        std::cerr << " Error: multi::array_ref can't be resized in h5data_proxy<>::read." << std::endl;
        std::cerr << dims[0] << " " << dims[1] << " " << ref.size(0) << " " << ref.size(1) << std::endl;
      }
      return false;
    }
    auto sz           = ref.num_elements();
    using iextensions = typename boost::multi::iextensions<1u>;
    boost::multi::array<T, 1> buf(iextensions{sz});
    auto ret = h5d_read(grp, aname, get_address(buf.data()), xfer_plist);
    device::copy_n(buf.data(), sz, ref.origin());
    return ret;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(to_address(ref.origin())), xfer_plist);
  }
};

template<typename T, unsigned RANK>
struct h5data_proxy<hyperslab_proxy<boost::multi::array<T, 2, device::device_allocator<T>>, RANK>>
{
  using CT        = boost::multi::array<T, 2, device::device_allocator<T>>;
  using data_type = hyperslab_proxy<CT, RANK>;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref.use_slab)
    {
      // later on specialize h5d_read for fancy pointers
      auto sz = ref.ref.num_elements();
      boost::multi::array<T, 1> buf(typename boost::multi::layout_t<1u>::extensions_type{sz});
      auto ret = h5d_read(grp, aname.c_str(), ref.slab_rank, ref.slab_dims.data(), ref.slab_dims_local.data(),
                          ref.slab_offset.data(), buf.origin(), xfer_plist);
      device::copy_n(buf.data(), sz, ref.ref.origin());
      return ret;
    }
    else
    {
      int rank = ref.slab_rank;
      if (!checkShapeConsistency<T>(grp, aname, rank, ref.slab_dims.data(), true))
      {
        std::cerr << " Disabled hyperslab resize with boost::multi::array<gpu_allocator>.\n";
        return false;
      }
      return h5d_read(grp, aname, ref.data(), xfer_plist);
    }
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    std::cerr << " Disabled hyperslab write with boost::multi::array<gpu_allocator>.\n";
    return false;
  }
};

template<typename T, unsigned RANK>
struct h5data_proxy<hyperslab_proxy<boost::multi::array_ref<T, 2, device::device_pointer<T>>, RANK>>
{
  using CT        = boost::multi::array_ref<T, 2, device::device_pointer<T>>;
  using data_type = hyperslab_proxy<CT, RANK>;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref.use_slab)
    {
      // later on specialize h5d_read for fancy pointers
      auto sz = ref.ref.num_elements();
      boost::multi::array<T, 1> buf(typename boost::multi::layout_t<1u>::extensions_type{sz});
      auto ret = h5d_read(grp, aname.c_str(), ref.slab_rank, ref.slab_dims.data(), ref.slab_dims_local.data(),
                          ref.slab_offset.data(), buf.origin(), xfer_plist);
      device::copy_n(buf.data(), sz, ref.ref.origin());
      return ret;
    }
    else
    {
      int rank = ref.slab_rank;
      if (!checkShapeConsistency<T>(grp, aname, rank, ref.slab_dims.data(), true))
      {
        std::cerr << " Disabled hyperslab resize with boost::multi::array_ref<gpu_ptr>.\n";
        return false;
      }
      return h5d_read(grp, aname, ref.data(), xfer_plist);
    }
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    std::cerr << " Disabled hyperslab write with boost::multi::array_ref<gpu_ptr>.\n";
    return false;
  }
};

#endif
#endif

} // namespace qmcplusplus
#endif
