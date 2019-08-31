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
 * @tparam CT container type, std::vector, Vector, Matrix, Array, boost::multi::array
 * @tparam DIM hyperslab user rank. The dimensions contributed by T are excluded.
 *
 * The container may get resized for sufficient space if
 * the template specialization of container_traits<CT> is available
 * additional restriction:
 * 1D containers can be resized to hold any multi-dimentional data
 * >1D containers can only be resize to hold data with matching dimensions.
 */
template<typename CT, unsigned DIM>
struct hyperslab_proxy
{
  ///user rank of a hyperslab
  static const unsigned int slab_rank = DIM;

  using value_type = typename CT::value_type;
  /** type alias for h5_space_type
   * encapsulates both DIM and ranks contributed by value_type
   * for constructing an hdf5 dataspace.
   */
  using SpaceType = h5_space_type<value_type, DIM>;
  ///global dimension of the hyperslab
  SpaceType file_space;
  ///local dimension of the hyperslab
  SpaceType selected_space;
  ///offset of the hyperslab
  std::array<hsize_t, SpaceType::size()> slab_offset;
  ///container reference
  CT& ref_;

  /** constructor
   * @tparam IT integer type
   * @param a data containter
   * @param dims_in dimension sizes of a dataset.
   * @param selected_in dimension sizes of the selected part of the dataset
   * @param offsets_in offsets of the selected part of the dataset
   *
   * value 0 in any dimension of dims_in or selected_in is allowed when reading a dataset.
   * The actual size is derived from file dataset.
   */
  template<typename IT>
  inline hyperslab_proxy(CT& a,
                         const std::array<IT, DIM>& dims_in,
                         const std::array<IT, DIM>& selected_in,
                         const std::array<IT, DIM>& offsets_in)
      : ref_(a)
  {
    for (int i = 0; i < slab_rank; ++i)
    {
      if (dims_in[i] < 0)
        throw std::runtime_error("Negative size detected in some dimensions of filespace input\n");
      file_space.dims[i] = static_cast<hsize_t>(dims_in[i]);
      if (selected_in[i] < 0)
        throw std::runtime_error("Negative size detected in some dimensions of selected space input\n");
      selected_space.dims[i] = static_cast<hsize_t>(selected_in[i]);
      if (offsets_in[i] < 0)
        throw std::runtime_error("Negative value detected in some dimensions of offset input\n");
      slab_offset[i] = static_cast<hsize_t>(offsets_in[i]);
    }

    /// value_type related dimensions always have offset 0
    for (int i = slab_rank; i < SpaceType::size(); ++i)
      slab_offset[i] = 0;
  }

  /** return the size of the i-th dimension of global space
   * @param i dimension
   */
  inline hsize_t size(int i) const { return (i > SpaceType::size()) ? 0 : file_space.dims[i]; }

  /** checks if file_space, elected_space and offset are self-consistent
   */
  inline void checkUserRankSizes()
  {
    if (std::any_of(file_space.dims, file_space.dims + slab_rank, [](int i) { return i <= 0; }))
      throw std::runtime_error("Non-positive size detected in some dimensions of filespace\n");
    if (std::any_of(selected_space.dims, selected_space.dims + slab_rank, [](int i) { return i <= 0; }))
      throw std::runtime_error("Non-positive size detected in some dimensions of selected filespace\n");
    for (int dim = 0; dim < slab_rank; dim++)
    {
      if (slab_offset[dim] < 0 || slab_offset[dim] > file_space.dims[dim])
        throw std::runtime_error("offset outsdie the bound of filespace");
      if (slab_offset[dim] + selected_space.dims[dim] > file_space.dims[dim])
      {
        std::ostringstream err_msg;
        err_msg << "dim " << dim << " offset " << slab_offset[dim] << " + selected_space size "
                << selected_space.dims[dim] << " outsdie the bound of filespace " << file_space.dims[dim] << std::endl;
        throw std::runtime_error(err_msg.str());
      }
    }
  }

  /** adjust file_space and selected_space shapes and resize the container when needed
   *  @tparam IT integer type
   *  @param sizes_file sizes of all the user dimensions of a dataset
   *
   * This function can only be used when reading a dataset not writing
   * if the size value of a dimension is 0, use the value based on the dataset on disk
   */
  template<typename IT>
  inline void adaptShape(const std::vector<IT>& sizes_file)
  {
    // validate user ranks
    if (sizes_file.size() != slab_rank)
      throw std::runtime_error("User specified and filespace dimensions mismatch!\n");

    hsize_t total_size = 1;
    for (int dim = 0; dim < slab_rank; dim++)
    {
      if (file_space.dims[dim] == 0)
        file_space.dims[dim] = sizes_file[dim];
      else if (file_space.dims[dim] != sizes_file[dim])
      {
        std::ostringstream err_msg;
        err_msg << "dim " << dim << " user specified size " << file_space.dims[dim] << " filespace size "
                << sizes_file[dim] << " mismatched!" << std::endl;
        throw std::runtime_error(err_msg.str());
      }

      if (selected_space.dims[dim] == 0)
        selected_space.dims[dim] = file_space.dims[dim];
      total_size *= selected_space.dims[dim];
    }

    if (total_size != container_traits::getSize(ref_))
      container_traits::resize(ref_, selected_space.dims, slab_rank);
  }
};

template<typename CT, unsigned DIM>
struct h5data_proxy<hyperslab_proxy<CT, DIM>>
{
  hyperslab_proxy<CT, DIM>& ref_;

  h5data_proxy(hyperslab_proxy<CT, DIM>& a) : ref_(a) {}

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    std::vector<hsize_t> sizes_file;
    getDataShape<typename CT::value_type>(grp, aname, sizes_file);
    ref_.adaptShape(sizes_file);
    ref_.checkUserRankSizes();
    return h5d_read(grp, aname.c_str(), ref_.slab_rank, ref_.file_space.dims, ref_.selected_space.dims,
                    ref_.slab_offset.data(), hyperslab_proxy<CT, DIM>::SpaceType::get_address(ref_.ref_.data()),
                    xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    ref_.checkUserRankSizes();
    return h5d_write(grp, aname.c_str(), ref_.slab_rank, ref_.file_space.dims, ref_.selected_space.dims,
                     ref_.slab_offset.data(), hyperslab_proxy<CT, DIM>::SpaceType::get_address(ref_.ref_.data()),
                     xfer_plist);
  }
};

#ifdef BUILD_AFQMC
#ifdef QMC_CUDA
template<typename T, unsigned DIM>
struct h5data_proxy<hyperslab_proxy<boost::multi::array<T, 2, qmc_cuda::cuda_gpu_allocator<T>>, DIM>>
{
  typedef boost::multi::array<T, 2, qmc_cuda::cuda_gpu_allocator<T>> CT;
  hyperslab_proxy<CT, DIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT, DIM>& a) : ref_(a) {}
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

template<typename T, unsigned DIM>
struct h5data_proxy<hyperslab_proxy<boost::multi::array_ref<T, 2, qmc_cuda::cuda_gpu_ptr<T>>, DIM>>
{
  typedef boost::multi::array_ref<T, 2, qmc_cuda::cuda_gpu_ptr<T>> CT;
  hyperslab_proxy<CT, DIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT, DIM>& a) : ref_(a) {}
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
