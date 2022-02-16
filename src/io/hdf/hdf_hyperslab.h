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

#include <array>
#include "type_traits/container_traits.h"
#include "hdf_datatype.h"
#include "hdf_dataspace.h"
#include "hdf_dataproxy.h"

namespace qmcplusplus
{
/** class to use file space hyperslab with a serialized container
 * @tparam CT container type, std::vector, Vector, Matrix, Array, boost::multi::array
 * @tparam RANK hyperslab user rank. The dimensions contributed by T are excluded.
 *
 * The container may get resized for sufficient space if
 * the template specialization of container_traits<CT> is available
 * additional restriction:
 * 1D containers can be resized to hold any multi-dimensional data
 * >1D containers can only be resize to hold data with matching dimensions.
 */
template<typename CT, unsigned RANK>
struct hyperslab_proxy
{
  ///user rank of a hyperslab
  static const unsigned int slab_rank = RANK;

  using element_type = typename container_traits<CT>::element_type;
  /** type alias for h5_space_type
   * encapsulates both RANK and ranks contributed by element_type
   * for constructing an hdf5 dataspace.
   */
  using SpaceType = h5_space_type<element_type, RANK>;
  ///global dimension of the hyperslab
  SpaceType file_space;
  ///local dimension of the hyperslab
  SpaceType selected_space;
  ///offset of the hyperslab
  std::array<hsize_t, SpaceType::rank> slab_offset;
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
                         const std::array<IT, RANK>& dims_in,
                         const std::array<IT, RANK>& selected_in,
                         const std::array<IT, RANK>& offsets_in)
      : ref_(a)
  {
    static_assert(std::is_unsigned<IT>::value, "only accept unsigned integer types like size_t");
    for (int i = 0; i < slab_rank; ++i)
    {
      file_space.dims[i]     = static_cast<hsize_t>(dims_in[i]);
      selected_space.dims[i] = static_cast<hsize_t>(selected_in[i]);
      slab_offset[i]         = static_cast<hsize_t>(offsets_in[i]);
    }

    /// element_type related dimensions always have offset 0
    for (int i = slab_rank; i < SpaceType::rank; ++i)
      slab_offset[i] = 0;
  }

  /** return the size of the i-th dimension of global space
   * @param i dimension
   */
  inline hsize_t size(int i) const { return (i > SpaceType::rank) ? 0 : file_space.dims[i]; }

  /** checks if file_space, elected_space and offset are self-consistent
   */
  inline void checkUserRankSizes() const
  {
    if (std::any_of(file_space.dims, file_space.dims + slab_rank, [](int i) { return i == 0; }))
      throw std::runtime_error("Zero size detected in some dimensions of filespace\n");
    if (std::any_of(selected_space.dims, selected_space.dims + slab_rank, [](int i) { return i == 0; }))
      throw std::runtime_error("Zero size detected in some dimensions of selected filespace\n");
    for (int dim = 0; dim < slab_rank; dim++)
    {
      if (slab_offset[dim] > file_space.dims[dim])
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

  /** check if the container is large enough for the selected space and resize if requested
   */
  inline bool checkContainerCapacity() const
  {
    hsize_t total_size = slab_rank > 0 ? 1 : 0;
    for (int dim = 0; dim < slab_rank; dim++)
      total_size *= selected_space.dims[dim];
    return total_size > container_traits<CT>::getSize(ref_) ? false : true;
  }

  /** adjust file_space and selected_space shapes based on sizes_file
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
    }
  }
};

template<typename CT, unsigned RANK>
struct h5data_proxy<hyperslab_proxy<CT, RANK>>
{
  using data_type = hyperslab_proxy<CT, RANK>;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    std::vector<hsize_t> sizes_file;
    getDataShape<typename data_type::element_type>(grp, aname, sizes_file);
    ref.adaptShape(sizes_file);
    ref.checkUserRankSizes();
    if (!ref.checkContainerCapacity())
      container_traits<CT>::resize(ref.ref_, ref.selected_space.dims, ref.slab_rank);
    return h5d_read(grp, aname.c_str(), ref.file_space.rank, ref.file_space.dims, ref.selected_space.dims,
                    ref.slab_offset.data(),
                    hyperslab_proxy<CT, RANK>::SpaceType::get_address(container_traits<CT>::getElementPtr(ref.ref_)),
                    xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    ref.checkUserRankSizes();
    if (!ref.checkContainerCapacity())
      throw std::runtime_error("Not large enough container capacity!\n");
    return h5d_write(grp, aname.c_str(), ref.file_space.rank, ref.file_space.dims, ref.selected_space.dims,
                     ref.slab_offset.data(),
                     hyperslab_proxy<CT, RANK>::SpaceType::get_address(container_traits<CT>::getElementPtr(ref.ref_)),
                     xfer_plist);
  }
};

} // namespace qmcplusplus
#endif
