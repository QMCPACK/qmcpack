/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_DOUBLE_HYPERSLAB_IO_H
#define QMCPLUSPLUS_HDF_DOUBLE_HYPERSLAB_IO_H
#include "type_traits/container_proxy.h"
#include "hdf_datatype.h"
#include "hdf_dataspace.h"
#include "hdf_dataproxy.h"
namespace qmcplusplus
{
/** class to use hyperslabs in both file and memory spaces
 *
 * container_proxy<CT> handles the size and datatype
 */
template<typename CT, unsigned MAXDIM>
struct double_hyperslab_proxy : public container_proxy<CT>
{
  ///determine the size of value_type
  enum
  {
    element_size = container_proxy<CT>::DIM
  };
  ///rank of hyperslab
  int slab_rank;
  int mem_rank;
  ///true, if hyperslab is used
  bool use_slab;
  ///global dimension of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> slab_dims;
  ///local dimension of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> slab_dims_local;
  ///offset of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> slab_offset;
  ///global dimension of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> mem_dims;
  ///local dimension of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> mem_dims_local;
  ///offset of the hyperslab
  TinyVector<hsize_t, MAXDIM + 1> mem_offset;
  ///1D
  double_hyperslab_proxy(CT& a)
      : container_proxy<CT>(a), slab_rank(a.slab_rank), slab_dims(a.slab_dims), slab_offset(a.slab_offset)
  {
    slab_dims_local = slab_dims;
    use_slab        = false;
  }

  template<typename IC>
  inline double_hyperslab_proxy(CT& a,
                                const IC& dims_in,
                                const IC& dims_loc,
                                const IC& offsets_in,
                                const IC& mem_dims_in,
                                const IC& mem_dims_loc,
                                const IC& mem_offsets_in)
      : container_proxy<CT>(a)
  {
    slab_rank = dims_in.size();
    for (int i = 0; i < dims_in.size(); ++i)
      slab_dims[i] = static_cast<hsize_t>(dims_in[i]);
    for (int i = 0; i < dims_loc.size(); ++i)
      slab_dims_local[i] = static_cast<hsize_t>(dims_loc[i]);
    for (int i = 0; i < dims_in.size(); ++i)
      slab_offset[i] = static_cast<hsize_t>(offsets_in[i]);

    mem_rank = mem_dims_in.size();
    for (int i = 0; i < mem_dims_in.size(); ++i)
      mem_dims[i] = static_cast<hsize_t>(mem_dims_in[i]);
    for (int i = 0; i < mem_dims_loc.size(); ++i)
      mem_dims_local[i] = static_cast<hsize_t>(mem_dims_loc[i]);
    for (int i = 0; i < mem_dims_in.size(); ++i)
      mem_offset[i] = static_cast<hsize_t>(mem_offsets_in[i]);
    if (element_size > 1)
    {
      slab_dims[slab_rank]       = element_size;
      slab_dims_local[slab_rank] = element_size;
      slab_offset[slab_rank]     = 0;
      slab_rank += 1;
      mem_dims[mem_rank]       = element_size;
      mem_dims_local[mem_rank] = element_size;
      mem_offset[mem_rank]     = 0;
      mem_rank += 1;
    }

    use_slab = true;
  }

  /** return the size of the i-th dimension
   * @param i dimension
   */
  inline hsize_t size(int i) const { return (i > MAXDIM) ? 0 : slab_dims[i]; }

  inline void change_shape() { this->resize(slab_dims.data(), MAXDIM + 1); }
};

template<typename CT, unsigned MAXDIM>
struct h5data_proxy<double_hyperslab_proxy<CT, MAXDIM>>
{
  double_hyperslab_proxy<CT, MAXDIM>& ref;
  h5data_proxy(double_hyperslab_proxy<CT, MAXDIM>& a) {}
  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref.use_slab)
    {
      return h5d_read(grp, aname.c_str(), ref.slab_rank, ref.slab_dims.data(), ref.slab_dims_local.data(),
                      ref.slab_offset.data(), ref.mem_rank, ref.mem_dims.data(), ref.mem_dims_local.data(),
                      ref.mem_offset.data(), ref.data(), xfer_plist);
    }
    else
    {
      int rank = ref.slab_rank;
      if (!get_space(grp, aname, rank, ref.slab_dims.data(), true))
      {
        ref.change_shape();
      }
      return h5d_read(grp, aname, ref.data(), xfer_plist);
    }
  }
  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref.use_slab)
    {
      return h5d_write(grp, aname.c_str(), ref.slab_rank, ref.slab_dims.data(), ref.slab_dims_local.data(),
                       ref.slab_offset.data(), ref.mem_rank, ref.mem_dims.data(), ref.mem_dims_local.data(),
                       ref.mem_offset.data(), ref.data(), xfer_plist);
    }
    else
    {
      return h5d_write(grp, aname.c_str(), ref.slab_rank, ref.slab_dims.data(), ref.data(), xfer_plist);
    }
  }
};
} // namespace qmcplusplus
#endif
