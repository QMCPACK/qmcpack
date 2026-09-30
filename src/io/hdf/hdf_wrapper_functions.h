//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_WRAPPER_FUNCTIONS_H
#define QMCPLUSPLUS_HDF_WRAPPER_FUNCTIONS_H

/**@file hdf_wrapper_functions.h
 * @brief free template functions wrapping HDF5 calls.
 */

#include <vector>
#include "hdf_datatype.h"
#include "hdf_dataspace.h"

namespace qmcplusplus
{
/** free template function to read the (user) dimension and shape of the dataset.
 * The dimensions contributed by T is excluded.
 * @tparam T data type supported by h5_space_type
 * @param grp group id
 * @param aname name of the dataspace
 * @param sizes_out sizes of each direction. For a scalar, sizes_out.size() == 0
 * @return true if sizes_out is extracted successfully
 *
 * For example, if the data on the file is Matrix<TinyVector<std::complex<double>, 3>>
 * The dataset on the file has a rank of 2 (matrix) + 1 (TinyVector) + 1 (std::complex) + 0 (double) = 4
 * getDataShape<TinyVector<std::complex<double>, 3>> only returns the first 2 dimension
 * getDataShape<std::complex<double>> only returns the first 3 dimension
 * getDataShape<double> returns all the 4 dimension
 */
template<typename T, typename IT>
inline bool getDataShape(hid_t grp, const std::string& aname, std::vector<IT>& sizes_out)
{
  using TSpaceType = h5_space_type<T, 0>;
  TSpaceType TSpace;

  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    return false;
  hid_t dataspace = H5Dget_space(h1);
  if (dataspace < 0)
  {
    H5Dclose(h1);
    return false;
  }
  int rank = H5Sget_simple_extent_ndims(dataspace);

  bool success = false;
  if (rank >= 0)
  {
    // check if the rank is sufficient to hold the data type
    if (rank < TSpaceType::rank)
    {
      H5Sclose(dataspace);
      H5Dclose(h1);
      throw std::runtime_error(aname + " dataset is too small for the requested data type");
    }
    else
    {
      std::vector<hsize_t> sizes_in(rank);
      if (H5Sget_simple_extent_dims(dataspace, sizes_in.data(), NULL) < 0)
      {
        H5Sclose(dataspace);
        H5Dclose(h1);
        return false;
      }

      // check if the lowest dimensions match the data type
      int user_rank   = rank - TSpaceType::added_rank();
      bool size_match = true;
      for (int dim = user_rank, dim_type = 0; dim < rank; dim++, dim_type++)
        if (sizes_in[dim] != TSpace.dims[dim_type])
          size_match = false;
      if (!size_match)
      {
        H5Sclose(dataspace);
        H5Dclose(h1);
        throw std::runtime_error("The lower dimensions (container element type) of " + aname +
                                 " dataset do not match the requested data type");
      }
      else
      {
        // save the sizes of each directions excluding dimensions contributed by the data type
        sizes_out.resize(user_rank);
        for (int dim = 0; dim < user_rank; dim++)
          sizes_out[dim] = static_cast<IT>(sizes_in[dim]);
        success = true;
      }
    }
  }

  H5Sclose(dataspace);
  H5Dclose(h1);
  return success;
}

/** free function to check dimension
 * @param grp group id
 * @param aname name of the dataspace
 * @param rank rank of the multi-dimensional array
 * @param dims[rank] size for each direction, return the actual size on file
 * @return true if the dims is the same as the dataspace
 */
template<typename T>
inline bool checkShapeConsistency(hid_t grp, const std::string& aname, int rank, hsize_t* dims)
{
  using TSpaceType = h5_space_type<T, 0>;

  std::vector<hsize_t> dims_in;
  if (getDataShape<T>(grp, aname, dims_in))
  {
    const int user_rank = rank - TSpaceType::added_rank();
    if (dims_in.size() != user_rank)
      throw std::runtime_error(aname + " dataspace rank does not match\n");

    bool is_same = true;
    for (int i = 0; i < user_rank; ++i)
    {
      is_same &= (dims_in[i] == dims[i]);
      dims[i] = dims_in[i];
    }
    return is_same;
  }
  else
    return false;
}

/** return true, if successful */
template<typename T>
inline bool h5d_read(hid_t grp, const std::string& aname, T* first, hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    return false;
  hid_t h5d_type_id = get_h5_datatype(*first);
  herr_t ret        = H5Dread(h1, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
  H5Dclose(h1);
  return ret != -1;
}

inline bool h5d_check_existence(hid_t grp, const std::string& aname)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    return false;
  H5Dclose(h1);
  return true;
}

template<typename T>
inline bool h5d_check_type(hid_t grp, const std::string& aname)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    throw std::runtime_error(aname + " dataset either does not exist or there was an error opening it.");
  T temp(0);
  hid_t h5d_type_id = get_h5_datatype(temp);
  hid_t datatype    = H5Dget_type(h1);
  if (datatype < 0)
  {
    H5Dclose(h1);
    throw std::runtime_error(aname + " dataset either does not exist or there was an error determining its type.");
  }
  htri_t equality_check = H5Tequal(datatype, h5d_type_id);
  H5Tclose(datatype);
  H5Dclose(h1);
  switch (equality_check)
  {
  case 1:
    return true;
  case 0:
    return false;
  default:
    throw std::runtime_error("Type comparison attempted with an invalid type or nonexistent dataset " + aname);
  }
}

template<typename T>
inline bool h5d_write(hid_t grp,
                      const std::string& aname,
                      hsize_t ndims,
                      const hsize_t* dims,
                      const T* first,
                      hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h5d_type_id = get_h5_datatype(*first);
  hid_t h1          = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  herr_t ret        = -1;
  if (h1 < 0) //missing create one
  {
    hid_t dataspace = H5Screate_simple(ndims, dims, NULL);
    if (dataspace < 0)
      return false;
    hid_t dataset = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0)
    {
      H5Sclose(dataspace);
      return false;
    }
    ret = H5Dwrite(dataset, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  else
  {
    ret = H5Dwrite(h1, h5d_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, first);
    H5Dclose(h1);
  }
  return ret >= 0;
}

// MAM: Make new h5d_read/write overloads that take more parameters which allow you to
// use a hyperslab on the memory space too. Then use it through a template specialization of
// hyperslap for multi::array that allows you to define the memory space hyperslab using
// shape and strides.

/** return true, if successful */
template<typename T>
bool h5d_read(hid_t grp,
              const std::string& aname,
              hsize_t ndims,
              const hsize_t* gcounts,
              const hsize_t* counts,
              const hsize_t* offsets,
              T* first,
              hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    return false;

  hid_t dataspace = H5Dget_space(h1);
  if (dataspace < 0)
  {
    H5Dclose(h1);
    return false;
  }

  hid_t memspace = H5Screate_simple(ndims, counts, NULL);
  if (memspace < 0)
  {
    H5Sclose(dataspace);
    H5Dclose(h1);
    return false;
  }

  // According to the HDF5 manual (https://support.hdfgroup.org/HDF5/doc/RM/H5S/H5Sselect_hyperslab.htm)
  // , the fifth argument (count) means the number of hyper-slabs to select along each dimensions
  // while the sixth argument (block) is the size of each hyper-slab.
  // To write a single hyper-slab of size counts in a dataset, we call
  // H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
  // The vector "ones" means we want to write one hyper-slab (block) along each dimensions.
  // The result is equivalent to calling
  // H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
  // , but it implies writing count hyper-slabs along each dimension and each hyper-slab is of size one.
  const std::vector<hsize_t> ones(ndims, 1);
  herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
  if (ret >= 0)
  {
    hid_t h5d_type_id = get_h5_datatype(*first);
    ret               = H5Dread(h1, h5d_type_id, memspace, dataspace, xfer_plist, first);
  }

  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(h1);
  return ret >= 0;
}


template<typename T>
inline bool h5d_write(hid_t grp,
                      const std::string& aname,
                      hsize_t ndims,
                      const hsize_t* gcounts,
                      const hsize_t* counts,
                      const hsize_t* offsets,
                      const T* first,
                      hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h5d_type_id = get_h5_datatype(*first);
  hid_t h1          = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  herr_t ret        = -1;

  const std::vector<hsize_t> ones(ndims, 1);
  if (h1 < 0) //missing create one
  {
    hid_t dataspace = H5Screate_simple(ndims, gcounts, NULL);
    if (dataspace < 0)
      return false;
    hid_t dataset = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0)
    {
      H5Sclose(dataspace);
      return false;
    }

    hid_t filespace = H5Dget_space(dataset);
    if (filespace < 0)
    {
      H5Dclose(dataset);
      H5Sclose(dataspace);
      return false;
    }
    ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
    if (ret >= 0)
    {
      hid_t memspace = H5Screate_simple(ndims, counts, NULL);
      if (memspace >= 0)
      {
        ret = H5Dwrite(dataset, h5d_type_id, memspace, filespace, xfer_plist, first);
        H5Dclose(memspace);
      }
      else
        ret = -1;
    }

    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  else
  {
    hid_t filespace = H5Dget_space(h1);
    if (filespace < 0)
    {
      H5Dclose(h1);
      return false;
    }
    ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
    if (ret >= 0)
    {
      hid_t memspace = H5Screate_simple(ndims, counts, NULL);
      if (memspace >= 0)
      {
        ret = H5Dwrite(h1, h5d_type_id, memspace, filespace, xfer_plist, first);
        H5Dclose(memspace);
      }
      else
        ret = -1;
    }

    H5Sclose(filespace);
    H5Dclose(h1);
  }
  return ret >= 0;
}

/** return true, if successful */
template<typename T>
bool h5d_read(hid_t grp,
              const std::string& aname,
              hsize_t ndims,
              const hsize_t* gcounts,
              const hsize_t* counts,
              const hsize_t* offsets,
              hsize_t mem_ndims,
              const hsize_t* mem_gcounts,
              const hsize_t* mem_counts,
              const hsize_t* mem_offsets,
              T* first,
              hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  if (h1 < 0)
    return false;

  hid_t dataspace = H5Dget_space(h1);
  if (dataspace < 0)
  {
    H5Dclose(h1);
    return false;
  }
  if (ndims != H5Sget_simple_extent_ndims(dataspace))
  {
    H5Sclose(dataspace);
    H5Dclose(h1);
    throw std::runtime_error(aname + " dataspace does not match ");
  }
  const std::vector<hsize_t> ones(std::max(ndims, mem_ndims), 1);
  herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
  if (ret >= 0)
  {
    hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
    if (memspace < 0)
    {
      H5Sclose(dataspace);
      H5Dclose(h1);
      return false;
    }
    herr_t mem_ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);
    if (mem_ret >= 0)
    {
      hid_t h5d_type_id = get_h5_datatype(*first);
      ret               = H5Dread(h1, h5d_type_id, memspace, dataspace, xfer_plist, first);
    }
    else
      ret = -1;
    H5Sclose(memspace);
  }

  H5Sclose(dataspace);
  H5Dclose(h1);
  return ret >= 0;
}

template<typename T>
inline bool h5d_write(hid_t grp,
                      const std::string& aname,
                      hsize_t ndims,
                      const hsize_t* gcounts,
                      const hsize_t* counts,
                      const hsize_t* offsets,
                      hsize_t mem_ndims,
                      const hsize_t* mem_gcounts,
                      const hsize_t* mem_counts,
                      const hsize_t* mem_offsets,
                      const T* first,
                      hid_t xfer_plist)
{
  if (grp < 0)
    return true;
  hid_t h5d_type_id = get_h5_datatype(*first);
  hid_t h1          = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  herr_t ret        = -1;

  const std::vector<hsize_t> ones(std::max(ndims, mem_ndims), 1);
  if (h1 < 0) //missing create one
  {
    hid_t dataspace = H5Screate_simple(ndims, gcounts, NULL);
    if (dataspace < 0)
      return false;
    hid_t dataset = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0)
    {
      H5Sclose(dataspace);
      return false;
    }

    hid_t filespace = H5Dget_space(dataset);
    if (filespace < 0)
    {
      H5Dclose(dataset);
      H5Sclose(dataspace);
      return false;
    }
    ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
    if (ret >= 0)
    {
      hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
      if (memspace >= 0)
      {
        herr_t mem_ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);
        if (mem_ret >= 0)
          ret = H5Dwrite(dataset, h5d_type_id, memspace, filespace, xfer_plist, first);
        else
          ret = -1;
        H5Dclose(memspace);
      }
      else
        ret = -1;
    }

    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  else
  {
    hid_t filespace = H5Dget_space(h1);
    if (filespace < 0)
    {
      H5Dclose(h1);
      return false;
    }
    ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
    if (ret >= 0)
    {
      hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
      if (memspace >= 0)
      {
        herr_t mem_ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);
        if (mem_ret >= 0)
          ret = H5Dwrite(h1, h5d_type_id, memspace, filespace, xfer_plist, first);
        else
          ret = -1;
        H5Dclose(memspace);
      }
      else
        ret = -1;
    }

    H5Sclose(filespace);
    H5Dclose(h1);
  }
  return ret >= 0;
}

/** Append to a hdf5 dataset
 *  @param[in]      grp         an invalid grp id will result in
 *                              a noop returning true
 *  @param[in]      aname       string identifier of dataset
 *  @param[in,out]  current     position in dataset, input is ignored if
 *                              dataset at aname doesn't exist
 *  @param[in]      ndims       number of dimensions
 *  @param[in]      dims        ptr to dims
 *  @param[in]      first       ptr to start of data
 *  @param[in]      chunk_size  chunk size
 *  @param[in]      xfer_plist  hid to hd5 plist
 *
 *  Assumption: we are appending slices in the first index.
 */
template<typename T>
inline bool h5d_append(hid_t grp,
                       const std::string& aname,
                       hsize_t& current,
                       hsize_t ndims,
                       const hsize_t* const dims,
                       const T* const first,
                       hsize_t chunk_size = 1,
                       hid_t xfer_plist   = H5P_DEFAULT)
{
  if (grp < 0)
    return true;
  hid_t h5d_type_id = get_h5_datatype(*first);
  hid_t dataspace;
  hid_t memspace;
  hid_t dataset = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  std::vector<hsize_t> max_dims(ndims);
  max_dims[0] = H5S_UNLIMITED;
  for (int d = 1; d < ndims; ++d)
    max_dims[d] = dims[d];
  herr_t ret = -1;
  if (dataset < 0) //missing create one
  {
    //set file pointer
    current = 0;
    // set max and chunk dims
    std::vector<hsize_t> chunk_dims(ndims);
    chunk_dims[0] = chunk_size;
    for (int d = 1; d < ndims; ++d)
      chunk_dims[d] = dims[d];
    // create a dataspace sized to the current buffer
    dataspace = H5Screate_simple(ndims, dims, max_dims.data());
    if (dataspace < 0)
      return false;
    // create dataset property list
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    if (p < 0)
    {
      H5Sclose(dataspace);
      return false;
    }
    // set layout (chunked, contiguous)
    if (H5Pset_layout(p, H5D_CHUNKED) < 0 || H5Pset_chunk(p, ndims, chunk_dims.data()) < 0)
    {
      H5Pclose(p);
      H5Sclose(dataspace);
      return false;
    }
    // create the dataset
    dataset = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, p, H5P_DEFAULT);
    if (dataset < 0)
    {
      H5Pclose(p);
      H5Sclose(dataspace);
      return false;
    }
    // create memory dataspace, size of current buffer
    memspace = H5Screate_simple(ndims, dims, NULL);
    if (memspace < 0)
    {
      H5Pclose(p);
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return false;
    }
    // write the data for the first time
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, xfer_plist, first);
    // update the "file pointer"
    current = dims[0];
    H5Pclose(p);
  }
  else
  {
    // new end of file
    std::vector<hsize_t> start(ndims);
    std::vector<hsize_t> end(ndims);
    for (int d = 1; d < ndims; ++d)
    {
      start[d] = 0;
      end[d]   = dims[d];
    }
    start[0] = current;
    end[0]   = start[0] + dims[0];
    //extend the dataset (file)
    if (H5Dset_extent(dataset, end.data()) < 0)
    {
      H5Dclose(dataset);
      return false;
    }
    //get the corresponding dataspace (filespace)
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
      H5Dclose(dataset);
      return false;
    }
    //set the extent
    if (H5Sset_extent_simple(dataspace, ndims, end.data(), max_dims.data()) < 0)
    {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return false;
    }
    //select hyperslab/slice of multidimensional data for appended write
    const std::vector<hsize_t> ones(ndims, 1);
    if (H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), NULL, ones.data(), dims) < 0)
    {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return false;
    }
    //create memory space describing current data block
    memspace = H5Screate_simple(ndims, dims, NULL);
    if (memspace < 0)
    {
      H5Sclose(dataspace);
      H5Dclose(dataset);
      return false;
    }
    //append the datablock to the dataset
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, H5P_DEFAULT, first);
    // update the "file pointer"
    current = end[0];
  }
  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return ret >= 0;
}

} // namespace qmcplusplus
#endif
