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

  hid_t h1        = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(h1);
  int rank        = H5Sget_simple_extent_ndims(dataspace);

  bool success = false;
  if (h1 >= 0 && dataspace >= 0 && rank >= 0)
  {
    // check if the rank is sufficient to hold the data type
    if (rank < TSpaceType::rank)
      throw std::runtime_error(aname + " dataset is too small for the requested data type");
    else
    {
      std::vector<hsize_t> sizes_in(rank);
      int status_n = H5Sget_simple_extent_dims(dataspace, sizes_in.data(), NULL);

      // check if the lowest dimensions match the data type
      int user_rank   = rank - TSpaceType::added_rank();
      bool size_match = true;
      for (int dim = user_rank, dim_type = 0; dim < rank; dim++, dim_type++)
        if (sizes_in[dim] != TSpace.dims[dim_type])
          size_match = false;
      if (!size_match)
        throw std::runtime_error("The lower dimensions (container element type) of " + aname +
                                 " dataset do not match the requested data type");
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
  {
    H5Dclose(h1);
    return false;
  }
  H5Dclose(h1);
  return true;
}

template<typename T>
inline bool h5d_check_type(hid_t grp, const std::string& aname)
{
  if (grp < 0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  T temp(0);
  hid_t h5d_type_id = get_h5_datatype(temp);
  hid_t datatype    = H5Dget_type(h1);
  if (datatype == H5I_INVALID_HID)
    throw std::runtime_error(aname + " dataset either does not exist or there was an error determining its type.");
  htri_t equality_check = H5Tequal(datatype, h5d_type_id);
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
    hid_t dataset   = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ret             = H5Dwrite(dataset, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  else
  {
    ret = H5Dwrite(h1, h5d_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, first);
  }
  H5Dclose(h1);
  return ret != -1;
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
  //herr_t ret = H5Dread(h1, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);

  hid_t dataspace = H5Dget_space(h1);
  hid_t memspace  = H5Screate_simple(ndims, counts, NULL);
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
  herr_t ret      = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

  hid_t h5d_type_id = get_h5_datatype(*first);
  ret               = H5Dread(h1, h5d_type_id, memspace, dataspace, xfer_plist, first);

  H5Sclose(dataspace);
  H5Sclose(memspace);

  H5Dclose(h1);
  return ret != -1;
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
  hid_t filespace, memspace;
  herr_t ret = -1;

  const std::vector<hsize_t> ones(ndims, 1);
  if (h1 < 0) //missing create one
  {
    hid_t dataspace = H5Screate_simple(ndims, gcounts, NULL);
    hid_t dataset   = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t filespace = H5Dget_space(dataset);
    ret             = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

    hid_t memspace = H5Screate_simple(ndims, counts, NULL);
    ret            = H5Dwrite(dataset, h5d_type_id, memspace, filespace, xfer_plist, first);

    H5Dclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  else
  {
    filespace = H5Dget_space(h1);
    ret       = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

    memspace = H5Screate_simple(ndims, counts, NULL);
    ret      = H5Dwrite(h1, h5d_type_id, memspace, filespace, xfer_plist, first);

    H5Sclose(filespace);
    H5Dclose(memspace);
  }
  H5Dclose(h1);
  return ret != -1;
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
  if (ndims != H5Sget_simple_extent_ndims(dataspace))
    throw std::runtime_error(aname + " dataspace does not match ");
  // check gcounts???
  const std::vector<hsize_t> ones(std::max(ndims, mem_ndims), 1);
  herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

  hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
  herr_t mem_ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);

  hid_t h5d_type_id = get_h5_datatype(*first);
  ret               = H5Dread(h1, h5d_type_id, memspace, dataspace, xfer_plist, first);

  H5Sclose(dataspace);
  H5Sclose(memspace);

  H5Dclose(h1);
  return ret != -1;
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
  std::cout << " h5d_write: " << mem_ndims << " " << *mem_gcounts << " " << *(mem_gcounts + 1) << " "
            << *(mem_gcounts + 2) << " " << *mem_counts << " " << *(mem_counts + 1) << " " << *(mem_counts + 2) << " "
            << *mem_offsets << " " << *(mem_offsets + 1) << " " << *(mem_offsets + 2) << " " << std::endl;
  hid_t h5d_type_id = get_h5_datatype(*first);
  hid_t h1          = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
  herr_t ret        = -1;

  const std::vector<hsize_t> ones(std::max(ndims, mem_ndims), 1);
  if (h1 < 0) //missing create one
  {
    hid_t dataspace = H5Screate_simple(ndims, gcounts, NULL);
    hid_t dataset   = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t filespace = H5Dget_space(dataset);
    ret             = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

    hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
    ret            = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);
    ret            = H5Dwrite(dataset, h5d_type_id, memspace, filespace, xfer_plist, first);

    H5Dclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  else
  {
    hid_t filespace = H5Dget_space(h1);
    ret             = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);

    hid_t memspace = H5Screate_simple(mem_ndims, mem_gcounts, NULL);
    ret            = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offsets, NULL, ones.data(), mem_counts);
    ret            = H5Dwrite(h1, h5d_type_id, memspace, filespace, xfer_plist, first);

    H5Sclose(filespace);
    H5Dclose(memspace);
  }
  H5Dclose(h1);
  return ret != -1;
}

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
  //app_log()<<omp_get_thread_num()<<"  h5d_append  group = "<<grp<<"  name = "<<aname.c_str()<< std::endl;
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
    // create dataset property list
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    // set layout (chunked, contiguous)
    hid_t sl = H5Pset_layout(p, H5D_CHUNKED);
    // set chunk size
    hid_t cs = H5Pset_chunk(p, ndims, chunk_dims.data());
    // create the dataset
    dataset = H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, p, H5P_DEFAULT);
    // create memory dataspace, size of current buffer
    memspace = H5Screate_simple(ndims, dims, NULL);
    // write the data for the first time
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, xfer_plist, first);
    // update the "file pointer"
    current = dims[0];

    //app_log()<<"  creating dataset"<< std::endl;
    //if(dataspace<0) app_log()<<"    dataspace could not be created"<< std::endl;
    //else            app_log()<<"    dataspace creation successful"<< std::endl;
    //if(p<0)         app_log()<<"    property list could not be created"<< std::endl;
    //else            app_log()<<"    property list creation successful"<< std::endl;
    //if(sl<0)        app_log()<<"    layout could not be set"<< std::endl;
    //else            app_log()<<"    layout set successfully"<< std::endl;
    //if(cs<0)        app_log()<<"    chunk size could not be set"<< std::endl;
    //else            app_log()<<"    chunk size set successfully"<< std::endl;
    //if(dataset<0)   app_log()<<"    dataset could not be created"<< std::endl;
    //else            app_log()<<"    dataset creation successful"<< std::endl;
    //if(memspace<0)  app_log()<<"    memspace could not be created"<< std::endl;
    //else            app_log()<<"    memspace creation successful"<< std::endl;
    //if(ret<0)       app_log()<<"    data could not be written"<< std::endl;
    //else            app_log()<<"    data write successful"<< std::endl;
    //H5Eprint(NULL);

    // close the property list
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
    herr_t he = H5Dset_extent(dataset, end.data());
    //get the corresponding dataspace (filespace)
    dataspace = H5Dget_space(dataset);
    //set the extent
    herr_t hse = H5Sset_extent_simple(dataspace, ndims, end.data(), max_dims.data());
    //select hyperslab/slice of multidimensional data for appended write
    const std::vector<hsize_t> ones(ndims, 1);
    herr_t hsh = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), NULL, ones.data(), dims);
    //create memory space describing current data block
    memspace = H5Screate_simple(ndims, dims, NULL);
    //append the datablock to the dataset
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, H5P_DEFAULT, first);
    // update the "file pointer"
    current = end[0];

    //app_log()<<"  appending to dataset"<< std::endl;
    //app_log()<<"      ndims = "<<ndims<< std::endl;
    //app_log()<<"      dims  = "<<dims[0]<<" "<<dims[1]<< std::endl;
    //app_log()<<"      start = "<<start[0]<<" "<<start[1]<< std::endl;
    //app_log()<<"      end   = "<<end[0]<<" "<<end[1]<< std::endl;
    //app_log()<<"      current = "<<current<<" "<<&current<< std::endl;
    //if(hse<0)       app_log()<<"    set_extent failed"<< std::endl;
    //else            app_log()<<"    set_extent successful"<< std::endl;
    //if(hsh<0)       app_log()<<"    select_hyperslab failed"<< std::endl;
    //else            app_log()<<"    select_hyperslab successful"<< std::endl;
    //if(he<0)        app_log()<<"    extend failed"<< std::endl;
    //else            app_log()<<"    extend successful"<< std::endl;
    //if(dataspace<0) app_log()<<"    dataspace could not be gotten"<< std::endl;
    //else            app_log()<<"    dataspace get successful"<< std::endl;
    //if(memspace<0)  app_log()<<"    memspace could not be created"<< std::endl;
    //else            app_log()<<"    memspace creation successful"<< std::endl;
    //if(ret<0)       app_log()<<"    data could not be written"<< std::endl;
    //else            app_log()<<"    data write successful"<< std::endl;
    //H5Eprint(NULL);
  }
  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return ret != -1;
}

} // namespace qmcplusplus
#endif
