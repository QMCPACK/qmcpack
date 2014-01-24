//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file hdf_dataproxy.h
 * @brief define h5data_proxy to handle i/o
 */
#ifndef QMCPLUSPLUS_HDF_H5_DATAPROXY_H
#define QMCPLUSPLUS_HDF_H5_DATAPROXY_H

#include <io/hdf_datatype.h>
#include <io/hdf_dataspace.h>

namespace qmcplusplus
{

/** free function to check dimension
 * @param grp group id
 * @param aname name of the dataspace
 * @param rank rank of the multi-dimensional array
 * @param dims[rank] size for each direction
 * @return true if the dims is the same as the dataspace
 */
inline bool get_space(hid_t grp, const std::string& aname, int rank, hsize_t* dims, bool fatalerr=true)
{
  hid_t h1 = H5Dopen(grp, aname.c_str());
  hid_t dataspace = H5Dget_space(h1);
  int rank_in = H5Sget_simple_extent_ndims(dataspace);
  if(rank < rank_in && fatalerr)
  {
    APP_ABORT(aname + " dataspace does not match ");
  }
  vector<hsize_t> dims_in(rank);
  int status_n = H5Sget_simple_extent_dims(dataspace, &dims_in[0], NULL);
  H5Dclose(h1);
  bool thesame=true;
  for(int i=0; i<rank; ++i)
  {
    thesame &= (dims_in[i]==dims[i]);
    dims[i]=dims_in[i];
  }
  return thesame;
}

/** return true, if successful */
template<typename T>
inline bool h5d_read(hid_t grp, const std::string& aname, T* first, hid_t xfer_plist)
{
  if(grp<0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str());
  if(h1<0)
    return false;
  hid_t h5d_type_id=get_h5_datatype(*first);
  herr_t ret = H5Dread(h1, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
  H5Dclose(h1);
  return ret != -1;
}

template<typename T>
inline bool h5d_write(hid_t grp, const std::string& aname, hsize_t ndims, const hsize_t* dims, const T* first
                      , hid_t xfer_plist)
{
  if(grp<0)
    return true;
  hid_t h5d_type_id=get_h5_datatype(*first);
  hid_t h1 = H5Dopen(grp, aname.c_str());
  herr_t ret=-1;
  if(h1<0) //missing create one
  {
    hid_t dataspace  = H5Screate_simple(ndims, dims, NULL);
    hid_t dataset =  H5Dcreate(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT);
    ret = H5Dwrite(dataset, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  else
  {
    ret = H5Dwrite(h1, h5d_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,first);
  }
  H5Dclose(h1);
  return ret != -1;
}

/** return true, if successful */
template<typename T>
bool h5d_read(hid_t grp, const std::string& aname, hsize_t ndims, 
    const hsize_t* gcounts, const hsize_t* counts, const hsize_t* offsets,
    T* first , hid_t xfer_plist)
{
  if(grp<0)
    return true;
  hid_t h1 = H5Dopen(grp, aname.c_str());
  if(h1<0)
    return false;
  //herr_t ret = H5Dread(h1, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
  
  hid_t dataspace = H5Dget_space(h1);
  hid_t memspace = H5Screate_simple(ndims, counts, NULL);
  herr_t ret = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offsets,NULL,counts,NULL);

  hid_t h5d_type_id=get_h5_datatype(*first);
  ret = H5Dread(h1, h5d_type_id, memspace, dataspace, xfer_plist, first);

  H5Sclose(dataspace);
  H5Sclose(memspace);

  H5Dclose(h1);
  return ret != -1;
}

template<typename T>
inline bool h5d_write(hid_t grp, const std::string& aname, hsize_t ndims, 
    const hsize_t* gcounts, const hsize_t* counts, const hsize_t* offsets,
    const T* first , hid_t xfer_plist)
{
  if(grp<0)
    return true;
  hid_t h5d_type_id=get_h5_datatype(*first);
  hid_t h1 = H5Dopen(grp, aname.c_str());
  herr_t ret=-1;
  if(h1<0) //missing create one
  {
    hid_t dataspace=H5Screate_simple(ndims,gcounts,NULL);
    hid_t dataset=H5Dcreate(grp, aname.c_str(),h5d_type_id, dataspace, H5P_DEFAULT);

    hid_t filespace=H5Dget_space(dataspace);
    ret=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offsets,NULL,counts,NULL); 

    hid_t memspace=H5Screate_simple(ndims,counts,NULL);
    ret=H5Dwrite(dataset,h5d_type_id,memspace,filespace,xfer_plist,first);

    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  H5Dclose(h1);
  return ret != -1;
}

template<typename T>
inline bool h5d_append(hid_t grp, const std::string& aname, hsize_t& current,
                       hsize_t ndims, const hsize_t* dims, const T* first,
                       hsize_t chunk_size=1,hid_t xfer_plist=H5P_DEFAULT)
{
  //app_log()<<omp_get_thread_num()<<"  h5d_append  group = "<<grp<<"  name = "<<aname.c_str()<<endl;
  if(grp<0)
    return true;
  hid_t h5d_type_id=get_h5_datatype(*first);
  hid_t dataspace;
  hid_t memspace;
  hid_t dataset = H5Dopen(grp, aname.c_str());
  hsize_t max_dims[ndims];
  max_dims[0]   = H5S_UNLIMITED;
  for(int d=1; d<ndims; ++d)
    max_dims[d] = dims[d];
  herr_t ret=-1;
  if(dataset<0) //missing create one
  {
    //set file pointer
    current = 0;
    // set max and chunk dims
    hsize_t chunk_dims[ndims];
    chunk_dims[0] = chunk_size;
    for(int d=1; d<ndims; ++d)
      chunk_dims[d] = dims[d];
    // create a dataspace sized to the current buffer
    dataspace = H5Screate_simple(ndims, dims, max_dims);
    // create dataset property list
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    // set layout (chunked, contiguous)
    hid_t sl = H5Pset_layout(p, H5D_CHUNKED);
    // set chunk size
    hid_t cs = H5Pset_chunk(p, ndims, chunk_dims);
    // create the dataset
    dataset = H5Dcreate2(grp, aname.c_str(), h5d_type_id, dataspace, H5P_DEFAULT, p, H5P_DEFAULT);
    // create memory dataspace, size of current buffer
    memspace = H5Screate_simple(ndims, dims, NULL);
    // write the data for the first time
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, xfer_plist, first);
    // update the "file pointer"
    current = dims[0];

    //app_log()<<"  creating dataset"<<endl;
    //if(dataspace<0) app_log()<<"    dataspace could not be created"<<endl;
    //else            app_log()<<"    dataspace creation successful"<<endl;
    //if(p<0)         app_log()<<"    property list could not be created"<<endl;
    //else            app_log()<<"    property list creation successful"<<endl;
    //if(sl<0)        app_log()<<"    layout could not be set"<<endl;
    //else            app_log()<<"    layout set successfully"<<endl;
    //if(cs<0)        app_log()<<"    chunk size could not be set"<<endl;
    //else            app_log()<<"    chunk size set successfully"<<endl;
    //if(dataset<0)   app_log()<<"    dataset could not be created"<<endl;
    //else            app_log()<<"    dataset creation successful"<<endl;
    //if(memspace<0)  app_log()<<"    memspace could not be created"<<endl;
    //else            app_log()<<"    memspace creation successful"<<endl;
    //if(ret<0)       app_log()<<"    data could not be written"<<endl;
    //else            app_log()<<"    data write successful"<<endl;
    //H5Eprint(NULL);

    // close the property list
    H5Pclose(p);
  }
  else
  {
    // new end of file
    hsize_t start[ndims];
    hsize_t end[ndims];
    for(int d=1; d<ndims; ++d)
    {
      start[d] = 0;
      end[d]   = dims[d];
    }
    start[0] = current;
    end[0] = start[0]+dims[0];
    //extend the dataset (file)
    herr_t he  = H5Dextend(dataset, end);
    //get the corresponding dataspace (filespace)
    dataspace  = H5Dget_space(dataset);
    //set the extent
    herr_t hse = H5Sset_extent_simple(dataspace, ndims, end, max_dims);
    //select hyperslab/slice of multidimensional data for appended write
    herr_t hsh = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, dims, NULL);
    //create memory space describing current data block
    memspace   = H5Screate_simple(ndims, dims, NULL);
    //append the datablock to the dataset
    ret = H5Dwrite(dataset, h5d_type_id, memspace, dataspace, H5P_DEFAULT, first);
    // update the "file pointer"
    current = end[0];

    //app_log()<<"  appending to dataset"<<endl;
    //app_log()<<"      ndims = "<<ndims<<endl;
    //app_log()<<"      dims  = "<<dims[0]<<" "<<dims[1]<<endl;
    //app_log()<<"      start = "<<start[0]<<" "<<start[1]<<endl;
    //app_log()<<"      end   = "<<end[0]<<" "<<end[1]<<endl;
    //app_log()<<"      current = "<<current<<" "<<&current<<endl;
    //if(hse<0)       app_log()<<"    set_extent failed"<<endl;
    //else            app_log()<<"    set_extent successful"<<endl;
    //if(hsh<0)       app_log()<<"    select_hyperslab failed"<<endl;
    //else            app_log()<<"    select_hyperslab successful"<<endl;
    //if(he<0)        app_log()<<"    extend failed"<<endl;
    //else            app_log()<<"    extend successful"<<endl;
    //if(dataspace<0) app_log()<<"    dataspace could not be gotten"<<endl;
    //else            app_log()<<"    dataspace get successful"<<endl;
    //if(memspace<0)  app_log()<<"    memspace could not be created"<<endl;
    //else            app_log()<<"    memspace creation successful"<<endl;
    //if(ret<0)       app_log()<<"    data could not be written"<<endl;
    //else            app_log()<<"    data write successful"<<endl;
    //H5Eprint(NULL);
  }
  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return ret != -1;
}


/** generic h5data_proxy<T> for hdf5 native type T
*/
template<typename T> struct h5data_proxy
    : public h5_space_type<T,1>
{
  typedef T data_type;
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=1;
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_read(grp,aname,get_address(&ref_),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname,hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(&ref_),xfer_plist);
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2764 $   $Date: 2008-06-26 10:21:31 -0500 (Thu, 26 Jun 2008) $
 * $Id: hdf_dataspace.h 2764 2008-06-26 15:21:31Z jnkim $
 ***************************************************************************/
