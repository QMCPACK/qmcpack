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

namespace qmcplusplus {

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
    for(int i=0; i<rank;++i)
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
      if(grp<0) return true;
      hid_t h1 = H5Dopen(grp, aname.c_str());
      if(h1<0) return false; 
      hid_t h5d_type_id=get_h5_datatype(*first);
      herr_t ret = H5Dread(h1, h5d_type_id, H5S_ALL, H5S_ALL, xfer_plist, first);
      H5Dclose(h1);
      return ret != -1;
    }

  template<typename T>
    inline bool h5d_write(hid_t grp, const std::string& aname, hsize_t ndims, const hsize_t* dims, const T* first
        , hid_t xfer_plist)
    {
      if(grp<0) return true;
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


  /** generic h5data_proxy<T> for hdf5 native type T
  */
  template<typename T> struct h5data_proxy
    : public h5_space_type<T,1>
  { 
    typedef T data_type;
    using h5_space_type<T,1>::dims;
    using h5_space_type<T,1>::get_address;
    data_type& ref_;

    inline h5data_proxy(data_type& a): ref_(a) { dims[0]=1; }

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
