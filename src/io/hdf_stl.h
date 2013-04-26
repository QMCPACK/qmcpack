//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_HDF_STL_INTERFACE_H
#define QMCPLUSPLUS_HDF_STL_INTERFACE_H
#include <vector>
#include <sstream>
#include <bitset>
#include <deque>
namespace qmcplusplus
{
/** specialization for vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T> struct h5data_proxy<vector<T> >
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef vector<T> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(dims[0]);
    return h5d_read(grp,aname,get_address(&ref_[0]),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(&ref_[0]),xfer_plist);
  }
};

/** specialization for bitset<N>
 */
template<std::size_t N>
struct h5data_proxy<std::bitset<N> >
{
  typedef std::bitset<N> ArrayType_t;
  ArrayType_t& ref;

  h5data_proxy<ArrayType_t>(ArrayType_t& a): ref(a)
  {
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    unsigned long c=ref.to_ulong();
    h5data_proxy<unsigned long> hc(c);
    return hc.write(grp,aname,xfer_plist);
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    unsigned long c=ref.to_ulong();
    h5data_proxy<unsigned long> hc(c);
    if(hc.read(grp,aname,xfer_plist))
    {
      ref=c;
      return true;
    }
    else
      return false;
  }
};


/** Specialization for string */
template<> struct h5data_proxy<std::string>
{

  typedef std::string ArrayType_t;
  ArrayType_t& ref;

  h5data_proxy<ArrayType_t>(ArrayType_t& a): ref(a) { }

  inline bool write(hid_t grp,const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    hid_t str80 = H5Tcopy(H5T_C_S1);
    H5Tset_size(str80,ref.size());
    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  H5Dcreate(grp, aname.c_str(), str80, dataspace, H5P_DEFAULT);
    herr_t ret = H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, xfer_plist,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return ret != -1;
  }

  inline bool read(hid_t grp,const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    hid_t dataset = H5Dopen(grp,aname.c_str());
    if (dataset > -1)
    {
      hid_t datatype=H5Dget_type(dataset);
      hsize_t dim_out;
      if(datatype == H5T_NATIVE_CHAR)
      {
        hid_t dataspace = H5Dget_space(dataset);
        hid_t status = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
        H5Sclose(dataspace);
      }
      else
      {
        dim_out=H5Tget_size(datatype);
      }
      ref.resize(dim_out);
      herr_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, xfer_plist,&(ref[0]));
      // Erase trailing null character
      ref.erase (dim_out-1, 1);
      H5Tclose(datatype);
      H5Dclose(dataset);
      return ret != -1;
    }
    return false;
  }
};

template<>
struct h5data_proxy<std::ostringstream>
{
  typedef std::ostringstream Data_t;
  Data_t& ref;

  h5data_proxy<Data_t>(Data_t& a): ref(a) { }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    std::string clone(ref.str());
    h5data_proxy<std::string> proxy(clone);
    return proxy.write(grp,aname);
  }

  inline bool read(hid_t grp, const char* name, hid_t xfer_plist=H5P_DEFAULT)
  {
    return false;
  }
};

///** i/o for deque<T>, internally use vector<T>
// */
//template<typename T> struct h5data_proxy<deque<T> >
//  : public h5_space_type<T,1>
//  {
//    using h5_space_type<T,1>::dims;
//    using h5_space_type<T,1>::get_address;
//    typedef deque<T> data_type;
//    data_type& ref_;

//    inline h5data_proxy(data_type& a): ref_(a) { dims[0]=ref_.size(); }

//    inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
//    {
//      vector<T> temp(ref_.size());
//      if(!h5d_getspace(grp,aname,temp.size(),dims)) {temp.resize(dims[0]);}
//      if(h5d_read(grp,aname,get_address(&temp[0]),xfer_plist))
//      {
//        ref_.resize(temp.size());
//        ref_.assign(temp.begin(),temp.end());
//        return true;
//      }
//      else
//        return false;
//    }

//    inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
//    {
//      vector<T> temp(ref_.size());
//      temp.assign(ref_.begin(),ref_.end());
//      return h5d_write(grp,aname.c_str(),temp.size(),dims,get_address(&temp[0]),xfer_plist);
//    }
//  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: hdf_stl.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/
