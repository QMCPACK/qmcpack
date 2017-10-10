//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_PETE_TRAITS_H
#define QMCPLUSPLUS_HDF_PETE_TRAITS_H

#include <io/hdf_dataproxy.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Utilities/PooledData.h>

namespace qmcplusplus
{

/** specialization for TinyVector<T,D>
 */
template<typename T, unsigned D> struct h5data_proxy<TinyVector<T,D> >
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef TinyVector<T,D> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=D;
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_read(grp,aname,get_address(ref_.data()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()),xfer_plist);
  }
};

/** specialization for Tensor<T,D>
 */
template<typename T, unsigned D> struct h5data_proxy<Tensor<T,D> >
    : public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef Tensor<T,D> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=D;
    dims[1]=D;
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_read(grp,aname,get_address(ref_.data()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()),xfer_plist);
  }
};


/** specialization for Vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T> struct h5data_proxy<Vector<T> >: public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef Vector<T> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(dims[0]);
    return h5d_read(grp,aname,get_address(ref_.data()),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()),xfer_plist);
  }
};


template<typename T>
struct h5data_proxy<Matrix<T> >: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef Matrix<T> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.rows();
    dims[1]=ref_.cols();
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(dims[0],dims[1]);
    return h5d_read(grp,aname,get_address(ref_.data()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()),xfer_plist);
  }
};


template<typename T, unsigned D>
struct h5data_proxy<Array<T,D> >: public h5_space_type<T,D>
{
  using h5_space_type<T,D>::dims;
  using h5_space_type<T,D>::get_address;
  typedef Array<T,D> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    for(int i=0; i<D; ++i)
      dims[i]=ref_.size(i);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(dims);
    return h5d_read(grp,aname,get_address(ref_.data()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()),xfer_plist);
  }
};

///** specialization for Vector<T>
// *
// * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
// */
//template<typename T> struct h5data_proxy<PooledData<T> >
//{
//  typedef PooledData<T> data_type;
//  std::vector<T>& ref_;
//  std::vector<hsize_t> dims;

//  inline h5data_proxy(data_type& a): ref_(a.myData){dims.resize(1,a.size());}

//  template<typename IC>
//  inline h5data_proxy(data_type& a, const IC& dims_in): ref_(a.myData)
//  {
//    dims.resize(dims_in.size());
//    for(int i=0; i<dims_in.size(); ++i) dims[i]=static_cast<hsize_t>(dims_in[i]);
//  }

//  template<typename IC>
//  inline h5data_proxy(std::vector<T>& a, const IC& dims_in): ref_(a)
//  {
//    dims.resize(dims_in.size());
//    for(int i=0; i<dims_in.size(); ++i) dims[i]=static_cast<hsize_t>(dims_in[i]);
//  }

//  inline hsize_t size(int i) const { return dims[i];}

//  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
//  {
//    const int rank=dims.size();
//    if(!h5d_getspace(grp,aname,rank,&dims[0]))
//    {
//      size_t ntot=dims[0];
//      for(int i=1;i<rank;++i) ntot*=dims[i];
//      ref_.resize(ntot);
//    }
//    return h5d_read(grp,aname,ref_.data(),xfer_plist);
//  }

//  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
//  {
//    return h5d_write(grp,aname.c_str(),dims.size(),&dims[0],ref_.data(),xfer_plist);
//  }
//};

}
#endif
