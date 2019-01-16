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

#ifndef QMCPLUSPLUS_HDF_MA_INTERFACE_H
#define QMCPLUSPLUS_HDF_MA_INTERFACE_H
#include <vector>
#include <sstream>
#include <bitset>
#include <deque>
#include "boost/multi_array.hpp"
namespace qmcplusplus
{
/** specialization for vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T> 
struct h5data_proxy<boost::multi_array<T,1>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi_array<T,1> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(boost::extents[dims[0]]);
    return h5d_read(grp,aname,get_address(ref_.origin()),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.origin()),xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi_array<T,2>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi_array<T,2> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.shape()[0];
    dims[1]=ref_.shape()[1];
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize(boost::extents[dims[0]][dims[1]]);
    return h5d_read(grp,aname,get_address(ref_.origin()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.origin()),xfer_plist);
  }
};

template<typename T> 
struct h5data_proxy<boost::multi_array_ref<T,1>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi_array_ref<T,1> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0] > 0) {  
        std::cerr<<" Error: multi_array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<ref_.shape()[0] <<std::endl; 
      }  
      return false;
    }
    return h5d_read(grp,aname,get_address(ref_.origin()),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.origin()),xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi_array_ref<T,2>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi_array_ref<T,2> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.shape()[0];
    dims[1]=ref_.shape()[1];
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0]*dims[1] > 0) {  
        std::cerr<<" Error: multi_array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<dims[1] <<" " <<ref_.shape()[0] <<" " <<ref_.shape()[1] <<std::endl; 
      }
      return false;
    }
    return h5d_read(grp,aname,get_address(ref_.origin()),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.origin()),xfer_plist);
  }
};

}
#endif
