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
#ifndef QMCPLUSPLUS_HDF_BOOST_SMVector_INTERFACE_H
#define QMCPLUSPLUS_HDF_BOOST_SMVector_INTERFACE_H
#include <vector>
#include <sstream>
#include <bitset>
#include <deque>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
namespace qmcplusplus
{
/** specialization for vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T> struct h5data_proxy<boost::interprocess::vector<T, boost::interprocess::allocator<T, boost::interprocess::managed_shared_memory::segment_manager> > >
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::interprocess::vector<T, boost::interprocess::allocator<T, boost::interprocess::managed_shared_memory::segment_manager> > data_type;
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

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: hdf_stl.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/
