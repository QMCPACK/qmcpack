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
#ifndef QMCPLUSPLUS_HDF_PETE_TRAITS_H
#define QMCPLUSPLUS_HDF_PETE_TRAITS_H

#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Utilities/PooledData.h>

namespace qmcplusplus {

  /** specialization for TinyVector<T,D>
   */
  template<typename T, unsigned D> struct HDFAttribIO<TinyVector<T,D> >
    : public h5_space_type<T,1>
  {
    using h5_space_type<T,1>::dims;
    using h5_space_type<T,1>::get_address;
    typedef TinyVector<T,D> data_type;
    data_type& ref_;
    inline HDFAttribIO(data_type& a): ref_(a) 
    {
      dims[0]=D;
    }
    inline void read(hid_t grp, const std::string& aname)
    {
      //probably want to check the dimension
      //if(get_dataspace(grp,aname))
        h5d_read(grp,aname,get_address(ref_.data()));
    }
    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()));
    }
  };

  /** specialization for Tensor<T,D>
   */
  template<typename T, unsigned D> struct HDFAttribIO<Tensor<T,D> >
    : public h5_space_type<T,2>
  {
    using h5_space_type<T,2>::dims;
    using h5_space_type<T,2>::get_address;
    typedef Tensor<T,D> data_type;
    data_type& ref_;
    inline HDFAttribIO(data_type& a): ref_(a) 
    {
      dims[0]=D; dims[1]=D;
    }
    inline void read(hid_t grp, const std::string& aname)
    {
      h5d_read(grp,aname,get_address(ref_.data()));
    }
    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()));
    }
  };


  /** specialization for Vector<T>
   *
   * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
   */ 
  template<typename T> struct HDFAttribIO<Vector<T> >: public h5_space_type<T,1>
  { 
    using h5_space_type<T,1>::dims;
    using h5_space_type<T,1>::get_address;
    typedef Vector<T> data_type;
    data_type& ref_;

    inline HDFAttribIO(data_type& a): ref_(a) { dims[0]=ref_.size(); }

    inline void read(hid_t grp, const std::string& aname)
    {
      if(!h5d_getspace(grp,aname,this->size(),dims)) ref_.resize(dims[0]);
      h5d_read(grp,aname,get_address(ref_.data()));
    }

    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()));
    }
  };


  template<typename T>
    struct HDFAttribIO<Matrix<T> >: public h5_space_type<T,2>
  {
    using h5_space_type<T,2>::dims;
    using h5_space_type<T,2>::get_address;
    typedef Matrix<T> data_type;
    data_type& ref_;
    inline HDFAttribIO(data_type& a): ref_(a) 
    {
      dims[0]=ref_.rows();
      dims[1]=ref_.cols();
    }
    inline void read(hid_t grp, const std::string& aname)
    {
      if(!h5d_getspace(grp,aname,this->size(),dims)) ref_.resize(dims[0],dims[1]);
      h5d_read(grp,aname,get_address(ref_.data()));
    }
    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()));
    }
  };


  template<typename T, unsigned D> 
    struct HDFAttribIO<Array<T,D> >: public h5_space_type<T,D>
  {
    using h5_space_type<T,D>::dims;
    using h5_space_type<T,D>::get_address;
    typedef Array<T,D> data_type;
    data_type& ref_;
    inline HDFAttribIO(data_type& a): ref_(a) 
    {
      for(int i=0; i<D; ++i) dims[i]=ref_.size(i);
    }
    inline void read(hid_t grp, const std::string& aname)
    {
      if(!h5d_getspace(grp,aname,this->size(),dims)) ref_.resize(dims);
      h5d_read(grp,aname,get_address(ref_.data()));
    }
    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.data()));
    }
  };

  /** specialization for Vector<T>
   *
   * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
   */ 
  template<typename T> struct HDFAttribIO<PooledData<T> >
  { 
    vector<hsize_t> dims;
    typedef PooledData<T> data_type;
    data_type& ref_;

    inline HDFAttribIO(data_type& a): ref_(a){dims.resize(1,a.size());}
    inline HDFAttribIO(data_type& a, const vector<hsize_t>& dims_in): ref_(a),dims(dims_in){}

    inline void read(hid_t grp, const std::string& aname)
    {
      const int rank=dims.size();
      if(!h5d_getspace(grp,aname,rank,&dims[0])) 
      {
        size_t ntot=dims[0]; 
        for(int i=1;i<rank;++i) ntot*=dims[i];
        ref_.resize(ntot);
      }
      h5d_read(grp,aname,ref_.data());
    }

    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),dims.size(),&dims[0],ref_.data());
    }
  };


}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2764 $   $Date: 2008-06-26 10:21:31 -0500 (Thu, 26 Jun 2008) $
 * $Id: HDFNumericAttrib.h 2764 2008-06-26 15:21:31Z jnkim $ 
 ***************************************************************************/
