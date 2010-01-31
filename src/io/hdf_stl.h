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
#ifndef QMCPLUSPLUS_HDFSTRINGATTRIB_H
#define QMCPLUSPLUS_HDFSTRINGATTRIB_H
#include <vector>
#include <sstream>        
#include <bitset>
namespace qmcplusplus
{
  /** specialization for vector<T>
   *
   * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
   */ 
  template<typename T> struct HDFAttribIO<vector<T> >
    : public h5_space_type<T,1>
  { 
    using h5_space_type<T,1>::dims;
    using h5_space_type<T,1>::get_address;
    typedef vector<T> data_type;
    data_type& ref_;

    inline HDFAttribIO(data_type& a): ref_(a) { dims[0]=ref_.size(); }

    inline void read(hid_t grp, const std::string& aname)
    {
      if(!h5d_getspace(grp,aname,this->size(),dims)) ref_.resize(dims[0]);
      h5d_read(grp,aname,get_address(&ref_[0]));
    }

    inline void write(hid_t grp, const std::string& aname)
    {
      h5d_write(grp,aname.c_str(),this->size(),dims,get_address(&ref_[0]));
    }
  };

  template<std::size_t N>
    struct HDFAttribIO<std::bitset<N> >
    {
      //NOTE: This specialization assumes each complex number was/will be saved
      // as a pair of doubles. This is checked.
      typedef std::bitset<N> ArrayType_t;
      ArrayType_t& ref;

      //Assumes complex stored as a pair of floats/doubles.
      HDFAttribIO<ArrayType_t>(ArrayType_t& a): ref(a)
      { 
      }

      inline void write(hid_t grp, const char* name) {
        unsigned long c=ref.to_ulong();
        HDFAttribIO<unsigned long> hc(c);
        hc.write(grp,name);
      }

      inline void read(hid_t grp, const char* name) {
        unsigned long c=ref.to_ulong();
        HDFAttribIO<unsigned long> hc(c);
        hc.read(grp,name);
        ref=c;
      }
    };


  /** Specialization for string */
  template<> struct HDFAttribIO<std::string>
    {

      typedef std::string ArrayType_t;
      ArrayType_t& ref;

      HDFAttribIO<ArrayType_t>(ArrayType_t& a): ref(a) { }

      inline void write(hid_t grp, const char* name) 
      {
        hid_t str80 = H5Tcopy(H5T_C_S1);
        H5Tset_size(str80,ref.size());
        hsize_t dim = 1;
        hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
        hid_t dataset =  H5Dcreate(grp, name, str80, dataspace, H5P_DEFAULT);
        hid_t ret = H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
        H5Sclose(dataspace);
        H5Dclose(dataset);
      }

      inline void read(hid_t grp, const char* name) 
      {
	// Turn off error printing
	H5E_auto_t func;
	void *client_data;
	H5Eget_auto (&func, &client_data);
	H5Eset_auto (NULL, NULL);
	
        hid_t dataset = H5Dopen(grp,name);
	if (dataset > -1) {
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
	  hid_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
	  // Erase trailing null character
	  ref.erase (dim_out-1, 1);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	}
	// Turn error printing back on
	H5Eset_auto (func, client_data);
      }
    };

  template<>
    struct HDFAttribIO<std::ostringstream>//: public HDFAttribIOBase 
    {
      typedef std::ostringstream Data_t;
      Data_t& ref;

      HDFAttribIO<Data_t>(Data_t& a): ref(a) { }

      inline void write(hid_t grp, const char* name) {
        herr_t status = H5Eset_auto(NULL, NULL);
        status = H5Gget_objinfo (grp, name, 0, NULL);
        hsize_t str80 = H5Tcopy(H5T_C_S1);
        H5Tset_size(str80,ref.str().size());
        if(status ==0)
        {
          hid_t dataset = H5Dopen(grp, name);
          hid_t ret = H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.str().c_str());
          H5Dclose(dataset);
        }
        else
        {
          hsize_t dim = 1;
          hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
          hid_t dataset =  H5Dcreate(grp, name, str80, dataspace, H5P_DEFAULT);
          hid_t ret = H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.str().c_str());
          H5Sclose(dataspace);
          H5Dclose(dataset);
        }
      }

      inline void read(hid_t grp, const char* name) {
      }
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: BoostRandom.h 894 2006-02-03 16:52:38Z jnkim $ 
 ***************************************************************************/
