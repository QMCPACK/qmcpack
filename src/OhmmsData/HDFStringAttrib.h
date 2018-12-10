//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_HDFSTRINGATTRIB_H
#define QMCPLUSPLUS_HDFSTRINGATTRIB_H
#include "OhmmsData/HDFAttribIO.h"
#include <sstream>
namespace qmcplusplus
{
/** Specialization for std::string */
template<>
struct HDFAttribIO<std::string>: public HDFAttribIOBase
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
struct HDFAttribIO<std::ostringstream>: public HDFAttribIOBase
{
  typedef std::ostringstream Data_t;
  Data_t& ref;

  HDFAttribIO<Data_t>(Data_t& a): ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
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

  inline void read(hid_t grp, const char* name)
  {
  }
};
}
#endif
