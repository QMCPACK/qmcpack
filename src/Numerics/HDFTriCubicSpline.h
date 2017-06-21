//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_HDF_TRICUBICSPLINE_H
#define OHMMS_HDF_TRICUBICSPLINE_H

#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/TriCubicSplineT.h"

/** Specialization for std::vector<double> */
template<>
struct HDFAttribIO<TriCubicSplineT<double> >: public HDFAttribIOBase
{

  typedef TriCubicSplineT<double> Data_t;
  Data_t&  ref;

  HDFAttribIO<Data_t>(Data_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim[4];
    dim[0]=ref.nX;
    dim[1]=ref.nY;
    dim[2]=ref.nZ;
    dim[3]=8;
    hid_t dataspace  = H5Screate_simple(4, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hsize_t dim[4];
    dim[0]=ref.nX;
    dim[1]=ref.nY;
    dim[2]=ref.nZ;
    dim[3]=8;
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};
#endif
