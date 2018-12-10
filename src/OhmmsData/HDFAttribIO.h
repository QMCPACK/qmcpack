//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_HDF_BASE_INTERAFCE_H
#define OHMMS_HDF_BASE_INTERAFCE_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#if defined(HAVE_LIBHDF5)
#include "hdf5.h"
#else
typedef int hid_t;
typedef std::size_t hsize_t;
const int H5P_DEFAULT=0;
#endif

#include <string>


namespace qmcplusplus
{
/*\class HDFParticleAttribBase
 *\brief Base class to write/read ParticleAttrib using hdf5
 */
struct HDFAttribIOBase
{

  ///hdf node to which a dataset belongs
  hid_t my_loc;

  ///default tranfer method
  hid_t xfer_plist;

  ///default constructor
  HDFAttribIOBase():my_loc(-1),xfer_plist(H5P_DEFAULT) {}

  virtual ~HDFAttribIOBase() { }

  inline void setTransferProperty(hid_t xfer_mode)
  {
    xfer_plist=xfer_mode;
  }

  //\fn void write(fileid, name)
  //\param fileid hid_t, file id
  //\param name cosnt char*, name of the attribute
  virtual void write(hid_t, const char*) = 0;

  //\fn void write(groupid, name)
  //\param groupid hid_t, group id
  //\param name cosnt char*, name of the attribute
  virtual void read(hid_t, const char*) = 0;
};


// generic templated class for type T and is not used in reality
template<class T>
struct HDFAttribIO: public HDFAttribIOBase
{

  void write(hid_t, const char*) { }
  void read(hid_t, const char*) { }
};

}
#endif
