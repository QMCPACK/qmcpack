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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
