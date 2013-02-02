//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file einspline_util.hpp
 * @brief utility functions for i/o and bcast of einspline objects
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_UTILITIES_H
#define QMCPLUSPLUS_EINSPLINE_UTILITIES_H

#include <Message/CommOperators.h>
#include <OhmmsData/FileUtility.h>
#include <io/hdf_archive.h>

namespace qmcplusplus
{
  ///handles i/o and bcast, testing for now
  template<typename T>
  inline void chunked_bcast(Communicate* comm, T* buffer, size_t ntot)
  {
    if(comm->size()==1) return;

    size_t chunk_size=(1<<30)/sizeof(T); //256 MB
    int n=static_cast<int>(ntot/chunk_size); 

    size_t offset=0;
    for(int i=0; i<n; ++i, offset+=chunk_size)
    {
      comm->bcast(buffer+offset,static_cast<int>(chunk_size));
    }

    if(offset<ntot)
    {
      comm->bcast(buffer+offset,static_cast<int>(ntot-offset));
    }
  }

  template<typename ENGT>
  inline void chunked_bcast(Communicate* comm, ENGT* buffer)
  {
    chunked_bcast(comm,buffer->coefs, buffer->coefs_size);
  }


  /** specialization of h5data_proxy for einspline_engine
   */
  template<typename ENGT>
    struct h5data_proxy<einspline_engine<ENGT> >
    : public h5_space_type<typename einspline_engine<ENGT>::value_type,4>
  {
    typedef typename einspline_engine<ENGT>::value_type value_type;
    using h5_space_type<value_type,4>::dims;
    using h5_space_type<value_type,4>::get_address;
    typedef einspline_engine<ENGT> data_type;

    data_type& ref_;

    inline h5data_proxy(data_type& a): ref_(a)
    {
      dims[0]=a.spliner->x_grid.num+3;
      dims[1]=a.spliner->y_grid.num+3;
      dims[2]=a.spliner->z_grid.num+3;
      dims[3]=a.spliner->z_stride;
    }

    inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      if(ref_.spliner) 
        return h5d_read(grp,aname,get_address(ref_.spliner->coefs),xfer_plist);
      else
        return false;
    }

    inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.spliner->coefs),xfer_plist);
    }
  };

  string make_spline_filename(const string& old, int twist, const TinyVector<int,3>& mesh)
  {
    string aname(old);
    if(getExtension(aname) == "h5")
    {
      aname.erase(aname.end()-3,aname.end());
    }
    ostringstream oo;
    oo<<".tw" << twist <<".g"<<mesh[0]<<"x"<<mesh[1]<<"x"<<mesh[2]<<".h5";
    aname+=oo.str();
    return aname;
  }

}
#endif
