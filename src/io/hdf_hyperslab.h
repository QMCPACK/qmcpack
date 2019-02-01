/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_HYPERSLAB_IO_H
#define QMCPLUSPLUS_HDF_HYPERSLAB_IO_H
#include <type_traits/container_proxy.h>
#include <io/hdf_datatype.h>
#include <io/hdf_dataspace.h>
#include <io/hdf_dataproxy.h>
namespace qmcplusplus
{
/** class to use hyperslab with a serialized container
 *
 * container_proxy<CT> handles the size and datatype
 */
template<typename CT, unsigned MAXDIM>
struct hyperslab_proxy: public container_proxy<CT>
{
  ///determine the size of value_type
  enum {element_size=container_proxy<CT>::DIM};
  ///rank of hyperslab
  int slab_rank;
  ///true, if hyperslab is used
  bool use_slab;
  ///global dimension of the hyperslab
  TinyVector<hsize_t,MAXDIM+1> slab_dims;
  ///local dimension of the hyperslab
  TinyVector<hsize_t,MAXDIM+1> slab_dims_local;
  ///offset of the hyperslab
  TinyVector<hsize_t,MAXDIM+1> slab_offset;
  ///1D
  hyperslab_proxy(CT& a): container_proxy<CT>(a), slab_rank(a.slab_rank),
    slab_dims(a.slab_dims), slab_offset(a.slab_offset)
  {
    slab_dims_local=slab_dims;
    use_slab=false;
  }

  template<typename IC>
  inline hyperslab_proxy(CT& a, const IC& dims_in):container_proxy<CT>(a)
  {
    use_slab=false;
    slab_rank=dims_in.size();
    for(int i=0; i<dims_in.size(); ++i)
      slab_dims[i]=static_cast<hsize_t>(dims_in[i]);
    slab_dims_local=slab_dims;
    if(element_size>1)
    {
      slab_dims[slab_rank]=element_size;
      slab_rank+=1;
    }
  }

  template<typename IC1, typename IC2, typename IC3>
  inline hyperslab_proxy(CT& a, const IC1& dims_in, const IC2& dims_loc, const IC3& offsets_in)
  :container_proxy<CT>(a)
  {
    slab_rank=dims_in.size();
    for(int i=0; i<dims_in.size(); ++i)
      slab_dims[i]=static_cast<hsize_t>(dims_in[i]);
    for(int i=0; i<dims_loc.size(); ++i)
      slab_dims_local[i]=static_cast<hsize_t>(dims_loc[i]);
    for(int i=0; i<dims_in.size(); ++i)
      slab_offset[i]=static_cast<hsize_t>(offsets_in[i]);
    if(element_size>1)
    {
      slab_dims[slab_rank]=element_size;
      slab_dims_local[slab_rank]=element_size;
      slab_offset[slab_rank]=0;
      slab_rank+=1;
    }
    use_slab=true;
  }

  /** return the size of the i-th dimension
   * @param i dimension
   */
  inline hsize_t size(int i) const
  {
    return (i>MAXDIM)?0:slab_dims[i];
  }

  inline void change_shape()
  {
    this->resize(slab_dims.data(),MAXDIM+1);
  }

};

template<typename CT, unsigned MAXDIM>
struct h5data_proxy<hyperslab_proxy<CT,MAXDIM> >
{
  hyperslab_proxy<CT,MAXDIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT,MAXDIM>& a): ref_(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(ref_.use_slab)
    {
      return h5d_read(grp,aname.c_str(),
          ref_.slab_rank,
          ref_.slab_dims.data(),
          ref_.slab_dims_local.data(),
          ref_.slab_offset.data(),
          ref_.data(),xfer_plist);
    }
    else
    {
      int rank=ref_.slab_rank;
      if(!get_space(grp,aname,rank,ref_.slab_dims.data(),true))
      {
        ref_.change_shape();
      }
      return h5d_read(grp,aname,ref_.data(),xfer_plist);
    }
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(ref_.use_slab)
    {
      return h5d_write(grp,aname.c_str(),
          ref_.slab_rank,
          ref_.slab_dims.data(),
          ref_.slab_dims_local.data(),
          ref_.slab_offset.data(),
          ref_.data(),xfer_plist);
    }
    else{
    return h5d_write(grp,aname.c_str(),ref_.slab_rank,ref_.slab_dims.data(),ref_.data(),xfer_plist);
    }
  }
};

#if 0
template<typename T, unsigned MAXDIM>
struct h5data_proxy<hyperslab_proxy<boost::multi::array<T,2,cuda::cuda_gpu_allocator<T>>,MAXDIM> >
{
  typedef boost::multi::array<T,2,cuda::cuda_gpu_allocator<T>> CT;
  hyperslab_proxy<CT,MAXDIM>& ref_;
  h5data_proxy(hyperslab_proxy<CT,MAXDIM>& a): ref_(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(ref_.use_slab)
    {
// later on specialize h5d_read for fancy pointers 
      std::size_t sz = ref_.ref.num_elements();
      boost::multi::array<T,1> buf( {sz} );
      auto ret = h5d_read(grp,aname.c_str(),
          ref_.slab_rank,
          ref_.slab_dims.data(),
          ref_.slab_dims_local.data(),
          ref_.slab_offset.data(),
          reinterpret_cast<T*>(buf.data()),xfer_plist);
      cuda::copy_n(buf.data(),sz,ref_.ref.origin());
      return ret;
    }
    else
    {
      int rank=ref_.slab_rank;
      if(!get_space(grp,aname,rank,ref_.slab_dims.data(),true))
      {
        std::cerr<<" Disabled hyperslab resize with boost::multi::array_ref.\n";
        return false;
      }
      return h5d_read(grp,aname,ref_.data(),xfer_plist);
    }
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    std::cerr<<" Disabled hyperslab write with boost::multi::array_ref.\n";
    return false;
  }
};
#endif
}
#endif
