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

#ifndef QMCPLUSPLUS_HDF_MULTI_INTERFACE_H
#define QMCPLUSPLUS_HDF_MULTI_INTERFACE_H
#include <io/hdf_dataproxy.h>
#include "multi/array.hpp"
#include "multi/array_ref.hpp"

#ifdef BUILD_AFQMC 
#ifdef QMC_CUDA
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#endif
#endif

namespace qmcplusplus
{

/** specialization for vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T, class Alloc>
struct h5data_proxy<boost::multi::array<T,1,Alloc>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi::array<T,1,Alloc> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    using extensions = typename boost::multi::layout_t<1u>::extensions_type;
    if(!get_space(grp,aname,this->size(),dims))
      ref_.reextent(extensions{dims[0]});
    return h5d_read(grp,aname,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
};

template<typename T, class Alloc>
struct h5data_proxy<boost::multi::array<T,2,Alloc>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi::array<T,2,Alloc> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size(0);
    dims[1]=ref_.size(1);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.reextent({dims[0],dims[1]});
    return h5d_read(grp,aname,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
};

template<typename T, class Ptr>
struct h5data_proxy<boost::multi::array_ref<T,1,Ptr>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi::array_ref<T,1,Ptr> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0] > 0) {
        std::cerr<<" Error: multi::array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<ref_.size(0) <<std::endl;
      }
      return false;
    }
    return h5d_read(grp,aname,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
};

template<typename T, class Ptr>
struct h5data_proxy<boost::multi::array_ref<T,2,Ptr>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi::array_ref<T,2,Ptr> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size(0);
    dims[1]=ref_.size(1);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0]*dims[1] > 0) {
        std::cerr<<" Error: multi::array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<dims[1] <<" " <<ref_.size(0) <<" " <<ref_.size(1) <<std::endl;
      }
      return false;
    }
    return h5d_read(grp,aname,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(std::addressof(*ref_.origin())),xfer_plist);
  }
};


#ifdef BUILD_AFQMC 
#ifdef QMC_CUDA
// Specializations for cuda_gpu_allocator 
// Need buffered I/O and copies to gpu
template<typename T>
struct h5data_proxy<boost::multi::array<T,1,qmc_cuda::cuda_gpu_allocator<T>>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi::array<T,1,qmc_cuda::cuda_gpu_allocator<T>> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.resize({dims[0]});
    std::size_t sz = ref_.num_elements();
    boost::multi::array<T,1> buf( {sz} );
    auto ret = h5d_read(grp,aname,get_address(buf.data()),xfer_plist);
    qmc_cuda::copy_n(buf.data(),sz,ref_.origin());
    return ret;
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(to_address(ref_.origin())),xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array<T,2,qmc_cuda::cuda_gpu_allocator<T>>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi::array<T,2,qmc_cuda::cuda_gpu_allocator<T>> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size(0);
    dims[1]=ref_.size(1);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims))
      ref_.reextent({dims[0],dims[1]});
    std::size_t sz = ref_.num_elements();
    using extensions = typename boost::multi::array<T,1>::extensions_type;
    boost::multi::array<T,1> buf( extensions{sz} );
    auto ret = h5d_read(grp,aname,get_address(buf.data()),xfer_plist);
    qmc_cuda::copy_n(buf.data(),sz,ref_.origin());
    return ret;
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(to_address(ref_.origin())),xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array_ref<T,1,qmc_cuda::cuda_gpu_ptr<T>>>
    : public h5_space_type<T,1>
{
  using h5_space_type<T,1>::dims;
  using h5_space_type<T,1>::get_address;
  typedef boost::multi::array_ref<T,1,qmc_cuda::cuda_gpu_ptr<T>> data_type;
  data_type& ref_;

  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.num_elements();
  }

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0] > 0) {
        std::cerr<<" Error: multi::array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<ref_.size(0) <<std::endl;
      }
      return false;
    }
    std::size_t sz = ref_.num_elements();
    using extensions = typename boost::multi::array<T,1>::extensions_type;
    boost::multi::array<T,1> buf( extensions{sz} );
    auto ret = h5d_read(grp,aname,get_address(buf.data()),xfer_plist);
    qmc_cuda::copy_n(buf.data(),sz,ref_.origin());
    return ret;
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(to_address(ref_.origin())),xfer_plist);
  }
};

template<typename T>
struct h5data_proxy<boost::multi::array_ref<T,2,qmc_cuda::cuda_gpu_ptr<T>>>: public h5_space_type<T,2>
{
  using h5_space_type<T,2>::dims;
  using h5_space_type<T,2>::get_address;
  typedef boost::multi::array_ref<T,2,qmc_cuda::cuda_gpu_ptr<T>> data_type;
  data_type& ref_;
  inline h5data_proxy(data_type& a): ref_(a)
  {
    dims[0]=ref_.size(0);
    dims[1]=ref_.size(1);
  }
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    if(!get_space(grp,aname,this->size(),dims)) {
      if(dims[0]*dims[1] > 0) {
        std::cerr<<" Error: multi::array_ref can't be resized in h5data_proxy<>::read." <<std::endl;
        std::cerr<<dims[0] <<" " <<dims[1] <<" " <<ref_.size(0) <<" " <<ref_.size(1) <<std::endl;
      }
      return false;
    }
    std::size_t sz = ref_.num_elements();
    boost::multi::array<T,1> buf( {sz} );
    auto ret = h5d_read(grp,aname,get_address(buf.data()),xfer_plist);
    qmc_cuda::copy_n(buf.data(),sz,ref_.origin());
    return ret;
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    throw std::runtime_error(" write from gpu not implemented yet.");
    return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(to_address(ref_.origin())),xfer_plist);
  }
};


#endif
#endif

}
#endif
