//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONTAINER_PROXY_MULTI_H
#define QMCPLUSPLUS_CONTAINER_PROXY_MULTI_H

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

namespace qmcplusplus
{

template<typename T, class Alloc>
struct container_proxy<boost::multi::array<T,2,Alloc> >
{
  enum {DIM=scalar_traits<T>::DIM};
  typedef typename container_proxy<T>::pointer pointer;
  boost::multi::array<T,2,Alloc>& ref;
  inline container_proxy(boost::multi::array<T,2,Alloc>& a):ref(a) {}
  inline size_t size() const
  {
    return ref.num_elements()*DIM;
  }
  inline pointer data()
  {
    //using detail::to_address;
    //return scalar_traits<T>::get_address(to_address(ref.origin()));
    return scalar_traits<T>::get_address(std::addressof(*ref.origin()));
  }
  inline void resize(size_t n)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array<T,D,Alloc> >. \n");
  }
  template<typename I>
  inline void resize(I* n, int d)
  {
    if(d < 2)
      APP_ABORT(" Error: Inconsistent dimension in container_proxy<boost::multi::array<T,D,Alloc> >::resize(I*,int). \n");
    ref.reextent({n[0],n[1]});
  }
};

template<typename T>
struct container_proxy<boost::multi::array_ref<T,2> >
{
  enum {DIM=scalar_traits<T>::DIM};
  typedef typename container_proxy<T>::pointer pointer;
  boost::multi::array_ref<T,2>& ref;
  inline container_proxy(boost::multi::array_ref<T,2>& a):ref(a) {}
  inline size_t size() const
  {
    return ref.num_elements()*DIM;
  }
  inline pointer data()
  {
    //using detail::to_address;
    //return scalar_traits<T>::get_address(to_address(ref.origin()));
    return scalar_traits<T>::get_address(std::addressof(*ref.origin()));
  }
  inline void resize(size_t n)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array_ref<T,D> >. \n");
  }

  template<typename I>
  inline void resize(I* n, int d)
  {
    APP_ABORT(" Error: Can not resize container_proxy<boost::multi::array_ref<T,D> >. \n");
  }

};

}

#endif
