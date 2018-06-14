//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file TensorSoaContainer.h
 *
 * Alternative to Container<Tensor<T,D>> to support SoA algorithms
 * In addition to AoS->SoA transformation, it exploits the symmetric
 * properties of Hessian, i.e., H(x,y)==H(y,x)
 */
#ifndef QMCPLUSPLUS_TENSOR_SOA_CONTAINER_H
#define QMCPLUSPLUS_TENSOR_SOA_CONTAINER_H
namespace qmcplusplus
{
  template<typename T, unsigned D> struct TensorSoaContainer {};

  /** SoA adaptor class for ParticleAttrib<TinyVector<T,3> >
   * @tparm T data type, float, double, complex<float>, complex<double>
   */
  template<typename T>
    struct TensorSoaContainer<T,3>
    {
      typedef T           Element_t;

      ///number of elements
      int nLocal;
      ///number of elements + padded
      int nGhosts;
      ///container
      aligned_vector<T> m_data;

      ///default constructor
      TensorSoaContainer():nLocal(0),nGhosts(0){}
      ///default copy constructor
      TensorSoaContainer(const TensorSoaContainer& in)=default;
      ///default copy operator
      TensorSoaContainer& operator=(const TensorSoaContainer& in) = default;
      /** constructor with size n  without initialization
       */
      explicit TensorSoaContainer(int n)
      { resize(n); }

      /** need A=0.0;
       */
      template<typename T1>
      TensorSoaContainer& operator=(T1 in)
      {
        std::fill(m_data.begin(),m_data.end(),static_cast<T>(in));
        return *this;
      }

      ~TensorSoaContainer()=default;

      __forceinline void resize(int n)
      {
        nLocal=n;
        nGhosts=getAlignedSize<T>(n);
        m_data.resize(nGhosts*6);
      }

      /** return TinyVector<T,3>
       */
      __forceinline Tensor<T,3> operator[](int i) const
      {
        const T* restrict b=m_data+i;
        T xx=*(b);
        T xy=*(b+1*nGhosts);
        T xz=*(b+2*nGhosts);
        T yy=*(b+3*nGhosts);
        T yz=*(b+4*nGhosts);
        T zz=*(b+5*nGhosts);
        return Tensor<T,3>(xx,xy,xz,xy,yy,yz,xz,yz,zz);
      }

      ///helper class for operator ()(int i) to assign a value
      struct Accessor
      {
        int M;
        T* _base;
        Accessor()=delete;
        Accessor(const Accessor& )=delete;
        __forceinline Accessor(T* a, int ng) : _base(a), M(ng){}

        template<unsigned D>
        __forceinline Accessor& operator=(const Tensor<T,D>& rhs)
        {
          *_base      =rhs(0); *(_base+  M)=rhs(1); *(_base+2*M)=rhs(2);
          *(_base+3*M)=rhs(4); *(_base+4*M)=rhs(5); *(_base+5*M)=rhs(8);
          return *this;
        }

        /** asign value */
        template<typename T1>
        __forceinline Accessor& operator=(T1 rhs)
        {
          *_base      =rhs; *(_base+  M)=rhs; *(_base+2*M)=rhs;
          *(_base+3*M)=rhs; *(_base+4*M)=rhs; *(_base+5*M)=rhs;
          return *this;
        }
      };

      /** access operator for assignment of the i-th value
       *
       * Use for (*this)[i]=Tensor<T,3>;
       */
      __forceinline Accessor operator()(int i) 
      {
        return Accessor(m_data.data()+i,nGhosts);
      }

      ///return the base
      __forceinline T* data() { return m_data.data();}
      ///return the base
      __forceinline const T* data() const { return m_data.data();}
      ///return the base of XX components
      __forceinline T* restrict data(int i, int j) 
      { 
        const int n=(i<j)? i*3+j:j*3+i;
        return m_data().data()+n*nGhosts;
      }
      ///return the base of XX components
      __forceinline const T* restrict data(int i, int j) const
      { 
        const int n=(i<j)? i*3+j:j*3+i;
        return m_data().data()+n*nGhosts;
      }

      /** serialization function */
      template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          //ar & m_data;
          ar & nLocal & nGhosts & m_data;
        }
    };

//Incorrect: provide wrapper class
//BOOST_CLASS_TRACKING(Pos3DSoA<double,3>, boost::serialization::track_never)
//BOOST_CLASS_TRACKING(Pos3DSoA<float,3>, boost::serialization::track_never)
}

#endif

