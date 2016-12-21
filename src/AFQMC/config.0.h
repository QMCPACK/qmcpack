#ifndef AFQMC_CONFIG_0_H 
#define AFQMC_CONFIG_0_H 

#include <string>
#include <algorithm>
#include<cstdlib>
#include<ctype.h>
#include <vector>
#include <map>
#include <complex>
#include <tuple>
#include <fstream>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>


//#include"../FCIQMC/config.0.h"

namespace qmcplusplus
{

//  template<typename spT> using ShmemAllocator = boost::interprocess::allocator<spT, boost::interprocess::managed_shared_memory::segment_manager>;
//  template<typename spT> using SMVector = boost::interprocess::vector<spT, ShmemAllocator<spT>>;

  typedef uint32_t                   IndexType;  
  typedef uint32_t                    OrbitalType;  
//  typedef uint16_t                    OrbitalType; 
  typedef double                         RealType;
#if defined(QMC_COMPLEX)
  typedef std::complex<RealType>  ValueType;
#else
  typedef RealType                       ValueType;
#endif
  typedef std::complex<RealType>         ComplexType;


#define HAVE_CPP11
#ifdef HAVE_CPP11

  template<typename T> using s1D = std::tuple<IndexType,T>;
  template<typename T> using s2D = std::tuple<IndexType,IndexType,T>;
  template<typename T> using s3D = std::tuple<IndexType,IndexType,IndexType,T>;
  template<typename T> using s4D = std::tuple<IndexType,IndexType,IndexType,IndexType,T>;

#else

  template<typename T> struct s1D { T value; IndexType i1; };
  template<typename T> struct s2D { T value; IndexType i1, i2; };
  template<typename T> struct s3D { T value; IndexType i1, i2, i3; };
  template<typename T> struct s4D { T value; IndexType i1, i2, i3, i4; };
  template<typename T> struct s5D { T value; IndexType i1, i2, i3, i4, i5; };
  template<typename T> struct s9D { T value; IndexType i1, i2, i3, i4, i5, i6, i7, i8, i9; };

#endif

template<typename T>
inline bool isComplex(const T& a) 
{
  return std::is_same<T,std::complex<RealType>>::value; 
}

template<typename T>
inline ComplexType toComplex(const T& a);

template<>
inline ComplexType toComplex(const RealType& a) 
{
  return ComplexType(a,0.0); 
}

template<>
inline ComplexType toComplex(const std::complex<RealType>& a) 
{
  return a; 
}

template<typename T>
inline void setImag(T& a, RealType b);

template<>
inline void setImag(RealType& a, RealType b)
{
}

template<>
inline void setImag(std::complex<RealType>& a, RealType b)
{
  a.imag(b);
}

template<typename T>
inline T myconj(const T& a) 
{
  return a;
}

template<typename T>
inline std::complex<T> myconj(const std::complex<T>& a) 
{
  return std::conj(a);
}

template<typename T>
inline RealType mynorm(const T& a)
{
  return a*a;
}

template<typename T>
inline RealType mynorm(const std::complex<T> &a)
{
  return std::norm(a);
}

template<typename T>
inline std::complex<T> operator*(const int &lhs, const std::complex<T> &rhs)
{
  return T(lhs) * rhs;
}

template<typename T>
inline std::complex<T> operator*(const std::complex<T> &lhs, const int &rhs)
{
  return lhs * T(rhs);
}

inline bool sortDecreasing (int i,int j) { return (i>j); }

struct _mySort_snD_ {
  bool operator() (const s1D<RealType>& lhs, const s1D<RealType>& rhs)
  { return (bool)(std::get<0>(lhs) < std::get<0>(rhs)); 
  } 
  bool operator() (const s2D<RealType>& lhs, const s2D<RealType>& rhs)
  { return (bool)(std::get<0>(lhs) < std::get<0>(rhs)) || 
           ( !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) && 
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) );  
  } 
  bool operator() (const s4D<RealType>& lhs, const s4D<RealType>& rhs)
  {
    return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
  }
/* I'm having issues with this stupid function. What's wrong???
  { return    (bool)(std::get<0>(lhs) < std::get<0>(rhs)) || 
         (   !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) && 
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) ||
          (  !(bool)(std::get<1>(rhs) < std::get<1>(lhs)) &&
              (bool)(std::get<2>(lhs) < std::get<2>(rhs)) ||
           ( !(bool)(std::get<2>(rhs) < std::get<2>(lhs)) &&
              (bool)(std::get<3>(lhs) < std::get<3>(rhs))))); 
  }
*/
  bool operator() (const s1D<std::complex<RealType> >& lhs, const s1D<std::complex<RealType> >& rhs)
  { return (bool)(std::get<0>(lhs) < std::get<0>(rhs)); 
  } 
  bool operator() (const s2D<std::complex<RealType> >& lhs, const s2D<std::complex<RealType> >& rhs)
  { return (bool)(std::get<0>(lhs) < std::get<0>(rhs)) || 
           ( !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) && 
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) );  
  } 
  bool operator() (const s4D<std::complex<RealType> >& lhs, const s4D<std::complex<RealType> >& rhs)
  {
    return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
  }
/*
  { return    (bool)(std::get<0>(lhs) < std::get<0>(rhs)) || 
         (   !(bool)(std::get<0>(rhs) < std::get<0>(lhs)) && 
              (bool)(std::get<1>(lhs) < std::get<1>(rhs)) ||
          (  !(bool)(std::get<1>(rhs) < std::get<1>(lhs)) &&
              (bool)(std::get<2>(lhs) < std::get<2>(rhs)) ||
           ( !(bool)(std::get<2>(rhs) < std::get<2>(lhs)) &&
              (bool)(std::get<3>(lhs) < std::get<3>(rhs))))); 
  }
*/
};

struct _myEqv_snD_ {
  // equivalence
  bool operator() (const s1D<RealType>& lhs, const s1D<RealType>& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs));
  }
  bool operator() (const s2D<RealType>& lhs, const s2D<RealType>& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs))
        && (bool)(std::get<1>(lhs) == std::get<1>(rhs));
  }
  bool operator() (const s4D<RealType>& lhs, const s4D<RealType>& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs))
        && (bool)(std::get<1>(lhs) == std::get<1>(rhs))
        && (bool)(std::get<2>(lhs) == std::get<2>(rhs))
        && (bool)(std::get<3>(lhs) == std::get<3>(rhs));
  }
  bool operator() (const s1D<std::complex<RealType> >& lhs, const s1D<std::complex<RealType> >& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs));
  }
  bool operator() (const s2D<std::complex<RealType> >& lhs, const s2D<std::complex<RealType> >& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs))
        && (bool)(std::get<1>(lhs) == std::get<1>(rhs));
  }
  bool operator() (const s4D<std::complex<RealType> >& lhs, const s4D<std::complex<RealType> >& rhs)
  { return (bool)(std::get<0>(lhs) == std::get<0>(rhs))
        && (bool)(std::get<1>(lhs) == std::get<1>(rhs))
        && (bool)(std::get<2>(lhs) == std::get<2>(rhs))
        && (bool)(std::get<3>(lhs) == std::get<3>(rhs));
  }
}; 


}


namespace std {
template<typename T>
inline bool operator<(const std::complex<T> &lhs, const std::complex<T> &rhs)
{
  if (lhs.real() != rhs.real())
  {
    return lhs.real() < rhs.real();
  }
  return lhs.imag() < rhs.imag();
}
}

#endif
