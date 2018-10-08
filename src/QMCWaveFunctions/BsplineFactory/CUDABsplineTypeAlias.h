//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
// File derived from sections of QMCWaveFunctions/EinsplineSet.h
//////////////////////////////////////////////////////////////////////////////////////


///\file CUDABsplineTypeAlias.h
#ifndef QMCPLUSPLUS_CUDABSPLINETYPEALIAS_H
#define QMCPLUSPLUS_CUDABSPLINETYPEALIAS_H

namespace qmcplusplus
{
namespace cudasoatemp
{

////////////////////////////////////////////////////////////////////
// These small template specializations relate the CUDA types to  //
// the storage type and dimensions                                //
// that we template particular Bsplines on the host               //
////////////////////////////////////////////////////////////////////
template<typename StorageType, int dim>  struct MultiOrbitalTraits {};

template<> struct MultiOrbitalTraits<double,2>
{
  typedef multi_UBspline_2d_d SplineType;
  typedef multi_UBspline_2d_d_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<std::complex<double>,2>
{
  typedef multi_UBspline_2d_z SplineType;
  typedef multi_UBspline_2d_z_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<float,2>
{
  typedef multi_UBspline_2d_s SplineType;
  typedef multi_UBspline_2d_s_cuda CudaSplineType;
};
  
template<> struct MultiOrbitalTraits<std::complex<float>,2>
{
  typedef multi_UBspline_2d_c SplineType;
  typedef multi_UBspline_2d_c_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<double,3>
{
  typedef multi_UBspline_3d_d SplineType;
  typedef BCtype_d BCType;
  typedef double DataType;
  typedef multi_UBspline_3d_d_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<std::complex<double>,3>
{
  typedef multi_UBspline_3d_z SplineType;
  typedef BCtype_z BCType;
  typedef std::complex<double> DataType;
  typedef multi_UBspline_3d_z_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<float,3>
{
  typedef multi_UBspline_3d_s SplineType;
  typedef BCtype_s BCType;
  typedef float DataType;
  typedef multi_UBspline_3d_s_cuda CudaSplineType;
};

template<> struct MultiOrbitalTraits<std::complex<float>,3>
{
  typedef multi_UBspline_3d_c SplineType;
  typedef BCtype_c BCType;
  typedef std::complex<float> DataType;
  typedef multi_UBspline_3d_c_cuda CudaSplineType;
};

template<typename StoreType, typename CudaPrec> struct StorageTypeConverter;
template<> struct StorageTypeConverter<double,double>
{
  typedef double CudaStorageType;
};
template<> struct StorageTypeConverter<double,float>
{
  typedef float CudaStorageType;
};
template<> struct StorageTypeConverter<std::complex<double>,float>
{
  typedef std::complex<float> CudaStorageType ;
};
template<> struct StorageTypeConverter<std::complex<double>,std::complex<double> >
{
  typedef std::complex<double> CudaStorageType;
};
template<> struct StorageTypeConverter<std::complex<double>,double>
{
  typedef std::complex<double> CudaStorageType;
};

}
}
#endif //QMCPLUSPLUS_CUDABSPLINETYPEALIAS_H
