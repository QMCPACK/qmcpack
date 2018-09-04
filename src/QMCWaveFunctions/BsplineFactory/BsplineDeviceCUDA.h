//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


///\file BsplineDeviceCUDA.h
#ifndef QMCPLUSPLUS_BSPLINEDEVICECUDA_H
#define QMCPLUSPLUS_BSPLINEDEVICECUDA_H

#include <iostream>
#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineDevice.h"
#include "QMCWaveFunctions/BsplineFactory/CUDABsplineConversions.h"
//#include "einspline/bspline_structs.h"
#include "QMCWaveFunctions/BsplineFactory/CUDABsplineTypeAlias.h"
#include "CUDA/gpu_vector.h"

namespace qmcplusplus
{

/** 
 * \class BsplineDeviceCUDA
 * 
 * Attempting to move EinsplineSet Cuda behavior here
 */
  
template<typename ST, unsigned D>
class BsplineDeviceCUDA : BsplineDevice<BsplineDeviceCUDA<ST,D>, ST, D>
{
public:
  //Code refactored from Einsplineset (still there for other functionality under AoS?)

  // I prefer this as a method of bringing in QMCTraits typedefs
  using QMCT = QMCTraits;
  using CudaRealType = QMCT::CudaRealType;
  using CudaPosType = QMCT::CudaPosType;
  //Pending BsplineReaderBase that doesn't pull in EinsplineSet.h
  using CudaStorageType = typename cudasoatemp::StorageTypeConverter<ST,CUDA_PRECISION>::CudaStorageType;
  using CudaSplineType =  typename cudasoatemp::MultiOrbitalTraits<CudaStorageType,OHMMS_DIM>::CudaSplineType;

  CudaSplineType* cuda_multi_bspline;
  gpu::device_vector<CudaStorageType> CudaValueVector, CudaGradLaplVector;
  gpu::device_vector<CudaStorageType*> CudaValuePointers, CudaGradLaplPointers;
  void resize_cuda(int numWalkers);
  // Cuda equivalent
  gpu::device_vector<int> CudaMakeTwoCopies;
  gpu::device_vector<int> CudaTwoCopiesIndex;
  // Cuda equivalent
  gpu::device_vector<TinyVector<CUDA_PRECISION,OHMMS_DIM > > CudakPoints,
      CudakPoints_reduced;
  void applyPhaseFactors (gpu::device_vector<CudaStorageType*> &storageVector,
                          gpu::device_vector<CudaRealType*> &phi);
  // Data for vectorized evaluations
  std::vector<CudaPosType> hostPos;
  gpu::host_vector<CudaPosType> NLhostPos;
  gpu::device_vector<CudaPosType> cudapos, NLcudapos;
  gpu::host_vector<CudaRealType> hostSign, NLhostSign;
  gpu::device_vector<CudaRealType> cudaSign, NLcudaSign;
  // This stores the inverse of the lattice vector matrix in
  // GPU memory.
  gpu::device_vector<CudaRealType> Linv_cuda, L_cuda;
  gpu::host_vector<CudaRealType> L_host, Linv_host;

public:
  
  void initDevice_imp(MultiBspline<ST>& multi_bspline)
  {
    app_log() << "Copying einspline orbitals to GPU.\n";
    create_multi_UBspline_3d_cuda(multi_bspline, cuda_multi_bspline);
    app_log() << "Successful copy.\n";
    // L_host.resize(9);
    // Linv_host.resize(9);
    // Linv_cuda.resize(9);
    // L_cuda.resize(9);
    // for (int i=0; i<3; i++)
    //   for (int j=0; j<3; j++)
    // 	{
    // 	  L_host[i*3+j]    = PrimLattice.R(i,j);
    // 	  Linv_host[i*3+j] = PrimLattice.G(i,j);
    // 	}
    // L_cuda    = L_host;
    // Linv_cuda = Linv_host;
  }
  
  template<typename GT, typename BCT>
  void createSpline(MultiBspline<ST>& multi_spline)
  {
    cudasoatemp::create_multi_UBspline_3d_cuda(multi_spline, cuda_multi_bspline);
    // Destroy original CPU spline
    // HACK HACK HACK
    //destroy_Bspline (MultiSpline);
    // L_host.resize(9);
    // Linv_host.resize(9);
    // Linv_cuda.resize(9);
    // L_cuda.resize(9);
    // for (int i=0; i<3; i++)
    //   for (int j=0; j<3; j++)
    // 	{
    // 	  L_host[i*3+j]    = PrimLattice.R(i,j);
    // 	  Linv_host[i*3+j] = PrimLattice.G(i,j);
    // 	}
    // L_cuda    = L_host;
    // Linv_cuda = Linv_host;
  }
  
  void implementation()
  {
    std::cout<< "implemented for CUDA\n";
  }
};

}
#endif
