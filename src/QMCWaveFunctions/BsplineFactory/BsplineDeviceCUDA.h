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
class BsplineDeviceCUDA : public BsplineDevice<BsplineDeviceCUDA<ST,D>, ST, D>
{
public:
  //Code refactored from Einsplineset (still there for other functionality under AoS?)
  TinyVector<int, 3> TileFactor;
  Tensor<int, D> TileMatrix;
  
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
  //void resize_cuda(int numWalkers);
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
    qmcplusplus::cudasoatemp::create_multi_UBspline_3d_cuda((multi_UBspline_3d_d*)multi_bspline.spline_m,
							    cuda_multi_bspline);
    app_log() << "Successful copy.\n";

    L_host.resize(9);
    Linv_host.resize(9);
    Linv_cuda.resize(9);
    L_cuda.resize(9);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
    	{
    	  L_host[i*3+j]    = this->PrimLattice.R(i,j);
    	  Linv_host[i*3+j] = this->PrimLattice.G(i,j);
    	}
    L_cuda    = L_host;
    Linv_cuda = Linv_host;
  }

  void createSpline_imp(CrystalLattice<ST, D>& prim_lattice, MultiBspline<ST>& multi_spline)
  {
    std::string something_string = "something";
    qmcplusplus::cudasoatemp::create_multi_UBspline_3d_cuda(multi_spline.spline_m,
							    cuda_multi_bspline);
    // Destroy original CPU spline
    // HACK HACK HACK
    //destroy_Bspline (MultiSpline);
    L_host.resize(9);
    Linv_host.resize(9);
    Linv_cuda.resize(9);
    L_cuda.resize(9);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
    	{
    	  L_host[i*3+j]    = prim_lattice.R(i,j);
    	  Linv_host[i*3+j] = prim_lattice.G(i,j);
    	}
    L_cuda    = L_host;
    Linv_cuda = Linv_host;
  }

  //template<typename ST> void
  //template<typename ST, 3>
  void resizeStorage_imp(size_t n, size_t nvals, int num_walkers)
  {
    CudaValuePointers.resize(num_walkers);
    CudaGradLaplPointers.resize(num_walkers);
    // int N = CudaMultiSpline->num_splines;
    // CudaValueVector.resize(N*num_walkers);
    // CudaGradLaplVector.resize(4*N*num_walkers);
    // gpu::host_vector<CudaStorageType*> hostValuePointers(num_walkers);
    // gpu::host_vector<CudaStorageType*> hostGradLaplPointers(num_walkers);
    // for (int i=0; i<num_walkers; i++)
    //   {
    // 	hostValuePointers[i]    = &(CudaValueVector.data()[i*N]);
    // 	hostGradLaplPointers[i] = &(CudaGradLaplVector.data()[4*i*N]);
    //   }
    // CudaValuePointers    = hostValuePointers;
    // CudaGradLaplPointers = hostGradLaplPointers;
    // int M = MakeTwoCopies.size();
    // CudaMakeTwoCopies.resize(M);
    // CudaTwoCopiesIndex.resize(M);
    // gpu::host_vector<int> hostMakeTwoCopies(M);
    // gpu::host_vector<int> hostTwoCopiesIndex(M);
    // int TwoCopiesIndexCounter = 0;
    // for (int i=0; i<M; i++)
    //   {
    // 	hostMakeTwoCopies[i] = MakeTwoCopies[i];
    // 	hostTwoCopiesIndex[i] = TwoCopiesIndexCounter;
    // 	TwoCopiesIndexCounter = MakeTwoCopies[i] ? TwoCopiesIndexCounter+2 : TwoCopiesIndexCounter+1;
    //   }
    // CudaMakeTwoCopies = hostMakeTwoCopies;
    // CudaTwoCopiesIndex = hostTwoCopiesIndex;
    // CudakPoints.resize(M);
    // CudakPoints_reduced.resize(M);
    // gpu::host_vector<TinyVector<CUDA_PRECISION,OHMMS_DIM> > hostkPoints(M),
    //   hostkPoints_reduced(M);
    // for (int i=0; i<M; i++)
    //   {
    // 	//      PosType k_red1 = PrimLattice.toCart(kPoints[i]);
    // 	PosType k_red2(dot(kPoints[i], PrimLattice.a(0)),
    // 		       dot(kPoints[i], PrimLattice.a(1)),
    // 		       dot(kPoints[i], PrimLattice.a(2)));
    // 	//       fprintf (stderr, "kred1 = %8.3f %8.3f %8.3f\n", k_red1[0], k_red1[1], k_red1[2]);
    // 	//       fprintf (stderr, "kred2 = %8.3f %8.3f %8.3f\n", k_red2[0], k_red2[1], k_red2[2]);
    // 	for (int j=0; j<OHMMS_DIM; j++)
    // 	  {
    // 	    hostkPoints[i][j]         = kPoints[i][j];
    // 	    hostkPoints_reduced[i][j] = k_red2[j];
    // 	  }
    //   }
    // CudakPoints = hostkPoints;
    // CudakPoints_reduced = hostkPoints_reduced;
  }

  void implementation()
  {
    std::cout<< "implemented for CUDA\n";
  }
};

}
#endif
