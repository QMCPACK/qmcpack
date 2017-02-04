//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Christos Kartsaklis, kartsaklisc@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/EinsplineSet.h"
#include <einspline/multi_bspline.h>
#include <einspline/multi_bspline_eval_cuda.h>
#include "Configuration.h"
#include "AtomicOrbitalCuda.h"
#include "PhaseFactors.h"
#ifdef HAVE_MKL
#include <mkl_vml.h>
#endif

namespace qmcplusplus
{
inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_s_cuda* &out)
{
  out = create_multi_UBspline_3d_s_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_d *in,
    multi_UBspline_3d_d_cuda * &out)
{
  out = create_multi_UBspline_3d_d_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_c_cuda* &out)
{
  out = create_multi_UBspline_3d_c_cuda_conv (in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_z_cuda * &out)
{
  out = create_multi_UBspline_3d_z_cuda(in);
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_d_cuda * &out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}

inline void create_multi_UBspline_3d_cuda (multi_UBspline_3d_z *in,
    multi_UBspline_3d_s_cuda * &out)
{
  app_error() << "Attempted to convert complex CPU spline into a real "
              << " GPU spline.\n";
  abort();
}

// Real evaluation functions
inline void
EinsplineMultiEval (multi_UBspline_3d_d *restrict spline,
                    TinyVector<double,3> r,
                    Vector<double> &psi)
{
  eval_multi_UBspline_3d_d (spline, r[0], r[1], r[2], psi.data());
}

inline void
EinsplineMultiEval (multi_UBspline_3d_d *restrict spline,
                    TinyVector<double,3> r,
                    std::vector<double> &psi)
{
  eval_multi_UBspline_3d_d (spline, r[0], r[1], r[2], &(psi[0]));
}


inline void
EinsplineMultiEval (multi_UBspline_3d_d *restrict spline,
                    TinyVector<double,3> r,
                    Vector<double> &psi,
                    Vector<TinyVector<double,3> > &grad,
                    Vector<Tensor<double,3> > &hess)
{
  eval_multi_UBspline_3d_d_vgh (spline, r[0], r[1], r[2],
                                psi.data(),
                                (double*)grad.data(),
                                (double*)hess.data());
}

//////////////////////////////////
// Complex evaluation functions //
//////////////////////////////////
inline void
EinsplineMultiEval (multi_UBspline_3d_z *restrict spline,
                    TinyVector<double,3> r,
                    Vector<std::complex<double> > &psi)
{
  eval_multi_UBspline_3d_z (spline, r[0], r[1], r[2], psi.data());
}


inline void
EinsplineMultiEval (multi_UBspline_3d_z *restrict spline,
                    TinyVector<double,3> r,
                    Vector<std::complex<double> > &psi,
                    Vector<TinyVector<std::complex<double>,3> > &grad,
                    Vector<Tensor<std::complex<double>,3> > &hess)
{
  eval_multi_UBspline_3d_z_vgh (spline, r[0], r[1], r[2],
                                psi.data(),
                                (std::complex<double>*)grad.data(),
                                (std::complex<double>*)hess.data());
}


inline void
eval_multi_multi_UBspline_3d_cuda (multi_UBspline_3d_s_cuda *spline,
                                   float *pos, float *sign, float *phi[], int N)
{
  eval_multi_multi_UBspline_3d_s_sign_cuda  (spline, pos, sign, phi, N);
}


inline void
eval_multi_multi_UBspline_3d_cuda (multi_UBspline_3d_d_cuda *spline,
                                   double *pos, double *sign, double *phi[], int N)
{
  eval_multi_multi_UBspline_3d_d_sign_cuda  (spline, pos, sign, phi, N);
}


inline void
eval_multi_multi_UBspline_3d_cuda (multi_UBspline_3d_d_cuda *spline,
                                   double *pos, double *phi[], int N)
{
  eval_multi_multi_UBspline_3d_d_cuda  (spline, pos, phi, N);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda
(multi_UBspline_3d_s_cuda *spline, float *pos, float *sign, float Linv[],
 float *phi[], float *grad_lapl[], int N, int row_stride)
{
  eval_multi_multi_UBspline_3d_s_vgl_sign_cuda
  (spline, pos, sign, Linv, phi, grad_lapl, N, row_stride);
}


inline void eval_multi_multi_UBspline_3d_vgl_cuda
(multi_UBspline_3d_d_cuda *spline, double *pos, double *sign, double Linv[],
 double *phi[], double *grad_lapl[], int N, int row_stride)
{
  eval_multi_multi_UBspline_3d_d_vgl_sign_cuda
  (spline, pos, sign, Linv, phi, grad_lapl, N, row_stride);
}


inline void eval_multi_multi_UBspline_3d_vgl_cuda
(multi_UBspline_3d_d_cuda *spline, double *pos, double Linv[],
 double *phi[], double *grad_lapl[], int N, int row_stride)
{
  eval_multi_multi_UBspline_3d_d_vgl_cuda
  (spline, pos, Linv, phi, grad_lapl, N, row_stride);
}

// Complex CUDA routines
inline void
eval_multi_multi_UBspline_3d_cuda (multi_UBspline_3d_c_cuda *spline,
                                   float *pos, std::complex<float> *phi[], int N)
{
  eval_multi_multi_UBspline_3d_c_cuda  (spline, pos, phi, N);
}

inline void
eval_multi_multi_UBspline_3d_cuda (multi_UBspline_3d_z_cuda *spline,
                                   double *pos, std::complex<double> *phi[],
                                   int N)
{
  eval_multi_multi_UBspline_3d_z_cuda  (spline, pos, phi, N);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda
(multi_UBspline_3d_c_cuda *spline, float *pos, float Linv[],
 std::complex<float> *phi[], std::complex<float> *grad_lapl[], int N, int row_stride)
{
  eval_multi_multi_UBspline_3d_c_vgl_cuda
  (spline, pos, Linv, phi, grad_lapl, N, row_stride);
}

inline void eval_multi_multi_UBspline_3d_vgl_cuda
(multi_UBspline_3d_z_cuda *spline, double *pos, double Linv[],
 std::complex<double> *phi[], std::complex<double> *grad_lapl[], int N, int row_stride)
{
  eval_multi_multi_UBspline_3d_z_vgl_cuda
  (spline, pos, Linv, phi, grad_lapl, N, row_stride);
}

//////////////////////////////////////////////
// Vectorized evaluation routines using GPU //
//////////////////////////////////////////////

template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, int iat,
 gpu::device_vector<CudaRealType*> &phi)
{
  // app_log() << "Start EinsplineSet CUDA evaluation\n";
  int N = walkers.size();
  CudaRealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
    hostSign.resize(N);
    cudaSign.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = walkers[iw]->R[iat];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i=0; i<OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int) img;
    }
    int sign = 0;
    for (int i=0; i<OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign&1];
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  cudaSign = hostSign;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*)(cudapos.data()), cudaSign.data(), phi.data(), N);
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, int iat,
 gpu::device_vector<CudaRealType*> &phi)
{
  //    app_log() << "Eval 1.\n";
  int N = walkers.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = walkers[iw]->R[iat];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*)cudapos.data(), CudaValuePointers.data(), N);
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = walkers[iw]->R[iat];
  cudapos = hostPos;
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       (CudaRealType*)cudapos.data(),
                       (CudaRealType**)CudaValuePointers.data(),
                       phi.data(), CudaMultiSpline->num_splines,
                       walkers.size());
}


#ifdef QMC_COMPLEX
template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, int iat,
 gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " 
              << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, int iat,
 gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double> >::evaluate at line " << __LINE__ 
              << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}
#endif


template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaRealType*> &phi)
{
  // app_log() << "Start EinsplineSet CUDA evaluation\n";
  int N = newpos.size();
  CudaRealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
    hostSign.resize(N);
    cudaSign.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i=0; i<OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int) img;
    }
    int sign = 0;
    for (int i=0; i<OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign&1];
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  cudaSign = hostSign;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*)(cudapos.data()), cudaSign.data(),
   phi.data(), N);
  //app_log() << "End EinsplineSet CUDA evaluation\n";
}

template<typename T> void
EinsplineSetExtended<T>::resize_cuda(int numWalkers)
{
  CudaValuePointers.resize(numWalkers);
  CudaGradLaplPointers.resize(numWalkers);
  int N = CudaMultiSpline->num_splines;
  CudaValueVector.resize(N*numWalkers);
  CudaGradLaplVector.resize(4*N*numWalkers);
  gpu::host_vector<CudaStorageType*> hostValuePointers(numWalkers);
  gpu::host_vector<CudaStorageType*> hostGradLaplPointers(numWalkers);
  for (int i=0; i<numWalkers; i++)
  {
    hostValuePointers[i]    = &(CudaValueVector.data()[i*N]);
    hostGradLaplPointers[i] = &(CudaGradLaplVector.data()[4*i*N]);
  }
  CudaValuePointers    = hostValuePointers;
  CudaGradLaplPointers = hostGradLaplPointers;
  int M = MakeTwoCopies.size();
  CudaMakeTwoCopies.resize(M);
  CudaTwoCopiesIndex.resize(M);
  gpu::host_vector<int> hostMakeTwoCopies(M);
  gpu::host_vector<int> hostTwoCopiesIndex(M);
  int TwoCopiesIndexCounter = 0;
  for (int i=0; i<M; i++)
  {
    hostMakeTwoCopies[i] = MakeTwoCopies[i];
    hostTwoCopiesIndex[i] = TwoCopiesIndexCounter;
    TwoCopiesIndexCounter = MakeTwoCopies[i] ? TwoCopiesIndexCounter+2 : TwoCopiesIndexCounter+1;
  }
  CudaMakeTwoCopies = hostMakeTwoCopies;
  CudaTwoCopiesIndex = hostTwoCopiesIndex;
  CudakPoints.resize(M);
  CudakPoints_reduced.resize(M);
  gpu::host_vector<TinyVector<CUDA_PRECISION,OHMMS_DIM> > hostkPoints(M),
      hostkPoints_reduced(M);
  for (int i=0; i<M; i++)
  {
    //      PosType k_red1 = PrimLattice.toCart(kPoints[i]);
    PosType k_red2(dot(kPoints[i], PrimLattice.a(0)),
                   dot(kPoints[i], PrimLattice.a(1)),
                   dot(kPoints[i], PrimLattice.a(2)));
//       fprintf (stderr, "kred1 = %8.3f %8.3f %8.3f\n", k_red1[0], k_red1[1], k_red1[2]);
//       fprintf (stderr, "kred2 = %8.3f %8.3f %8.3f\n", k_red2[0], k_red2[1], k_red2[2]);
    for (int j=0; j<OHMMS_DIM; j++)
    {
      hostkPoints[i][j]         = kPoints[i][j];
      hostkPoints_reduced[i][j] = k_red2[j];
    }
  }
  CudakPoints = hostkPoints;
  CudakPoints_reduced = hostkPoints_reduced;
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaRealType*> &phi)
{
  //    app_log() << "Eval 2.\n";
  int N = walkers.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*)cudapos.data(), CudaValuePointers.data(), N);
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       (CudaRealType*)cudapos.data(),
                       (CudaRealType**)CudaValuePointers.data(),
                       phi.data(), CudaMultiSpline->num_splines,
                       walkers.size());
}


#ifdef QMC_COMPLEX
template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " 
              << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double>>::evaluate at line " << __LINE__ 
              << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}
#endif


template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaRealType*> &phi, gpu::device_vector<CudaRealType*> &grad_lapl,
 int row_stride)
{
  int N = walkers.size();
  CudaRealType plus_minus[2] = {1.0, -1.0};
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
    hostSign.resize(N);
    cudaSign.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i=0; i<OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int) img;
    }
    int sign = 0;
    for (int i=0; i<OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    hostSign[iw] = plus_minus[sign&1];
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  cudaSign = hostSign;
  eval_multi_multi_UBspline_3d_vgl_cuda
  (CudaMultiSpline, (CudaRealType*)cudapos.data(), cudaSign.data(),
   Linv_cuda.data(), phi.data(), grad_lapl.data(), N, row_stride);
  //gpu::host_vector<CudaRealType*> pointers;
  //pointers = phi;
  //CudaRealType data[N];
  //cudaMemcpy (data, pointers[0], N*sizeof(CudaRealType), cudaMemcpyDeviceToHost);
  //for (int i=0; i<N; i++)
  //  fprintf (stderr, "%1.12e\n", data[i]);
}


template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaRealType*> &phi, gpu::device_vector<CudaRealType*> &grad_lapl,
 int row_stride)
{
  //    app_log() << "Eval 3.\n";
  int N = walkers.size();
  int M = CudaMultiSpline->num_splines;
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  eval_multi_multi_UBspline_3d_vgl_cuda
  (CudaMultiSpline, (CudaRealType*)cudapos.data(),  Linv_cuda.data(), CudaValuePointers.data(),
   CudaGradLaplPointers.data(), N, CudaMultiSpline->num_splines);
  // DEBUG
  //  TinyVector<double,OHMMS_DIM> r(hostPos[0][0], hostPos[0][1], hostPos[0][2]);
  //  Vector<std::complex<double > > psi(M);
  //  Vector<TinyVector<std::complex<double>,3> > grad(M);
  //  Vector<Tensor<std::complex<double>,3> > hess(M);
  //  EinsplineMultiEval (MultiSpline, r, psi, grad, hess);
  //  // std::complex<double> cpuSpline[M];
  //  // TinyVector<std::complex<double>,OHMMS_DIM> std::complex<double> cpuGrad[M];
  //  // Tensor cpuHess[M];
  //  // eval_multi_UBspline_3d_z_vgh (MultiSpline, hostPos[0][0], hostPos[0][1], hostPos[0][2],
  //  // 				  cpuSpline);
  //  gpu::host_vector<CudaStorageType*> pointers;
  //  pointers = CudaGradLaplPointers;
  //  std::complex<CudaRealType> gpuSpline[4*M];
  //  cudaMemcpy(gpuSpline, pointers[10], 4*M * sizeof(std::complex<CudaRealType>), cudaMemcpyDeviceToHost);
  //  for (int i=0; i<M; i++)
  //    fprintf (stderr, "real: %10.6e %10.6e %10.6e , imag: %10.6e %10.6e %10.6e .\n",
  //             trace(hess[i],GGt).real(), gpuSpline[3*M+i].real(), trace(hess[i],GGt).real() - gpuSpline[3*M+i].real(),
  //             trace(hess[i], GGt).imag(), gpuSpline[3*M+i].imag(), trace(hess[i], GGt).imag() - gpuSpline[3*M+i].imag());
  //  fprintf (stderr, "\n");
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;

  /* Original implementation
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       (CudaRealType*)cudapos.data(),
                       (CudaRealType**)CudaValuePointers.data(), phi.data(),
                       (CudaRealType**)CudaGradLaplPointers.data(), grad_lapl.data(),
                       CudaMultiSpline->num_splines,  walkers.size(), row_stride);
  */
  // Ye: optimized memory access.
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       CudaTwoCopiesIndex.data(),
                       (CudaRealType*)cudapos.data(),
                       (CudaRealType**)CudaValuePointers.data(), phi.data(),
                       (CudaRealType**)CudaGradLaplPointers.data(), grad_lapl.data(),
                       CudaMultiSpline->num_splines,  walkers.size(), row_stride);
}


#ifdef QMC_COMPLEX
template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaComplexType*> &phi,
 gpu::device_vector<CudaComplexType*> &grad_lapl,
 int row_stride)
{
  app_error() << "Code should not arrive at this point: "
              << "EinsplineSetExtended<double>::evaluate at line " 
              << __LINE__ << " in file " << __FILE__ << "\n";
  abort();
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
 gpu::device_vector<CudaComplexType*> &phi,
 gpu::device_vector<CudaComplexType*> &grad_lapl,
 int row_stride)
{
  int N = walkers.size();
  int M = CudaMultiSpline->num_splines;
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = newpos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  eval_multi_multi_UBspline_3d_vgl_cuda
  (CudaMultiSpline, (CudaRealType*)cudapos.data(),  Linv_cuda.data(), CudaValuePointers.data(),
   CudaGradLaplPointers.data(), N, CudaMultiSpline->num_splines);
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       (CudaRealType*) cudapos.data(),
                       (CudaValueType**) CudaValuePointers.data(),
                       (CudaValueType**) phi.data(),
                       (CudaValueType**) CudaGradLaplPointers.data(),
                       (CudaValueType**) grad_lapl.data(),
                       CudaMultiSpline->num_splines, walkers.size(), row_stride);
}
#endif


template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  int N = pos.size();
  CudaRealType plus_minus[2] = {1.0, -1.0};
  if (NLcudapos.size() < N)
  {
    NLhostPos.resize(N);
    NLcudapos.resize(N);
    NLhostSign.resize(N);
    NLcudaSign.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    int image[OHMMS_DIM];
    for (int i=0; i<OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      image[i] = (int) img;
    }
    int sign = 0;
    for (int i=0; i<OHMMS_DIM; i++)
      sign += HalfG[i] * image[i];
    NLhostSign[iw] = plus_minus[sign&1];
    NLhostPos[iw] = ru;
  }
  NLcudapos  = NLhostPos;
  NLcudaSign = NLhostSign;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*)(NLcudapos.data()),
   NLcudaSign.data(), phi.data(), N);
}

template<> void
EinsplineSetExtended<double>::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetExtended<std::complex<double> >::evaluate at " 
              << __LINE__ << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
}



template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  //    app_log() << "Eval 4.\n";
  int N = pos.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*) cudapos.data(),
   CudaValuePointers.data(), N);
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;
  apply_phase_factors ((CudaRealType*) CudakPoints.data(),
                       CudaMakeTwoCopies.data(),
                       (CudaRealType*)cudapos.data(),
                       (CudaRealType**)CudaValuePointers.data(),
                       phi.data(), CudaMultiSpline->num_splines, N);
}

template<> void
EinsplineSetExtended<std::complex<double> >::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
#ifdef QMC_COMPLEX
  int N = pos.size();
  if (CudaValuePointers.size() < N)
    resize_cuda(N);
  if (cudapos.size() < N)
  {
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
  {
    PosType r = pos[iw];
    PosType ru(PrimLattice.toUnit(r));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= std::floor (ru[i]);
    hostPos[iw] = ru;
  }
  cudapos = hostPos;
// AT debug:
//  std::cout << "# splines: " << CudaMultiSpline->num_splines << "\n";
  eval_multi_multi_UBspline_3d_cuda
  (CudaMultiSpline, (CudaRealType*) cudapos.data(),
   CudaValuePointers.data(), N);
  // Now, add on phases
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;
  apply_phase_factors((CudaRealType*) CudakPoints.data(),
                      (CudaRealType*) cudapos.data(),
                      CudaValuePointers.data(),
                      phi.data(),
                      CudaMultiSpline->num_splines, N);
// AT debug:
/*  gpu::host_vector<CudaValueType*> pointers;
  pointers = CudaValuePointers;
  CudaValueType data[N], data_new[N];
  cudaMemcpy (data, pointers[0], N*sizeof(CudaValueType), cudaMemcpyDeviceToHost);
  pointers = phi;
  cudaMemcpy (data_new, pointers[0], N*sizeof(CudaValueType), cudaMemcpyDeviceToHost);
  std::cout << "CudaValuePointers -> phi (# splines: " << CudaMultiSpline->num_splines << "):\n";
  for (int i=0; i<N; i++)
    std::cout << i << ": " << data[i].real() << " + " << data[i].imag() << "i -> " << data_new[i].real() << " + " << data_new[i].imag() << "i\n";
  std::cout.flush();
  abort();*/
#else
  app_error() << "EinsplineSetExtended<std::complex<double> >::evaluate at " 
              << __LINE__ << " in file " << __FILE__
              << " not yet implemented.\n";
  abort();
#endif

}


////////////////////////////////
// Hybrid evaluation routines //
////////////////////////////////
// template<typename StorageType> void
// EinsplineSetHybrid<StorageType>::sort_electrons(std::vector<PosType> &pos)
// {
//   int nw = pos.size();
//   if (nw > CurrentWalkers)
//     resize_cuda(nw);

//   AtomicPolyJobs_CPU.clear();
//   AtomicSplineJobs_CPU.clear();
//   rhats_CPU.clear();
//   int numAtomic = 0;

//   // First, sort electrons into three categories:
//   // 1) Interstitial region with 3D-Bsplines
//   // 2) Atomic region near origin:      polynomial  radial functions
//   // 3) Atomic region not near origin:  1D B-spline radial functions

//   for (int i=0; i<newpos.size(); i++) {
//     PosType r = newpos[i];
//     // Note: this assumes that the atomic radii are smaller than the simulation cell radius.
//     for (int j=0; j<AtomicOrbitals.size(); j++) {
// 	AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[j];
// 	PosType dr = r - orb.Pos;
// 	PosType u = PrimLattice.toUnit(dr);
// 	for (int k=0; k<OHMMS_DIM; k++)
// 	  u[k] -= round(u[k]);
// 	dr = PrimLattice.toCart(u);
// 	RealType dist2 = dot (dr,dr);
// 	if (dist2 < orb.PolyRadius * orb.PolyRadius) {
// 	  AtomicPolyJob<CudaRealType> job;
// 	  RealType dist = std::sqrt(dist2);
// 	  job.dist = dist;
// 	  RealType distInv = 1.0/dist;
// 	  for (int k=0; k<OHMMS_DIM; k++) {
// 	    CudaRealType x = distInv*dr[k];
// 	    job.rhat[k] = distInv * dr[k];
// 	    rhats_CPU.push_back(x);
// 	  }
// 	  job.lMax = orb.lMax;
// 	  job.YlmIndex  = i;
// 	  job.PolyOrder = orb.PolyOrder;
// 	  //job.PolyCoefs = orb.PolyCoefs;
// 	  AtomicPolyJobs_CPU.push_back(job);
// 	  numAtomic++;
// 	}
// 	else if (dist2 < orb.CutoffRadius*orb.CutoffRadius) {
// 	  AtomicSplineJob<CudaRealType> job;
// 	   RealType dist = std::sqrt(dist2);
// 	  job.dist = dist;
// 	  RealType distInv = 1.0/dist;
// 	  for (int k=0; k<OHMMS_DIM; k++) {
// 	    CudaRealType x = distInv*dr[k];
// 	    job.rhat[k] = distInv * dr[k];
// 	    rhats_CPU.push_back(x);
// 	  }
// 	  job.lMax      = orb.lMax;
// 	  job.YlmIndex  = i;
// 	  job.phi       = phi[i];
// 	  job.grad_lapl = grad_lapl[i];
// 	  //job.PolyCoefs = orb.PolyCoefs;
// 	  AtomicSplineJobs_CPU.push_back(job);
// 	  numAtomic++;
// 	}
// 	else { // Regular 3D B-spline job
// 	  BsplinePos_CPU.push_back (r);

// 	}
//     }
//   }


// }

template<> void
EinsplineSetHybrid<double>::resize_cuda(int numwalkers)
{
  EinsplineSetExtended<double>::resize_cuda(numwalkers);
  CurrentWalkers = numwalkers;
  // Resize Ylm temporaries
  // Find lMax;
  lMax=-1;
  for (int i=0; i<AtomicOrbitals.size(); i++)
    lMax = std::max(AtomicOrbitals[i].lMax, lMax);
  numlm = (lMax+1)*(lMax+1);
  Ylm_BS = ((numlm+15)/16) * 16;
  Ylm_GPU.resize(numwalkers*Ylm_BS*3);
  Ylm_ptr_GPU.resize        (numwalkers);
  Ylm_ptr_CPU.resize       (numwalkers);
  dYlm_dtheta_ptr_GPU.resize(numwalkers);
  dYlm_dtheta_ptr_CPU.resize(numwalkers);
  dYlm_dphi_ptr_GPU.resize  (numwalkers);
  dYlm_dphi_ptr_CPU.resize  (numwalkers);
  rhats_CPU.resize(OHMMS_DIM*numwalkers);
  rhats_GPU.resize(OHMMS_DIM*numwalkers);
  HybridJobs_GPU.resize(numwalkers);
  HybridData_GPU.resize(numwalkers);
  for (int iw=0; iw<numwalkers; iw++)
  {
    Ylm_ptr_CPU[iw]         = Ylm_GPU.data() + (3*iw+0)*Ylm_BS;
    dYlm_dtheta_ptr_CPU[iw] = Ylm_GPU.data() + (3*iw+1)*Ylm_BS;
    dYlm_dphi_ptr_CPU[iw]   = Ylm_GPU.data() + (3*iw+2)*Ylm_BS;
  }
  Ylm_ptr_GPU         = Ylm_ptr_CPU;
  dYlm_dtheta_ptr_GPU = dYlm_dtheta_ptr_CPU;
  dYlm_dphi_ptr_GPU   = dYlm_dphi_ptr_CPU;
  // Resize AtomicJob temporaries
  // AtomicPolyJobs_GPU.resize(numwalkers);
  // AtomicSplineJobs_GPU.resize(numwalkers);
}


template<> void
EinsplineSetHybrid<std::complex<double> >::resize_cuda(int numwalkers)
{
  EinsplineSetExtended<std::complex<double> >::resize_cuda(numwalkers);
  CurrentWalkers = numwalkers;
  // Resize Ylm temporaries
  // Find lMax;
  lMax=-1;
  for (int i=0; i<AtomicOrbitals.size(); i++)
    lMax = std::max(AtomicOrbitals[i].lMax, lMax);
  numlm = (lMax+1)*(lMax+1);
  Ylm_BS = ((2*numlm+15)/16) * 16;
  Ylm_GPU.resize(numwalkers*Ylm_BS*3);
  Ylm_ptr_GPU.resize        (numwalkers);
  Ylm_ptr_CPU.resize       (numwalkers);
  dYlm_dtheta_ptr_GPU.resize(numwalkers);
  dYlm_dtheta_ptr_CPU.resize(numwalkers);
  dYlm_dphi_ptr_GPU.resize  (numwalkers);
  dYlm_dphi_ptr_CPU.resize  (numwalkers);
  rhats_CPU.resize(OHMMS_DIM*numwalkers);
  rhats_GPU.resize(OHMMS_DIM*numwalkers);
  HybridJobs_GPU.resize(numwalkers);
  HybridData_GPU.resize(numwalkers);
  for (int iw=0; iw<numwalkers; iw++)
  {
    Ylm_ptr_CPU[iw]         = Ylm_GPU.data() + (3*iw+0)*Ylm_BS;
    dYlm_dtheta_ptr_CPU[iw] = Ylm_GPU.data() + (3*iw+1)*Ylm_BS;
    dYlm_dphi_ptr_CPU[iw]   = Ylm_GPU.data() + (3*iw+2)*Ylm_BS;
  }
  Ylm_ptr_GPU         = Ylm_ptr_CPU;
  dYlm_dtheta_ptr_GPU = dYlm_dtheta_ptr_CPU;
  dYlm_dphi_ptr_GPU   = dYlm_dphi_ptr_CPU;
  // Resize AtomicJob temporaries
  // AtomicPolyJobs_GPU.resize(numwalkers);
  // AtomicSplineJobs_GPU.resize(numwalkers);
}


// Vectorized evaluation functions
template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,
                                      gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << " gpu::device_vector<CudaRealType*> &phi) not implemented.\n";
  abort();
}


template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,
                                      gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << " gpu::device_vector<CudaComplexType*> &phi) not implemented.\n";
  abort();
}


template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                                      gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << " (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CudaRealType*> &phi) not implemented.\n";
  abort();
}

template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                                      gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "  (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,\n"
              << "   gpu::device_vector<CudaComplexType*> &phi) not implemented.\n";
  abort();
}



template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                                      gpu::device_vector<CudaRealType*> &phi,
                                      gpu::device_vector<CudaRealType*> &grad_lapl,
                                      int row_stride)
{
  int N = newpos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;
  // hostPos = cudapos;
  // for (int i=0; i<newpos.size(); i++)
  //   std::cerr << "newPos[" << i << "] = " << newpos[i] << std::endl;
  // gpu::host_vector<CudaRealType> IonPos_CPU(IonPos_GPU.size());
  // IonPos_CPU = IonPos_GPU;
  // for (int i=0; i<IonPos_CPU.size()/3; i++)
  //   fprintf (stderr, "ion[%d] = [%10.6f %10.6f %10.6f]\n",
  // 	       i, IonPos_CPU[3*i+0], IonPos_CPU[3*i+2], IonPos_CPU[3*i+2]);
  // std::cerr << "cudapos.size()        = " << cudapos.size() << std::endl;
  // std::cerr << "IonPos.size()         = " << IonPos_GPU.size() << std::endl;
  // std::cerr << "PolyRadii.size()      = " << PolyRadii_GPU.size() << std::endl;
  // std::cerr << "CutoffRadii.size()    = " << CutoffRadii_GPU.size() << std::endl;
  // std::cerr << "AtomicOrbitals.size() = " << AtomicOrbitals.size() << std::endl;
  // std::cerr << "L_cuda.size()         = " << L_cuda.size() << std::endl;
  // std::cerr << "Linv_cuda.size()      = " << Linv_cuda.size() << std::endl;
  // std::cerr << "HybridJobs_GPU.size() = " << HybridJobs_GPU.size() << std::endl;
  // std::cerr << "rhats_GPU.size()      = " << rhats_GPU.size() << std::endl;

  MakeHybridJobList<CudaRealType> ((CudaRealType*)cudapos.data(), N, IonPos_GPU.data(),
                                   PolyRadii_GPU.data(), CutoffRadii_GPU.data(),
                                   AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(),
                                   HybridData_GPU.data());

  CalcYlmRealCuda<CudaRealType> (rhats_GPU.data(), HybridJobs_GPU.data(),
                                 Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
                                 dYlm_dphi_ptr_GPU.data(), lMax, newpos.size());

  evaluate3DSplineReal<CudaRealType> (HybridJobs_GPU.data(), (CudaRealType*)cudapos.data(),
                                      (CudaRealType*)CudakPoints_reduced.data(),CudaMultiSpline,
                                      Linv_cuda.data(), phi.data(), grad_lapl.data(),
                                      row_stride, NumOrbitals, newpos.size());

  evaluateHybridSplineReal<CudaRealType> (HybridJobs_GPU.data(), rhats_GPU.data(), 
                                          Ylm_ptr_GPU.data(),
                                          dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(),
                                          AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
                                          (CudaRealType*)CudakPoints_reduced.data(),
                                          phi.data(), grad_lapl.data(),
                                          row_stride, NumOrbitals, newpos.size(), lMax);
#ifdef HYBRID_DEBUG
  gpu::host_vector<CudaRealType*> phi_CPU (phi.size()), grad_lapl_CPU(phi.size());
  phi_CPU = phi;
  grad_lapl_CPU = grad_lapl;
  gpu::host_vector<CudaRealType> vals_CPU(NumOrbitals), GL_CPU(4*row_stride);
  gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  gpu::host_vector<HybridData<CudaRealType> > HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  rhats_CPU = rhats_GPU;
  for (int iw=0; iw<newpos.size(); iw++)
    if (false && HybridJobs_CPU[iw] == ATOMIC_POLY_JOB)
    {
      ValueVector_t CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector_t CPUgrad(NumOrbitals);
      HybridData<CudaRealType> &d = HybridData_CPU[iw];
      AtomicOrbital<double> &atom = AtomicOrbitals[d.ion];
      atom.evaluate (newpos[iw], CPUvals, CPUgrad, CPUlapl);
      cudaMemcpy (&vals_CPU[0], phi_CPU[iw], NumOrbitals*sizeof(float),
                  cudaMemcpyDeviceToHost);
      cudaMemcpy (&GL_CPU[0], grad_lapl_CPU[iw], 4*row_stride*sizeof(float),
                  cudaMemcpyDeviceToHost);
      // fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
      // 	 iw, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
      double mindist = 1.0e5;
      for (int ion=0; ion<AtomicOrbitals.size(); ion++)
      {
        PosType disp = newpos[iw] - AtomicOrbitals[ion].Pos;
        PosType u = PrimLattice.toUnit(disp);
        PosType img;
        for (int i=0; i<OHMMS_DIM; i++)
          u[i] -= round(u[i]);
        disp = PrimLattice.toCart(u);
        double dist = std::sqrt(dot(disp,disp));
        if (dist < AtomicOrbitals[ion].CutoffRadius)
          mindist = dist;
      }
      if (std::fabs (mindist - d.dist) > 1.0e-3)
        fprintf (stderr, "CPU dist = %1.8f  GPU dist = %1.8f\n",
                 mindist, d.dist);
      for (int j=0; j<NumOrbitals; j++)
      {
        //	  if (isnan(vals_CPU[j])) {
        if (true || isnan(GL_CPU[0*row_stride+j]))
        {
          std::cerr << "iw = " << iw << std::endl;
          fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f] dist = %10.6e\n",
                   iw, rhats_CPU[3*iw+0], rhats_CPU[3*iw+1], rhats_CPU[3*iw+2],
                   d.dist);
          fprintf (stderr, "val[%2d]  = %10.5e %10.5e\n",
                   j, vals_CPU[j], CPUvals[j]);
          fprintf (stderr, "grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j,
                   GL_CPU[0*row_stride+j], CPUgrad[j][0],
                   GL_CPU[1*row_stride+j], CPUgrad[j][1],
                   GL_CPU[2*row_stride+j], CPUgrad[j][2]);
          fprintf (stderr, "lapl[%2d] = %10.5e %10.5e\n",
                   j, GL_CPU[3*row_stride+j], CPUlapl[j]);
        }
      }
    }
    else if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
    {
      std::cerr << "HalfG = " << HalfG << std::endl;
      ValueVector_t CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector_t CPUgrad(NumOrbitals);
      PosType ru(PrimLattice.toUnit(newpos[iw]));
      PosType img;
      int sign = 0;
      for (int i=0; i<3; i++)
      {
        img[i] = std::floor(ru[i]);
        ru[i] -= img[i];
        sign += HalfG[i]*(int)img[i];
      }
      EinsplineMultiEval (MultiSpline, ru, CPUvals, CPUgrad,
                          StorageHessVector);
      cudaMemcpy (&vals_CPU[0], phi_CPU[iw], NumOrbitals*sizeof(float),
                  cudaMemcpyDeviceToHost);
      cudaMemcpy (&GL_CPU[0], grad_lapl_CPU[iw], 4*row_stride*sizeof(float),
                  cudaMemcpyDeviceToHost);
      for (int j=0; j<NumOrbitals; j++)
      {
        CPUgrad[j] = dot (PrimLattice.G, CPUgrad[j]);
        CPUlapl[j] = trace (StorageHessVector[j], GGt);
        if (sign & 1)
        {
          CPUvals[j] *= -1.0;
          CPUgrad[j] *= -1.0;
          CPUlapl[j] *= -1.0;
        }
        fprintf (stderr, "\nGPU=%10.6f  %10.6f %10.6f %10.6f  %10.6f\n",
                 vals_CPU[j], GL_CPU[0*row_stride+j], GL_CPU[1*row_stride+j],
                 GL_CPU[2*row_stride+j], GL_CPU[3*row_stride+j]);
        fprintf (stderr, "CPU=%10.6f  %10.6f %10.6f %10.6f  %10.6f sign = %d\n",
                 CPUvals[j], CPUgrad[j][0], CPUgrad[j][1], CPUgrad[j][2],
                 CPUlapl[j], sign);
        if (std::isnan(GL_CPU[0*row_stride+j]))
        {
          std::cerr << "r[" << iw << "] = " << newpos[iw] << std::endl;
          std::cerr << "iw = " << iw << std::endl;
          fprintf (stderr, "3D val[%2d]  = %10.5e %10.5e\n",
                   j, vals_CPU[j], CPUvals[j]);
          fprintf (stderr, "3D grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j,
                   GL_CPU[0*row_stride+j], CPUgrad[j][0],
                   GL_CPU[1*row_stride+j], CPUgrad[j][1],
                   GL_CPU[2*row_stride+j], CPUgrad[j][2]);
          fprintf (stderr, "3D lapl[%2d] = %10.5e %10.5e\n",
                   j, GL_CPU[3*row_stride+j], CPUlapl[j]);
        }
      }
    }
  gpu::host_vector<float> Ylm_CPU(Ylm_GPU.size());
  Ylm_CPU = Ylm_GPU;
  rhats_CPU = rhats_GPU;
  for (int i=0; i<rhats_CPU.size()/3; i++)
    fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n",
             i, rhats_CPU[3*i+0], rhats_CPU[3*i+1], rhats_CPU[3*i+2]);
  //    gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  //    gpu::host_vector<HybridDataFloat> HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  std::cerr << "Before loop.\n";
  for (int i=0; i<newpos.size(); i++)
    if (HybridJobs_CPU[i] != BSPLINE_3D_JOB)
    {
      std::cerr << "Inside if.\n";
      PosType rhat(rhats_CPU[3*i+0], rhats_CPU[3*i+1], rhats_CPU[3*i+2]);
      AtomicOrbital<double> &atom = AtomicOrbitals[HybridData_CPU[i].ion];
      int numlm = (atom.lMax+1)*(atom.lMax+1);
      std::vector<double> Ylm(numlm), dYlm_dtheta(numlm), dYlm_dphi(numlm);
      atom.CalcYlm (rhat, Ylm, dYlm_dtheta, dYlm_dphi);
      for (int lm=0; lm < numlm; lm++)
      {
        fprintf (stderr, "lm=%3d  Ylm_CPU=%8.5f  Ylm_GPU=%8.5f\n",
                 lm, Ylm[lm], Ylm_CPU[3*i*Ylm_BS+lm]);
      }
    }
  fprintf (stderr, " N  img      dist    ion    lMax\n");
  for (int i=0; i<HybridData_CPU.size(); i++)
  {
    HybridData<CudaRealType> &d = HybridData_CPU[i];
    fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
             i, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
  }
#endif
  // int N = newpos.size();
  // if (N > CurrentWalkers)
  //   resize_cuda(N);
  // AtomicPolyJobs_CPU.clear();
  // AtomicSplineJobs_CPU.clear();
  // rhats_CPU.clear();
  // BsplinePos_CPU.clear();
  // BsplineVals_CPU.clear();
  // BsplineGradLapl_CPU.clear();
  // int numAtomic = 0;
  // // First, sort electrons into three categories:
  // // 1) Interstitial region with 3D B-splines
  // // 2) Atomic region near origin:      polynomial  radial functions
  // // 3) Atomic region not near origin:  1D B-spline radial functions
  // for (int i=0; i<newpos.size(); i++) {
  //   PosType r = newpos[i];
  //   // Note: this assumes that the atomic radii are smaller than the simulation cell radius.
  //   for (int j=0; j<AtomicOrbitals.size(); j++) {
  // 	AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[j];
  // 	PosType dr = r - orb.Pos;
  // 	PosType u = PrimLattice.toUnit(dr);
  // 	for (int k=0; k<OHMMS_DIM; k++)
  // 	  u[k] -= round(u[k]);
  // 	dr = PrimLattice.toCart(u);
  // 	RealType dist2 = dot (dr,dr);
  // 	if (dist2 < orb.PolyRadius * orb.PolyRadius) {
  // 	  AtomicPolyJob<CudaRealType> job;
  // 	  RealType dist = std::sqrt(dist2);
  // 	  job.dist = dist;
  // 	  RealType distInv = 1.0/dist;
  // 	  for (int k=0; k<OHMMS_DIM; k++) {
  // 	    CudaRealType x = distInv*dr[k];
  // 	    job.rhat[k] = distInv * dr[k];
  // 	    rhats_CPU.push_back(x);
  // 	  }
  // 	  job.lMax = orb.lMax;
  // 	  job.YlmIndex  = i;
  // 	  job.PolyOrder = orb.PolyOrder;
  // 	  //job.PolyCoefs = orb.PolyCoefs;
  // 	  AtomicPolyJobs_CPU.push_back(job);
  // 	  numAtomic++;
  // 	}
  // 	else if (dist2 < orb.CutoffRadius*orb.CutoffRadius) {
  // 	  AtomicSplineJob<CudaRealType> job;
  // 	   RealType dist = std::sqrt(dist2);
  // 	  job.dist = dist;
  // 	  RealType distInv = 1.0/dist;
  // 	  for (int k=0; k<OHMMS_DIM; k++) {
  // 	    CudaRealType x = distInv*dr[k];
  // 	    job.rhat[k] = distInv * dr[k];
  // 	    rhats_CPU.push_back(x);
  // 	  }
  // 	  job.lMax      = orb.lMax;
  // 	  job.YlmIndex  = i;
  // 	  job.phi       = phi[i];
  // 	  job.grad_lapl = grad_lapl[i];
  // 	  //job.PolyCoefs = orb.PolyCoefs;
  // 	  AtomicSplineJobs_CPU.push_back(job);
  // 	  numAtomic++;
  // 	}
  // 	else { // Regular 3D B-spline job
  // 	  BsplinePos_CPU
  // 	}
  //   }
  // }
  // //////////////////////////////////
  // // First, evaluate 3D B-splines //
  // //////////////////////////////////
  // int N = newpos.size();
  // CudaRealType plus_minus[2] = {1.0, -1.0};
  // if (cudapos.size() < N) {
  //   hostPos.resize(N);
  //   cudapos.resize(N);
  //   hostSign.resize(N);
  //   cudaSign.resize(N);
  // }
  // for (int iw=0; iw < N; iw++) {
  //   PosType r = newpos[iw];
  //   PosType ru(PrimLattice.toUnit(r));
  //   int image[OHMMS_DIM];
  //   for (int i=0; i<OHMMS_DIM; i++) {
  // 	RealType img = std::floor(ru[i]);
  // 	ru[i] -= img;
  // 	image[i] = (int) img;
  //   }
  //   int sign = 0;
  //   for (int i=0; i<OHMMS_DIM; i++)
  // 	sign += HalfG[i] * image[i];
  //   hostSign[iw] = plus_minus[sign&1];
  //   hostPos[iw] = ru;
  // }
  // cudapos = hostPos;
  // cudaSign = hostSign;
  // eval_multi_multi_UBspline_3d_cuda
  //   (CudaMultiSpline, (CudaRealType*)(cudapos.data()), cudaSign.data(),
  //    phi.data(), N);
  // ////////////////////////////////////////////////////////////
  // // Next, evaluate spherical harmonics for atomic orbitals //
  // ////////////////////////////////////////////////////////////
  // // Evaluate Ylms
  // if (rhats_CPU.size()) {
  //   rhats_GPU = rhats_CPU;
  //   CalcYlmComplexCuda(rhats_GPU.data(), Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
  // 			 dYlm_dphi_ptr_GPU.data(), lMax, numAtomic);
  //   std::cerr << "Calculated Ylms.\n";
  // }
  // ///////////////////////////////////////
  // // Next, evaluate 1D spline orbitals //
  // ///////////////////////////////////////
  // if (AtomicSplineJobs_CPU.size()) {
  //   AtomicSplineJobs_GPU = AtomicSplineJobs_CPU;
  // }
  // ///////////////////////////////////////////
  // // Next, evaluate 1D polynomial orbitals //
  // ///////////////////////////////////////////
  // if (AtomicSplineJobs_CPU.size()) {
  //   AtomicPolyJobs_GPU = AtomicPolyJobs_CPU;
  // }
}

template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                                      gpu::device_vector<CudaComplexType*> &phi,
                                      gpu::device_vector<CudaComplexType*> &grad_lapl,
                                      int row_stride)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CudaComplexType*> &phi,\n"
              << " gpu::device_vector<CudaComplexType*> &grad_lapl, int row_stride)\n"
              << "     is not yet implemented.\n";
  abort();
}

template<> void
EinsplineSetHybrid<double>::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  int N = pos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CudaRealType> ((CudaRealType*)cudapos.data(), N, IonPos_GPU.data(),
                                   PolyRadii_GPU.data(), CutoffRadii_GPU.data(),
                                   AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(),
                                   HybridData_GPU.data());

  CalcYlmRealCuda<CudaRealType> (rhats_GPU.data(), HybridJobs_GPU.data(),
                                 Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
                                 dYlm_dphi_ptr_GPU.data(), lMax, pos.size());

  evaluateHybridSplineReal<CudaRealType> (HybridJobs_GPU.data(), Ylm_ptr_GPU.data(),
                                          AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
                                          (CudaRealType*)CudakPoints_reduced.data(),
                                          phi.data(), NumOrbitals, pos.size(), lMax);

  evaluate3DSplineReal<CudaRealType> (HybridJobs_GPU.data(), (CudaRealType*)cudapos.data(),
                                      (CudaRealType*)CudakPoints_reduced.data(),CudaMultiSpline,
                                      Linv_cuda.data(), phi.data(), NumOrbitals, pos.size());
}

template<> void
EinsplineSetHybrid<double>::evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<double>::evaluate \n"
              << "(std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)\n"
              << "     is not yet implemented.\n";
  abort();
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, int iat,
    gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << "			                            gpu::device_vector<CudaRealType*> &phi)\n"
              << "not yet implemented.\n";
}


template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, int iat,
    gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, int iat,\n"
              << "			                            gpu::device_vector<CudaComplexType*> &phi)\n"
              << "not yet implemented.\n";
}


template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
    gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,\n"
              << "			                            gpu::device_vector<CudaRealType*> &phi)\n"
              << "not yet implemented.\n";
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
    gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> ,\n"
              << "			                            gpu::device_vector<CudaComplexType*> &phi)\n"
              << "not yet implemented.\n";
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate
(std::vector<Walker_t*> &walkers,  std::vector<PosType> &newpos,
 gpu::device_vector<CudaRealType*> &phi, gpu::device_vector<CudaRealType*> &grad_lapl,
 int row_stride)
{
  static int numAtomic=0;
  static int num3D=0;
  int N = newpos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = newpos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CudaRealType> ((CudaRealType*)cudapos.data(), N, IonPos_GPU.data(),
                                   PolyRadii_GPU.data(), CutoffRadii_GPU.data(),
                                   AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(),
                                   HybridData_GPU.data());

  CalcYlmComplexCuda<CudaRealType> (rhats_GPU.data(), HybridJobs_GPU.data(),
                                    Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
                                    dYlm_dphi_ptr_GPU.data(), lMax, newpos.size());

  evaluate3DSplineComplexToReal<CudaRealType> (HybridJobs_GPU.data(),
                                               (CudaRealType*)cudapos.data(),
                                               (CudaRealType*)CudakPoints.data(),
                                               CudaMakeTwoCopies.data(), CudaMultiSpline,
                                               Linv_cuda.data(),
                                               phi.data(), grad_lapl.data(),
                                               row_stride, CudaMakeTwoCopies.size(),
                                               newpos.size());

  evaluateHybridSplineComplexToReal<CudaRealType>
                      (HybridJobs_GPU.data(), rhats_GPU.data(),
                       Ylm_ptr_GPU.data(),
                       dYlm_dtheta_ptr_GPU.data(), dYlm_dphi_ptr_GPU.data(),
                       AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
                       (CudaRealType*)CudakPoints_reduced.data(),
                       CudaMakeTwoCopies.data(), (CudaRealType**)phi.data(),
                       grad_lapl.data(), row_stride,
                       CudaMakeTwoCopies.size(), newpos.size(), lMax);

#ifdef HYBRID_DEBUG
  // gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  // HybridJobs_CPU = HybridJobs_GPU;
  // int M = MakeTwoCopies.size();
  // // ComplexValueVector_t CPUzvals(M), CPUzlapl(M);
  // // ComplexGradVector_t CPUzgrad(M);
  // for (int iw=0; iw<newpos.size(); iw++) {
  //   if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  // 	num3D++;
  //   else
  // 	numAtomic++;
  //   // bool atomic=false;
  //   // for (int ion=0; ion<AtomicOrbitals.size(); ion++)
  //   // 	if (AtomicOrbitals[ion].evaluate
  //   // 	    (newpos[iw], CPUzvals, CPUzgrad, CPUzlapl)) {
  //   // 	  atomic = true;
  //   // 	  break;
  //   // 	}
  //   // if (atomic)
  //   // 	numAtomic++;
  //   // else
  //   // 	num3D++;
  //   // if (atomic && HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   // 	cerr << "Error!  Used BSPLINE_3D when I should have used atomic.\n";
  //   // else if (!atomic && HybridJobs_CPU[iw] != BSPLINE_3D_JOB)
  //   // 	cerr << "Error!  Used atomic when I should have used 3D.\n";
  // }
  // // fprintf (stderr, "Num atomic = %d  Num 3D = %d\n",
  // // 	     numAtomic, num3D);
  // if (numAtomic + num3D > 100000) {
  //   fprintf (stderr, "Num atomic = %d  Num 3D = %d\n",
  // 	       numAtomic, num3D);
  //   fprintf (stderr, "Percent atomic = %1.5f\%\n",
  // 	       100.0*(double)numAtomic / (double)(numAtomic+num3D));
  //   numAtomic = num3D = 0;
  // }
  gpu::host_vector<CudaRealType*> phi_CPU (phi.size()), grad_lapl_CPU(phi.size());
  phi_CPU = phi;
  grad_lapl_CPU = grad_lapl;
  gpu::host_vector<CudaRealType> vals_CPU(NumOrbitals), GL_CPU(4*row_stride);
  gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  HybridJobs_CPU = HybridJobs_GPU;
  gpu::host_vector<HybridData<CudaRealType> > HybridData_CPU(HybridData_GPU.size());
  HybridData_CPU = HybridData_GPU;
  rhats_CPU = rhats_GPU;
  // for (int iw=0; iw<newpos.size(); iw++)
  //   fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n",
  // 	       iw, rhats_CPU[3*iw+0], rhats_CPU[3*iw+1], rhats_CPU[3*iw+2]);
  for (int iw=0; iw<newpos.size(); iw++)
    if (false && HybridJobs_CPU[iw] == ATOMIC_POLY_JOB)
    {
      //if (HybridJobs_CPU[iw] != BSPLINE_3D_JOB && std::abs(rhats_CPU[3*iw+2]) < 1.0e-6) {
      int M = MakeTwoCopies.size();
      ComplexValueVector_t CPUzvals(M), CPUzlapl(M);
      ComplexGradVector_t CPUzgrad(M);
      ValueVector_t CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector_t CPUgrad(NumOrbitals);
      HybridData<CudaRealType> &d = HybridData_CPU[iw];
      AtomicOrbital<std::complex<double> > &atom = AtomicOrbitals[d.ion];
      atom.evaluate (newpos[iw], CPUzvals, CPUzgrad, CPUzlapl);
      int index=0;
      for (int i=0; i<StorageValueVector.size(); i++)
      {
        CPUvals[index] = CPUzvals[i].real();
        CPUlapl[index] = CPUzlapl[i].real();
        for (int j=0; j<OHMMS_DIM; j++)
          CPUgrad[index][j] = CPUzgrad[i][j].real();
        index++;
        if (MakeTwoCopies[i])
        {
          CPUvals[index] = CPUzvals[i].imag();
          CPUlapl[index] = CPUzlapl[i].imag();
          for (int j=0; j<OHMMS_DIM; j++)
            CPUgrad[index][j] = CPUzgrad[i][j].imag();
          index++;
        }
      }
      cudaMemcpy (&vals_CPU[0], phi_CPU[iw], NumOrbitals*sizeof(float),
                  cudaMemcpyDeviceToHost);
      cudaMemcpy (&GL_CPU[0], grad_lapl_CPU[iw], 4*row_stride*sizeof(float),
                  cudaMemcpyDeviceToHost);
      // fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
      // 	 iw, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
      double mindist = 1.0e5;
      for (int ion=0; ion<AtomicOrbitals.size(); ion++)
      {
        PosType disp = newpos[iw] - AtomicOrbitals[ion].Pos;
        PosType u = PrimLattice.toUnit(disp);
        PosType img;
        for (int i=0; i<OHMMS_DIM; i++)
          u[i] -= round(u[i]);
        disp = PrimLattice.toCart(u);
        double dist = std::sqrt(dot(disp,disp));
        if (dist < AtomicOrbitals[ion].CutoffRadius)
          mindist = dist;
      }
      if (std::fabs (mindist - d.dist) > 1.0e-3)
        fprintf (stderr, "CPU dist = %1.8f  GPU dist = %1.8f\n",
                 mindist, d.dist);
      fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f] dist = %10.6e\n",
               iw, rhats_CPU[3*iw+0], rhats_CPU[3*iw+1], rhats_CPU[3*iw+2],
               d.dist);
      for (int j=0; j<NumOrbitals; j++)
      {
        //	  if (isnan(vals_CPU[j])) {
        if (true || isnan(GL_CPU[0*row_stride+j]))
        {
          std::cerr << "iw = " << iw << std::endl;
          fprintf (stderr, "val[%2d]  = %10.5e %10.5e\n",
                   j, vals_CPU[j], CPUvals[j]);
          fprintf (stderr, "grad[%2d] = %10.5e %10.5e  %10.5e %10.5e  %10.5e %10.5e\n", j,
                   GL_CPU[0*row_stride+j], CPUgrad[j][0],
                   GL_CPU[1*row_stride+j], CPUgrad[j][1],
                   GL_CPU[2*row_stride+j], CPUgrad[j][2]);
          fprintf (stderr, "lapl[%2d] = %10.5e %10.5e\n",
                   j, GL_CPU[3*row_stride+j], CPUlapl[j]);
        }
      }
    }
    else if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
    {
      ComplexValueVector_t CPUzvals(NumOrbitals), CPUzlapl(NumOrbitals);
      ComplexGradVector_t CPUzgrad(NumOrbitals);
      ValueVector_t CPUvals(NumOrbitals), CPUlapl(NumOrbitals);
      GradVector_t CPUgrad(NumOrbitals);
      PosType ru(PrimLattice.toUnit(newpos[iw]));
      for (int i=0; i<3; i++)
        ru[i] -= std::floor(ru[i]);
      EinsplineMultiEval (MultiSpline, ru, CPUzvals, CPUzgrad,
                          StorageHessVector);
      for (int j=0; j<MakeTwoCopies.size(); j++)
      {
        CPUzgrad[j] = dot (PrimLattice.G, CPUzgrad[j]);
        CPUzlapl[j] = trace (StorageHessVector[j], GGt);
      }
      // Add e^-ikr phase to B-spline orbitals
      std::complex<double> eye(0.0, 1.0);
      for (int j=0; j<MakeTwoCopies.size(); j++)
      {
        std::complex<double> u = CPUzvals[j];
        TinyVector<std::complex<double>,OHMMS_DIM> gradu = CPUzgrad[j];
        std::complex<double> laplu = CPUzlapl[j];
        PosType k = kPoints[j];
        TinyVector<std::complex<double>,OHMMS_DIM> ck;
        for (int n=0; n<OHMMS_DIM; n++)	  ck[n] = k[n];
        double s,c;
        double phase = -dot(newpos[iw], k);
        sincos (phase, &s, &c);
        std::complex<double> e_mikr (c,s);
        CPUzvals[j]   = e_mikr*u;
        CPUzgrad[j]  = e_mikr*(-eye*u*ck + gradu);
        CPUzlapl[j]  = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
      }
      int index=0;
      for (int i=0; i<MakeTwoCopies.size(); i++)
      {
        CPUvals[index] = CPUzvals[i].real();
        CPUlapl[index] = CPUzlapl[i].real();
        for (int j=0; j<OHMMS_DIM; j++)
          CPUgrad[index][j] = CPUzgrad[i][j].real();
        index++;
        if (MakeTwoCopies[i])
        {
          CPUvals[index] = CPUzvals[i].imag();
          CPUlapl[index] = CPUzlapl[i].imag();
          for (int j=0; j<OHMMS_DIM; j++)
            CPUgrad[index][j] = CPUzgrad[i][j].imag();
          index++;
        }
      }
      cudaMemcpy (&vals_CPU[0], phi_CPU[iw], NumOrbitals*sizeof(float),
                  cudaMemcpyDeviceToHost);
      cudaMemcpy (&GL_CPU[0], grad_lapl_CPU[iw], 4*row_stride*sizeof(float),
                  cudaMemcpyDeviceToHost);
      // for (int i=0; i<4*row_stride; i++)
      //   fprintf (stderr, "%d %10.5e\n", i, GL_CPU[i]);
      static long int numgood=0, numbad=0;
      for (int j=0; j<NumOrbitals; j++)
      {
        double lap_ratio = GL_CPU[3*row_stride+j] /  CPUlapl[j];
        if (std::abs(GL_CPU[3*row_stride+j] - CPUlapl[j]) > 1.0e-4)
        {
          fprintf (stderr, "Error:  CPU laplacian = %1.8e  GPU = %1.8e\n",
                   CPUlapl[j],  GL_CPU[3*row_stride+j]);
          fprintf (stderr, "        CPU value     = %1.8e  GPU = %1.8e\n",
                   CPUvals[j],  vals_CPU[j]);
          fprintf (stderr, "u = %1.8f %1.8f %1.8f \n", ru[0], ru[1], ru[2]);
          std::cerr << "iw = " << iw << std::endl;
          numbad++;
        }
        else
          numgood++;
        if (numbad + numgood >= 100000)
        {
          double percent_bad = 100.0*(double)numbad / (double)(numbad + numgood);
          fprintf (stderr, "Percent bad = %1.8f\n", percent_bad);
          numbad = numgood = 0;
        }
        // if (true || isnan(GL_CPU[0*row_stride+j])) {
        //   // std::cerr << "r[" << iw << "] = " << newpos[iw] << std::endl;
        //   // std::cerr << "iw = " << iw << std::endl;
        //   fprintf (stderr, "\n3D      val[%2d]  = %10.5e %10.5e\n",
        // 	     j, vals_CPU[j], CPUvals[j]);
        //   fprintf (stderr, "3D GPU grad[%2d] = %10.5e %10.5e %10.5e\n", j,
        // 	     GL_CPU[0*row_stride+j],
        // 	     GL_CPU[1*row_stride+j],
        // 	     GL_CPU[2*row_stride+j]);
        //   fprintf (stderr, "3D CPU grad[%2d] = %10.5e %10.5e %10.5e\n", j,
        // 	     CPUgrad[j][0],
        // 	     CPUgrad[j][1],
        // 	     CPUgrad[j][2]);
        //   fprintf (stderr, "3D     lapl[%2d] = %10.5e %10.5e\n",
        // 	     j, GL_CPU[3*row_stride+j], CPUlapl[j]);
        // }
      }
    }
  gpu::host_vector<float> Ylm_CPU(Ylm_GPU.size());
  Ylm_CPU = Ylm_GPU;
  rhats_CPU = rhats_GPU;
  for (int i=0; i<rhats_CPU.size()/3; i++)
    fprintf (stderr, "rhat[%d] = [%10.6f %10.6f %10.6f]\n",
             i, rhats_CPU[3*i+0], rhats_CPU[3*i+1], rhats_CPU[3*i+2]);
  //    gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  //    HybridJobs_CPU = HybridJobs_GPU;
  //    gpu::host_vector<HybridDataFloat> HybridData_CPU(HybridData_GPU.size());
  //    HybridData_CPU = HybridData_GPU;
  std::cerr << "Before loop.\n";
  for (int i=0; i<newpos.size(); i++)
    if (HybridJobs_CPU[i] != BSPLINE_3D_JOB)
    {
      std::cerr << "Inside if.\n";
      PosType rhat(rhats_CPU[3*i+0], rhats_CPU[3*i+1], rhats_CPU[3*i+2]);
      AtomicOrbital<std::complex<double> > &atom = AtomicOrbitals[HybridData_CPU[i].ion];
      int numlm = (atom.lMax+1)*(atom.lMax+1);
      std::vector<double> Ylm(numlm), dYlm_dtheta(numlm), dYlm_dphi(numlm);
      atom.CalcYlm (rhat, Ylm, dYlm_dtheta, dYlm_dphi);
      for (int lm=0; lm < numlm; lm++)
      {
        fprintf (stderr, "lm=%3d  Ylm_CPU=%8.5f  Ylm_GPU=%8.5f\n",
                 lm, Ylm[lm], Ylm_CPU[3*i*Ylm_BS+lm]);
      }
    }
  fprintf (stderr, " N  img      dist    ion    lMax\n");
  for (int i=0; i<HybridData_CPU.size(); i++)
  {
    HybridData<CudaRealType> &d = HybridData_CPU[i];
    fprintf (stderr, " %d %2.0f %2.0f %2.0f  %8.5f  %d %d\n",
             i, d.img[0], d.img[1], d.img[2], d.dist, d.ion, d.lMax);
  }
#endif
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
    gpu::device_vector<CudaComplexType*> &phi,
    gpu::device_vector<CudaComplexType*> &grad_lapl,
    int row_stride)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate \n"
              << "(std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos, \n"
              << " gpu::device_vector<CudaComplexType*> &phi,\n"
              << " gpu::device_vector<CudaComplexType*> &grad_lapl,\n"
              << " int row_stride)\n"
              << "not yet implemented.\n";
  abort();
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  int N = pos.size();
  if (cudapos.size() < N)
  {
    resize_cuda(N);
    hostPos.resize(N);
    cudapos.resize(N);
  }
  for (int iw=0; iw < N; iw++)
    hostPos[iw] = pos[iw];
  cudapos = hostPos;

  MakeHybridJobList<CudaRealType> ((CudaRealType*)cudapos.data(), N, IonPos_GPU.data(),
                                   PolyRadii_GPU.data(), CutoffRadii_GPU.data(),
                                   AtomicOrbitals.size(), L_cuda.data(), Linv_cuda.data(),
                                   HybridJobs_GPU.data(), rhats_GPU.data(),
                                   HybridData_GPU.data());

  CalcYlmComplexCuda<CudaRealType> (rhats_GPU.data(), HybridJobs_GPU.data(),
                                    Ylm_ptr_GPU.data(), dYlm_dtheta_ptr_GPU.data(),
                                    dYlm_dphi_ptr_GPU.data(), lMax, pos.size());

  evaluate3DSplineComplexToReal<CudaRealType>
  (HybridJobs_GPU.data(), (CudaRealType*)cudapos.data(),
   (CudaRealType*)CudakPoints.data(),CudaMakeTwoCopies.data(),
   CudaMultiSpline, Linv_cuda.data(),
   phi.data(), CudaMakeTwoCopies.size(), pos.size());

  evaluateHybridSplineComplexToReal<CudaRealType> //NLPP
  (HybridJobs_GPU.data(),
   Ylm_ptr_GPU.data(),
   AtomicOrbitals_GPU.data(), HybridData_GPU.data(),
   (CudaRealType*)CudakPoints_reduced.data(),
   CudaMakeTwoCopies.data(), (CudaRealType**)phi.data(),
   CudaMakeTwoCopies.size(), pos.size(), lMax);

  // gpu::host_vector<HybridJobType> HybridJobs_CPU(HybridJobs_GPU.size());
  // HybridJobs_CPU = HybridJobs_GPU;
  // int M = CudaMakeTwoCopies.size();
  // ComplexValueVector_t CPUzvals(M);
  // for (int iw=0; iw<pos.size(); iw++)
  //   if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  // 	cerr << "Error:  used BSPLINE_3D for PP eval.  Walker = "
  // 	     << iw << "\n";
  // std::cerr << "pos.size() = " << pos.size() << std::endl;
  // for (int iw=0; iw<pos.size(); iw++) {
  //   // if (HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   // 	num3D++;
  //   // else
  //   // 	numAtomic++;
  //   bool atomic=false;
  //   for (int ion=0; ion<AtomicOrbitals.size(); ion++)
  //   	if (AtomicOrbitals[ion].evaluate (pos[iw], CPUzvals)) {
  //   	  atomic = true;
  //   	  break;
  //   	}
  //   // if (atomic)
  //   // 	numAtomic++;
  //   // else
  //   // 	num3D++;
  //   if (atomic && HybridJobs_CPU[iw] == BSPLINE_3D_JOB)
  //   	cerr << "Error!  Used BSPLINE_3D when I should have used atomic.\n";
  //   else if (!atomic && HybridJobs_CPU[iw] != BSPLINE_3D_JOB)
  //   	cerr << "Error!  Used atomic when I should have used 3D.\n";
  // }
  //////////////////////////
  // This must be tested! //
  //////////////////////////
}

template<> void
EinsplineSetHybrid<std::complex<double> >::evaluate
(std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "EinsplineSetHybrid<std::complex<double> >::evaluate \n"
              << "(std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)\n"
              << "not yet implemented.\n";
}

template<> void
EinsplineSetExtended<double>::initGPU()
{
  app_log() << "Copying einspline orbitals to GPU.\n";
  create_multi_UBspline_3d_cuda
  (MultiSpline, CudaMultiSpline);
  app_log() << "Successful copy.\n";
  // Destroy original CPU spline
  // HACK HACK HACK
  //destroy_Bspline (MultiSpline);
  gpu::host_vector<CudaRealType> L_host(9), Linv_host(9);
  Linv_cuda.resize(9);
  L_cuda.resize(9);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      L_host[i*3+j]    = PrimLattice.R(i,j);
      Linv_host[i*3+j] = PrimLattice.G(i,j);
    }
  L_cuda    = L_host;
  Linv_cuda = Linv_host;
}

template<> void
EinsplineSetExtended<std::complex<double> >::initGPU()
{
  app_log() << "Copying einspline orbitals to GPU.\n";
  create_multi_UBspline_3d_cuda
  (MultiSpline, CudaMultiSpline);
  app_log() << "Successful copy.\n";
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
      L_host[i*3+j]    = PrimLattice.R(i,j);
      Linv_host[i*3+j] = PrimLattice.G(i,j);
    }
  L_cuda    = L_host;
  Linv_cuda = Linv_host;
}



template<> void
EinsplineSetHybrid<double>::initGPU()
{
  EinsplineSetExtended<double>::initGPU();
  // Setup B-spline Acuda matrix in constant memory
  init_atomic_cuda();
  gpu::host_vector<AtomicOrbitalCuda<CudaRealType> > AtomicOrbitals_CPU;
  const int BS=16;
  NumOrbitals = getOrbitalSetSize();
  // Bump up the stride to be a multiple of 512-bit bus width
  int lm_stride = ((NumOrbitals+BS-1)/BS)*BS;
  AtomicSplineCoefs_GPU.resize(AtomicOrbitals.size());
  AtomicPolyCoefs_GPU.resize(AtomicOrbitals.size());
  std::vector<CudaRealType> IonPos_CPU, CutoffRadii_CPU, PolyRadii_CPU;
  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
  {
    app_log() << "Copying real atomic orbitals for ion " << iat << " to GPU memory.\n";
    AtomicOrbital<double> &atom = AtomicOrbitals[iat];
    for (int i=0; i<OHMMS_DIM; i++)
      IonPos_CPU.push_back(atom.Pos[i]);
    CutoffRadii_CPU.push_back(atom.CutoffRadius);
    PolyRadii_CPU.push_back(atom.PolyRadius);
    AtomicOrbitalCuda<CudaRealType> atom_cuda;
    atom_cuda.lMax = atom.lMax;
    int numlm = (atom.lMax+1)*(atom.lMax+1);
    atom_cuda.lm_stride = lm_stride;
    atom_cuda.spline_stride = numlm * lm_stride;
    AtomicOrbital<double>::SplineType &cpu_spline =
      *atom.get_radial_spline();
    atom_cuda.spline_dr_inv = cpu_spline.x_grid.delta_inv;
    int Ngrid = cpu_spline.x_grid.num;
    int spline_size = 2*atom_cuda.spline_stride * (Ngrid+2);
    gpu::host_vector<CudaRealType> spline_coefs(spline_size);
    AtomicSplineCoefs_GPU[iat].resize(spline_size);
    atom_cuda.spline_coefs = AtomicSplineCoefs_GPU[iat].data();
    // Reorder and copy splines to GPU memory
    for (int igrid=0; igrid<Ngrid; igrid++)
      for (int lm=0; lm<numlm; lm++)
        for (int orb=0; orb<NumOrbitals; orb++)
        {
          // Convert splines to Hermite spline form, because
          // B-spline form in single precision with dense grids
          // leads to large truncation error
          double u = 1.0/6.0*
                     (1.0*cpu_spline.coefs[(igrid+0)*cpu_spline.x_stride + orb*numlm +lm] +
                      4.0*cpu_spline.coefs[(igrid+1)*cpu_spline.x_stride + orb*numlm +lm] +
                      1.0*cpu_spline.coefs[(igrid+2)*cpu_spline.x_stride + orb*numlm +lm]);
          double d2u = cpu_spline.x_grid.delta_inv * cpu_spline.x_grid.delta_inv *
                       (1.0*cpu_spline.coefs[(igrid+0)*cpu_spline.x_stride + orb*numlm +lm] +
                        -2.0*cpu_spline.coefs[(igrid+1)*cpu_spline.x_stride + orb*numlm +lm] +
                        1.0*cpu_spline.coefs[(igrid+2)*cpu_spline.x_stride + orb*numlm +lm]);
          spline_coefs[(2*igrid+0)*atom_cuda.spline_stride +
                       lm   *atom_cuda.lm_stride + orb] = u;
          spline_coefs[(2*igrid+1)*atom_cuda.spline_stride +
                       lm   *atom_cuda.lm_stride + orb] = d2u;
        }
    AtomicSplineCoefs_GPU[iat] = spline_coefs;
    atom_cuda.poly_stride = numlm*atom_cuda.lm_stride;
    atom_cuda.poly_order = atom.PolyOrder;
    int poly_size = (atom.PolyOrder+1)*atom_cuda.poly_stride;
    gpu::host_vector<CudaRealType> poly_coefs(poly_size);
    AtomicPolyCoefs_GPU[iat].resize(poly_size);
    atom_cuda.poly_coefs = AtomicPolyCoefs_GPU[iat].data();
    for (int lm=0; lm<numlm; lm++)
      for (int n=0; n<atom.PolyOrder; n++)
        for (int orb=0; orb<NumOrbitals; orb++)
          poly_coefs[n*atom_cuda.poly_stride + lm*atom_cuda.lm_stride + orb]
            = atom.get_poly_coefs()(n,orb,lm);
    AtomicPolyCoefs_GPU[iat] = poly_coefs;
    AtomicOrbitals_CPU.push_back(atom_cuda);
  }
  AtomicOrbitals_GPU = AtomicOrbitals_CPU;
  IonPos_GPU      = IonPos_CPU;
  CutoffRadii_GPU = CutoffRadii_CPU;
  PolyRadii_GPU   = PolyRadii_CPU;
}

template<> void
EinsplineSetHybrid<std::complex<double> >::initGPU()
{
  EinsplineSetExtended<std::complex<double> >::initGPU();
  // Setup B-spline Acuda matrix in constant memory
  init_atomic_cuda();
  gpu::host_vector<AtomicOrbitalCuda<CudaRealType> > AtomicOrbitals_CPU;
  const int BS=16;
  NumOrbitals = getOrbitalSetSize();
  // Bump up the stride to be a multiple of 512-bit bus width
  int lm_stride = ((2*NumOrbitals+BS-1)/BS)*BS;
  AtomicSplineCoefs_GPU.resize(AtomicOrbitals.size());
  AtomicPolyCoefs_GPU.resize(AtomicOrbitals.size());
  std::vector<CudaRealType> IonPos_CPU, CutoffRadii_CPU, PolyRadii_CPU;
  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
  {
    app_log() << "Copying atomic orbitals for ion " << iat << " to GPU memory.\n";
    AtomicOrbital<std::complex<double> > &atom = AtomicOrbitals[iat];
    for (int i=0; i<OHMMS_DIM; i++)
      IonPos_CPU.push_back(atom.Pos[i]);
    CutoffRadii_CPU.push_back(atom.CutoffRadius);
    PolyRadii_CPU.push_back(atom.PolyRadius);
    AtomicOrbitalCuda<CudaRealType> atom_cuda;
    atom_cuda.lMax = atom.lMax;
    int numlm = (atom.lMax+1)*(atom.lMax+1);
    atom_cuda.lm_stride = lm_stride;
    atom_cuda.spline_stride = numlm * lm_stride;
    AtomicOrbital<std::complex<double> >::SplineType &cpu_spline =
      *atom.get_radial_spline();
    atom_cuda.spline_dr_inv = cpu_spline.x_grid.delta_inv;
    int Ngrid = cpu_spline.x_grid.num;
    int spline_size = 2*atom_cuda.spline_stride * (Ngrid+2);
    gpu::host_vector<CudaRealType> spline_coefs(spline_size);
    AtomicSplineCoefs_GPU[iat].resize(spline_size);
    atom_cuda.spline_coefs = AtomicSplineCoefs_GPU[iat].data();
    // Reorder and copy splines to GPU memory
    for (int igrid=0; igrid<Ngrid; igrid++)
      for (int lm=0; lm<numlm; lm++)
        for (int orb=0; orb<NumOrbitals; orb++)
        {
          // Convert splines to Hermite spline form, because
          // B-spline form in single precision with dense grids
          // leads to large truncation error
          std::complex<double> u = 1.0/6.0*
                              (1.0*cpu_spline.coefs[(igrid+0)*cpu_spline.x_stride + orb*numlm +lm] +
                               4.0*cpu_spline.coefs[(igrid+1)*cpu_spline.x_stride + orb*numlm +lm] +
                               1.0*cpu_spline.coefs[(igrid+2)*cpu_spline.x_stride + orb*numlm +lm]);
          std::complex<double> d2u = cpu_spline.x_grid.delta_inv * cpu_spline.x_grid.delta_inv *
                                (1.0*cpu_spline.coefs[(igrid+0)*cpu_spline.x_stride + orb*numlm +lm] +
                                 -2.0*cpu_spline.coefs[(igrid+1)*cpu_spline.x_stride + orb*numlm +lm] +
                                 1.0*cpu_spline.coefs[(igrid+2)*cpu_spline.x_stride + orb*numlm +lm]);
          spline_coefs[((2*igrid+0)*atom_cuda.spline_stride + lm*atom_cuda.lm_stride + 2*orb)+0] =
            u.real();
          spline_coefs[((2*igrid+0)*atom_cuda.spline_stride + lm*atom_cuda.lm_stride + 2*orb)+1] =
            u.imag();
          spline_coefs[((2*igrid+1)*atom_cuda.spline_stride + lm*atom_cuda.lm_stride + 2*orb)+0] =
            d2u.real();
          spline_coefs[((2*igrid+1)*atom_cuda.spline_stride + lm*atom_cuda.lm_stride + 2*orb)+1] =
            d2u.imag();
        }
    AtomicSplineCoefs_GPU[iat] = spline_coefs;
    atom_cuda.poly_stride = numlm*atom_cuda.lm_stride;
    atom_cuda.poly_order = atom.PolyOrder;
    int poly_size = (atom.PolyOrder+1)*atom_cuda.poly_stride;
    gpu::host_vector<CudaRealType> poly_coefs(poly_size);
    AtomicPolyCoefs_GPU[iat].resize(poly_size);
    atom_cuda.poly_coefs = &AtomicPolyCoefs_GPU[iat].data()[0];
    for (int lm=0; lm<numlm; lm++)
      for (int n=0; n<atom.PolyOrder; n++)
        for (int orb=0; orb<NumOrbitals; orb++)
        {
          poly_coefs
          [n*atom_cuda.poly_stride + lm*atom_cuda.lm_stride + 2*orb+0] =
            atom.get_poly_coefs()(n,orb,lm).real();
          poly_coefs
          [n*atom_cuda.poly_stride + lm*atom_cuda.lm_stride + 2*orb+1] =
            atom.get_poly_coefs()(n,orb,lm).imag();
        }
    AtomicPolyCoefs_GPU[iat] = poly_coefs;
    AtomicOrbitals_CPU.push_back(atom_cuda);
  }
  AtomicOrbitals_GPU = AtomicOrbitals_CPU;
  IonPos_GPU      = IonPos_CPU;
  CutoffRadii_GPU = CutoffRadii_CPU;
  PolyRadii_GPU   = PolyRadii_CPU;
}



}
