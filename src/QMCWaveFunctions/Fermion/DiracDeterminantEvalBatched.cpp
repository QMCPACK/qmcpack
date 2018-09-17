
#include "QMCWaveFunctions/Fermion/DiracDeterminantEvalBatched.h"

namespace qmcplusplus
{

  DiracDeterminantEval<Batching::BATCHED>::DiracDeterminantEval() :
    UpdateJobList_d("DiracDeterminantBase::UpdateJobList_d"),
    srcList_d("DiracDeterminantBase::srcList_d"),
    destList_d("DiracDeterminantBase::destList_d"),
    AList_d("DiracDeterminantBase::AList_d"),
    AinvList_d("DiracDeterminantBase::AinvList_d"),
    newRowList_d("DiracDeterminantBase::newRowList_d"),
    AinvDeltaList_d("DiracDeterminantBase::AinvDeltaList_d"),
    AinvColkList_d("DiracDeterminantBase::AinvColkList_d"),
    gradLaplList_d("DiracDeterminantBase::gradLaplList_d"),
    newGradLaplList_d("DiracDeterminantBase::newGradLaplList_d"),
    AWorkList_d("DiracDeterminantBase::AWorkList_d"),
    AinvWorkList_d("DiracDeterminantBase::AinvWorkList_d"),
    PivotArray_d("DiracDeterminantBase::PivotArray_d"),
    infoArray_d("DiracDeterminantBase::infoArray_d"),
    GLList_d("DiracDeterminantBase::GLList_d"),
    ratio_d("DiracDeterminantBase::ratio_d"),
    gradLapl_d("DiracDeterminantBase::gradLapl_d"),
    iatList_d("DiracDeterminantBase::iatList_d"),
    NLrowBuffer_d("DiracDeterminantBase::NLrowBuffer_d"),
    SplineRowList_d("DiracDeterminantBase::SplineRowList_d"),
    RatioRowList_d("DiracDeterminantBase::RatioRowList_d"),
    NLposBuffer_d("DiracDeterminantBase::NLposBuffer_d"),
    NLAinvList_d("DiracDeterminantBase::NLAinvList_d"),
    NLnumRatioList_d("DiracDeterminantBase::NLnumRatioList_d"),
    NLelecList_d("DiracDeterminantBase::NLelecList_d"),
    NLratioList_d("DiracDeterminantBase::NLratioList_d")
  {
     for(int i = 0; i < 2; ++i)
      NLratios_d[i] = gpu::device_vector<CudaValueType>("DiracDeterminantBase::NLratios_d");
  }
  
  void DiracDeterminantEval<Batching::BATCHED>
  ::reserve (PointerPool<gpu::device_vector<QMCT::CudaValueType>> &pool,
	     SPOSetBatched& Phi, int num_particles, int num_orbitals)
{
    RowStride = ((num_orbitals + 31)/32) * 32;
    AOffset           = pool.reserve((size_t)    num_particles * RowStride);
    AinvOffset        = pool.reserve((size_t)    num_particles * RowStride);
    gradLaplOffset    = pool.reserve((size_t)4 * num_particles * RowStride);
    newRowOffset      = pool.reserve((size_t)1            * RowStride);
    AinvDeltaOffset   = pool.reserve((size_t)1            * RowStride);
    AinvColkOffset    = pool.reserve((size_t)1            * RowStride);
    newGradLaplOffset = pool.reserve((size_t)4            * RowStride);
    if (typeid(CudaRealType) == typeid(float))
    {
      AWorkOffset       = pool.reserve((size_t)2 * num_particles * RowStride);
      AinvWorkOffset    = pool.reserve((size_t)2 * num_particles * RowStride);
    }
    else if (typeid(CudaRealType) == typeid(double))
    {
      AWorkOffset       = pool.reserve((size_t)    num_particles * RowStride);
      AinvWorkOffset    = 0;                  // not needed for inversion
    }
    Phi.reserve(pool);
}

void
DiracDeterminantEval<Batching::BATCHED>::addLog (MCWalkerConfiguration &W, std::vector<QMCT::RealType> &logPsi)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(),NumOrbitals);
  std::vector<QMCT::PosType> R(walkers.size());
  // Fill in the A matrix row by row
  for (int iat=FirstIndex; iat<LastIndex; iat++)
  {
    int off = (iat-FirstIndex)*RowStride;
    for (int iw=0; iw<walkers.size(); iw++)
    {
      Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
      newRowList[iw]    =  &(data.data()[AOffset+off]);
      gradLaplList[iw]  =  &(data.data()[gradLaplOffset+4*off]);
      R[iw] = walkers[iw]->R[iat];
    }
    newRowList_d = newRowList;
    gradLaplList_d = gradLaplList;
    Phi->evaluate (walkers, R, newRowList_d, gradLaplList_d, RowStride);
  }
  // Now, compute determinant
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    gpu::host_vector<CudaValueType> host_data;
    Vector<CudaValueType> A(NumPtcls*NumOrbitals);
    host_data = data;
    for (int i=0; i<NumPtcls; i++)
      for (int j=0; j<NumOrbitals; j++)
        A[i*NumOrbitals+j] = host_data[AOffset+i*RowStride+j];
    logPsi[iw] += std::log(std::abs(Invert(A.data(), NumPtcls, NumOrbitals)));
    int N = NumPtcls;
    bool passed = true;
    for (int i=0; i<NumPtcls; i++)
    {
      for (int j=0; j<NumOrbitals; j++)
        host_data[AinvOffset+i*RowStride+j] = A[i*NumOrbitals + j];
      for (int j=NumOrbitals; j<RowStride; j++)
        host_data[AinvOffset+i*RowStride+j] = 0.0;
    }
    data = host_data;
    // for (int i=0; i<N; i++)
    // 	for (int j=0; j<N; j++) {
    // 	  double val = 0.0;
    // 	  for (int k=0; k<N; k++) {
    // 	    double aval = host_data[AOffset+i*N+k];
    // 	    double ainv = host_data[AinvOffset+k*N+j];
    // 	    val += aval * ainv;
    // 	  }
    // 	  if (i == j) {
    // 	    if (std::abs(val - 1.0) > 1.0e-2) {
    // 	      app_error() << "Error in inverse, (i,j) = " << i << ", " << j << ".\n";
    // 	      passed = false;
    // 	    }
    // 	  }
    // 	  else
    // 	    if (std::abs(val) > 1.0e-2) {
    // 	      app_error() << "Error in inverse, (i,j) = " << i << ", " << j << ".\n";
    // 	      passed = false;
    // 	    }
    // 	}
    // if (!passed)
    // 	app_log() << (passed ? "Passed " : "Failed " ) << "inverse test.\n";
  }
}

void
DiracDeterminantEval<Batching::BATCHED>::calcGradient(MCWalkerConfiguration &W, int iat,
                                   std::vector<GradType> &grad)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(), NumOrbitals);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]        =  &(data.data()[AinvOffset]);
    GLList[iw]    =  &(data.data()[gradLaplOffset]);
  }
  AinvList_d.asyncCopy(AinvList);
  GLList_d.asyncCopy(GLList);
  calc_gradient (AinvList_d.data(), GLList_d.data(),
                 ratio_d.data(), NumOrbitals, RowStride,
                 iat-FirstIndex, walkers.size());
  gpu::streamsSynchronize();
  ratio_host.asyncCopy(ratio_d);
  cudaEventRecord(gpu::gradientSyncDiracEvent, gpu::memoryStream);
}

void
DiracDeterminantEval<Batching::BATCHED>::addGradient(MCWalkerConfiguration &W, int iat,
                                  std::vector<GradType> &grad)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncDiracEvent);
  for (int iw=0; iw<walkers.size(); iw++)
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      grad[iw][dim] += ratio_host[3*iw+dim];
    }
#ifdef CUDA_DEBUG3
  if (NumOrbitals == 31)
  {
    gpu::host_vector<CudaRealType> host_data;
    std::vector<CudaRealType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw=0; iw<walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb=0; iorb<NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset+RowStride*iorb+iat-FirstIndex] *
                          host_data[gradLaplOffset + 4*RowStride*(iat-FirstIndex) + iorb + RowStride];
      }
      fprintf (stderr, "CPU grad = %10.6e   GPU grad = %10.6e\n",
               cpu_ratios[iw], grad[iw][1]);
    }
  }
#endif
}

void DiracDeterminantEval<Batching::BATCHED>::ratio (MCWalkerConfiguration &W,
                                  int iat,
                                  std::vector<ValueType> &psi_ratios)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(), NumOrbitals);
  // First evaluate orbitals
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]      =  &(data.data()[AinvOffset]);
    newRowList[iw]    =  &(data.data()[newRowOffset]);
  }
  newRowList_d = newRowList;
  AinvList_d   = AinvList;
  Phi->evaluate (walkers, W.Rnew, newRowList_d);
  // Now evaluate ratios
  determinant_ratios_cuda
  (AinvList_d.data(), newRowList_d.data(), ratio_d.data(),
   NumPtcls, RowStride, iat-FirstIndex, walkers.size());
  // Copy back to host
  ratio_host = ratio_d;
  for (int iw=0; iw<psi_ratios.size(); iw++)
    psi_ratios[iw] *= ratio_host[iw];
}

DiracDeterminantEval<Batching::BATCHED>::ValueType
DiracDeterminantEval<Batching::BATCHED>::ratio(ParticleSet& P, int iat)
{
  APP_ABORT("This should not be called on DiracDeterminantEval<Batching::BATCHED>\n");
}


void DiracDeterminantEval<Batching::BATCHED>::ratio (MCWalkerConfiguration &W, int iat,
                                  std::vector<ValueType> &psi_ratios,
                                  std::vector<GradType>  &grad)
{
}


// The gradient is (\nabla psi(Rnew))/psi(Rnew)
// The laplaican is (\nabla^2 psi(Rnew))/psi(Rew)
void DiracDeterminantEval<Batching::BATCHED>::ratio (MCWalkerConfiguration &W, int iat,
                                  std::vector<ValueType> &psi_ratios,
                                  std::vector<GradType>  &grad,
                                  std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(),NumOrbitals);
  //    if (iat-FirstIndex == 0) {
  // First evaluate orbitals
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]        =  &(data.data()[AinvOffset]);
    newRowList[iw]      =  &(data.data()[newRowOffset]);
    newGradLaplList[iw] =  &(data.data()[newGradLaplOffset]);
    if (((unsigned long)AinvList[iw] %64) ||
        ((unsigned long)newRowList[iw] % 64) ||
        ((unsigned long)newGradLaplList[iw] %64))
      app_log() << "**** CUDA misalignment!!!! ***\n";
  }
  newRowList_d.asyncCopy(newRowList);
  newGradLaplList_d.asyncCopy(newGradLaplList);
  AinvList_d.asyncCopy(AinvList);
  //    }
  Phi->evaluate (walkers, W.Rnew, newRowList_d, newGradLaplList_d, RowStride);
#ifdef CUDA_DEBUG2
  Vector<ValueType> testPhi(NumOrbitals), testLapl(NumOrbitals);
  Vector<GradType> testGrad(NumOrbitals);
  ParticleSet P;
  P.R.resize(NumPtcls);
  gpu::host_vector<CudaValueType> host_vec;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    host_vec = walkers[iw]->cuda_DataSet;
    P.R[iat-FirstIndex] = W.Rnew[iw];
    //      Phi->evaluate(P, iat-FirstIndex, testPhi);
    Phi->evaluate(P, iat-FirstIndex, testPhi, testGrad, testLapl);
    for (int iorb=0; iorb<NumOrbitals; iorb++)
    {
      //if (std::abs(host_vec[newRowOffset+iorb]-testPhi[iorb]) > 1.0e-6)
      //   fprintf (stderr, "CUDA = %1.8e    CPU = %1.8e\n",
      // 	   host_vec[newRowOffset+iorb], testPhi[iorb]);
      fprintf (stderr, "CUDA = %1.8e    CPU = %1.8e\n",
               host_vec[newGradLaplOffset+2*NumOrbitals+iorb], testGrad[iorb][2]);
    }
  }
#endif
  determinant_ratios_grad_lapl_cuda
  (AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(),
   ratio_d.data(), NumPtcls, RowStride, iat-FirstIndex, walkers.size());
  // Copy back to host
  ratio_host = ratio_d;
#ifdef CUDA_DEBUG
  // Now, check against CPU
  gpu::host_vector<CudaValueType> host_data;
  std::vector<CudaValueType> cpu_ratios(walkers.size(), 0.0f);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    host_data = walkers[iw]->cuda_DataSet;
    for (int iorb=0; iorb<NumOrbitals; iorb++)
    {
      cpu_ratios[iw] += host_data[AinvOffset+RowStride*iorb+iat-FirstIndex] *
                        host_data[newRowOffset + iorb];
    }
    fprintf (stderr, "CPU ratio = %10.6e   GPU ratio = %10.6e\n",
             cpu_ratios[iw], ratio_host[5*iw+0]);
  }
#endif
  // Calculate ratio, gradient and laplacian
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= ratio_host[5*iw+0];
    GradType g(ratio_host[5*iw+1],
               ratio_host[5*iw+2],
               ratio_host[5*iw+3]);
    grad[iw] += g;
    lapl[iw] += ratio_host[5*iw+0];
    // grad[iw] += g / ratio_host[5*iw+0];
    // lapl[iw] += (ratio_host[5*iw+4] - dot(g,g)) / ratio_host[5*iw+0];
  }
#ifdef CUDA_DEBUG
  if (NumOrbitals == 31)
  {
    gpu::host_vector<CudaValueType> host_data;
    std::vector<CudaValueType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw=0; iw<walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb=0; iorb<NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset+RowStride*iorb+iat-FirstIndex] *
                          host_data[newGradLaplOffset + iorb] / ratio_host[5*iw+0];
      }
      fprintf (stderr, "ratio CPU grad = %10.6e   GPU grad = %10.6e\n",
               cpu_ratios[iw], grad[iw][0]);
    }
  }
#endif
}
void DiracDeterminantEval<Batching::BATCHED>::calcRatio (MCWalkerConfiguration &W, int iat,
                                      std::vector<ValueType> &psi_ratios,
                                      std::vector<GradType>  &grad,
                                      std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(), NumOrbitals);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]        =  &(data.data()[AinvOffset]);
    newRowList[iw]      =  &(data.data()[newRowOffset]);
    newGradLaplList[iw] =  &(data.data()[newGradLaplOffset]);
    if (((unsigned long)AinvList[iw] %64) ||
        ((unsigned long)newRowList[iw] % 64) ||
        ((unsigned long)newGradLaplList[iw] %64))
      app_log() << "**** CUDA misalignment!!!! ***\n";
  }
  newRowList_d.asyncCopy(newRowList);
  newGradLaplList_d.asyncCopy(newGradLaplList);
  AinvList_d.asyncCopy(AinvList);
  Phi->evaluate (walkers, W.Rnew, newRowList_d, newGradLaplList_d, RowStride);
#ifdef CUDA_DEBUG2
  Vector<ValueType> testPhi(NumOrbitals), testLapl(NumOrbitals);
  Vector<GradType> testGrad(NumOrbitals);
  ParticleSet P;
  P.R.resize(NumPtcls);
  gpu::host_vector<CudaValueType> host_vec;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    host_vec = walkers[iw]->cuda_DataSet;
    P.R[iat-FirstIndex] = W.Rnew[iw];
    Phi->evaluate(P, iat-FirstIndex, testPhi, testGrad, testLapl);
    for (int iorb=0; iorb<NumOrbitals; iorb++)
    {
      fprintf (stderr, "CUDA = %1.8e    CPU = %1.8e\n",
               host_vec[newGradLaplOffset+2*NumOrbitals+iorb], testGrad[iorb][2]);
    }
  }
#endif
  determinant_ratios_grad_lapl_cuda
  (AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(),
   ratio_d.data(), NumPtcls, RowStride, iat-FirstIndex, walkers.size());
  gpu::streamsSynchronize();
  // Copy back to host
  ratio_host.asyncCopy(ratio_d);
  cudaEventRecord(gpu::ratioSyncDiracEvent, gpu::memoryStream);
#ifdef CUDA_DEBUG
  // Now, check against CPU
  gpu::host_vector<CudaValueType> host_data;
  std::vector<CudaValueType> cpu_ratios(walkers.size(), 0.0f);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    host_data = walkers[iw]->cuda_DataSet;
    for (int iorb=0; iorb<NumOrbitals; iorb++)
    {
      cpu_ratios[iw] += host_data[AinvOffset+RowStride*iorb+iat-FirstIndex] *
                        host_data[newRowOffset + iorb];
    }
    fprintf (stderr, "CPU ratio = %10.6e   GPU ratio = %10.6e\n",
             cpu_ratios[iw], ratio_host[5*iw+0]);
  }
#endif
}

void DiracDeterminantEval<Batching::BATCHED>::addRatio (MCWalkerConfiguration &W, int iat,
                                     std::vector<ValueType> &psi_ratios,
                                     std::vector<GradType>  &grad,
                                     std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::ratioSyncDiracEvent);
  // Calculate ratio, gradient and laplacian
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= ratio_host[5*iw+0];
    GradType g(ratio_host[5*iw+1],
               ratio_host[5*iw+2],
               ratio_host[5*iw+3]);
    grad[iw] += g;
    lapl[iw] += ratio_host[5*iw+0];
  }
#ifdef CUDA_DEBUG
  if (NumOrbitals == 31)
  {
    gpu::host_vector<CudaValueType> host_data;
    std::vector<CudaValueType> cpu_ratios(walkers.size(), 0.0f);
    for (int iw=0; iw<walkers.size(); iw++)
    {
      host_data = walkers[iw]->cuda_DataSet;
      for (int iorb=0; iorb<NumOrbitals; iorb++)
      {
        cpu_ratios[iw] += host_data[AinvOffset+RowStride*iorb+iat-FirstIndex] *
                          host_data[newGradLaplOffset + iorb] / ratio_host[5*iw+0];
      }
      fprintf (stderr, "ratio CPU grad = %10.6e   GPU grad = %10.6e\n",
               cpu_ratios[iw], grad[iw][0]);
    }
  }
#endif
}




// The gradient is (\nabla psi(Rnew))/psi(Rnew)
// The laplaican is (\nabla^2 psi(Rnew))/psi(Rew)
void DiracDeterminantEval<Batching::BATCHED>::ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iat_list,
                                  std::vector<PosType> &rNew,
                                  std::vector<ValueType> &psi_ratios,
                                  std::vector<GradType>  &grad,
                                  std::vector<ValueType> &lapl)
{
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(), NumOrbitals);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]        =  &(data.data()[AinvOffset]);
    newRowList[iw]      =  &(data.data()[newRowOffset]);
    newGradLaplList[iw] =  &(data.data()[newGradLaplOffset]);
    if (((unsigned long)AinvList[iw] %64) ||
        ((unsigned long)newRowList[iw] % 64) ||
        ((unsigned long)newGradLaplList[iw] %64))
      app_log() << "**** CUDA misalignment!!!! ***\n";
    iatList[iw] = iat_list[iw] - FirstIndex;
  }
  newRowList_d = newRowList;
  newGradLaplList_d = newGradLaplList;
  AinvList_d   = AinvList;
  iatList_d    = iatList;
  Phi->evaluate (walkers, rNew, newRowList_d, newGradLaplList_d, RowStride);
  determinant_ratios_grad_lapl_cuda
  (AinvList_d.data(), newRowList_d.data(), newGradLaplList_d.data(),
   ratio_d.data(), NumPtcls, RowStride, iatList_d.data(), walkers.size());
  // Copy back to host
  ratio_host = ratio_d;
  // Calculate ratio, gradient and laplacian
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= ratio_host[5*iw+0];
    GradType g(ratio_host[5*iw+1],
               ratio_host[5*iw+2],
               ratio_host[5*iw+3]);
    grad[iw] += g;
    lapl[iw] += ratio_host[5*iw+0];
  }
}



void
DiracDeterminantEval<Batching::BATCHED>::gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                                ValueMatrix_t &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (AList.size() < walkers.size())
    resizeLists(walkers.size(), NumOrbitals);
  // First evaluate orbitals
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t::cuda_Buffer_t& data = walkers[iw]->cuda_DataSet;
    AinvList[iw]        =  &(data.data()[AinvOffset]);
    gradLaplList[iw]    =  &(data.data()[gradLaplOffset]);
    newGradLaplList[iw] =  &(gradLapl_d.data()[4*RowStride*iw]);
  }
  AinvList_d      = AinvList;
  gradLaplList_d = gradLaplList;
  newGradLaplList_d = newGradLaplList;
  calc_grad_lapl (AinvList_d.data(), gradLaplList_d.data(),
                  newGradLaplList_d.data(), NumOrbitals,
                  RowStride, walkers.size());
  // Now copy data into the output matrices
  gradLapl_host = gradLapl_d;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    for(int iat=0; iat < NumPtcls; iat++)
    {
      GradType g(gradLapl_host[4*(iw*RowStride + iat)+0],
                 gradLapl_host[4*(iw*RowStride + iat)+1],
                 gradLapl_host[4*(iw*RowStride + iat)+2]);
      grads(iw,iat+FirstIndex) += g;
#ifdef QMC_COMPLEX
      //YingWai's trial fix; need to revise the dot product
      // should be converted into "#if defined QMC_COMPLEX" later
      lapl(iw,iat+FirstIndex)  += gradLapl_host[4*(iw*RowStride + iat)+3] - (CudaValueType)dot(g,g);
      ValueType lapl_temp = lapl(iw,iat+FirstIndex);
      if ( std::isnan(std::real(lapl_temp)) || std::isnan(std::imag(lapl_temp)) )
#else
      lapl(iw,iat+FirstIndex)  += gradLapl_host[4*(iw*RowStride + iat)+3] - dot(g,g);
      if (std::isnan(lapl(iw,iat+FirstIndex)))
#endif
      {
        fprintf (stderr, "Offending walker = %d\n", iw);
        char name[1000];
        gethostname(name, 1000);
        fprintf (stderr, "Offending hostname = %s\n", name);
        int dev;
        cudaGetDevice(&dev);
        fprintf (stderr, "Offending device = %d\n", dev);
        gpu::host_vector<CudaValueType> host_data;
        host_data = walkers[iw]->cuda_DataSet;
        FILE *Amat = fopen ("Amat.dat", "w");
        FILE *Ainv = fopen ("Ainv.dat", "w");
        FILE *Lmat = fopen ("Alapl.dat", "w");
        FILE *Gmat = fopen ("Agrad.dat", "w");
        for (int i=0; i<NumPtcls; i++)
        {
          for (int j=0; j<NumPtcls; j++)
          {
#ifdef QMC_COMPLEX
            fprintf (Amat, "%14.8e+%14.8ei ", host_data[AOffset+i*RowStride+j].real(), host_data[AOffset+i*RowStride+j].imag());
            fprintf (Ainv, "%14.8e+%14.8ei ", host_data[AinvOffset+i*RowStride+j].real(), host_data[AinvOffset+i*RowStride+j].imag());
            fprintf (Lmat, "%14.8e+%14.8ei ", host_data[gradLaplOffset+(4*i+3)*RowStride+j].real(), host_data[gradLaplOffset+(4*i+3)*RowStride+j].imag());
            for (int k=0; k<3; k++)
              fprintf (Gmat, "%14.8e+%14.8ei ", host_data[gradLaplOffset+(4*i+k)*RowStride+j].real(), host_data[gradLaplOffset+(4*i+k)*RowStride+j].imag());
#else
            fprintf (Amat, "%14.8e ", host_data[AOffset+i*RowStride+j]);
            fprintf (Ainv, "%14.8e ", host_data[AinvOffset+i*RowStride+j]);
            fprintf (Lmat, "%14.8e ", host_data[gradLaplOffset+(4*i+3)*RowStride+j]);
            for (int k=0; k<3; k++)
              fprintf (Gmat, "%14.8e ", host_data[gradLaplOffset+(4*i+k)*RowStride+j]);
#endif
          }
          fprintf (Amat, "\n");
          fprintf (Ainv, "\n");
          fprintf (Lmat, "\n");
          fprintf (Gmat, "\n");
        }
        fclose (Amat);
        fclose (Ainv);
        fclose (Lmat);
        fclose (Gmat);
        abort();
// 	  if (std::isnan(host_data[AinvOffset+i*RowStride+j]))
// 	    std::cerr << "NAN in inverse at (" << i << "," << j << ")\n";
        std::cerr << "NAN in walker " << iw << ", iat " << iat + FirstIndex
             << "  grad = " << grads(iw,iat+FirstIndex)
             << "  lapl = " << gradLapl_host[4*(iw*RowStride + iat)+3] << std::endl;
        fprintf (stderr, "grad-lapl row:\n");
        fprintf (stderr, "r = %1.8f %1.8f %1.8f\n",
                 walkers[iw]->R[iat+FirstIndex][0],
                 walkers[iw]->R[iat+FirstIndex][1],
                 walkers[iw]->R[iat+FirstIndex][2]);
        for (int orb=0; orb<NumPtcls; orb++)
#ifdef QMC_COMPLEX
          fprintf (stderr, "%1.10e+%1.10ei %1.10e+%1.10ei %1.10e+%1.10ei %1.10e+%1.10ei \n",
                   host_data[gradLaplOffset +(4*iat+0)*RowStride+orb].real(),
                   host_data[gradLaplOffset +(4*iat+0)*RowStride+orb].imag(),
                   host_data[gradLaplOffset +(4*iat+1)*RowStride+orb].real(),
                   host_data[gradLaplOffset +(4*iat+1)*RowStride+orb].imag(),
                   host_data[gradLaplOffset +(4*iat+2)*RowStride+orb].real(),
                   host_data[gradLaplOffset +(4*iat+2)*RowStride+orb].imag(),
                   host_data[gradLaplOffset +(4*iat+3)*RowStride+orb].real(),
                   host_data[gradLaplOffset +(4*iat+3)*RowStride+orb].imag());
#else
          fprintf (stderr, "%1.10e %1.10e %1.10e %1.10e \n",
                   host_data[gradLaplOffset +(4*iat+0)*RowStride+orb],
                   host_data[gradLaplOffset +(4*iat+1)*RowStride+orb],
                   host_data[gradLaplOffset +(4*iat+2)*RowStride+orb],
                   host_data[gradLaplOffset +(4*iat+3)*RowStride+orb]);
#endif
      }
    }
  }
#ifdef CUDA_DEBUG
  // Now do it on the CPU
  gpu::host_vector<CudaValueType> host_data;
  GradMatrix_t cpu_grads(grads.rows(), grads.cols());
  ValueMatrix_t cpu_lapl(grads.rows(), grads.cols());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    host_data = walkers[iw]->cuda_DataSet;
    for (int iat=0; iat < NumPtcls; iat++)
    {
      cpu_grads(iw,iat+FirstIndex) = GradType();
      cpu_lapl (iw,iat+FirstIndex) = ValueType();
      for (int iorb=0; iorb<NumOrbitals; iorb++)
      {
        cpu_lapl(iw,iat+FirstIndex) += host_data[AinvOffset+NumPtcls*iorb+iat] *
                                       host_data[gradLaplOffset+(4*iat+3)*NumOrbitals + iorb];
      }
      fprintf (stderr, "CPU lapl = %10.6e   GPU lapl = %10.6e\n",
               cpu_lapl(iw,iat+FirstIndex), lapl(iw,iat+FirstIndex));
    }
  }
#endif
}

void
DiracDeterminantEval<Batching::BATCHED>::NLratios (MCWalkerConfiguration &W,
                                std::vector<NLjob> &jobList,
                                std::vector<PosType> &quadPoints,
                                std::vector<ValueType> &psi_ratios)
{
  // DEBUG!
  // std::vector<ValueType> cpu_ratios;
  // cpu_ratios = psi_ratios;
  // NLratios_CPU (W, jobList, quadPoints, cpu_ratios);
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int posIndex=0, numJobs=0;
  std::vector<PosType> posBuffer[2];
  int rowIndex = 0;
  std::vector<ValueType*> ratio_pointers[2];
  bool hasResults = false;
  for(int i = 0; i < 2; ++i)
  {
    posBuffer[i].clear();
    ratio_pointers[i].clear();
    NLAinvList_host[i].clear();
    NLnumRatioList_host[i].clear();
    NLelecList_host[i].clear();
    NLratioList_host[i].clear();
    RatioRowList_host[i].clear();
  }
  int ijob = 0;
  int counter = 0;
  //First batch of data (PT)
  for (; ijob < jobList.size(); ++ijob)
  {
    NLjob &job = jobList[ijob];
    int numQuad = job.numQuadPoints;
    int elec    = job.elec;
    // Check if this electron belongs to this determinant
    if (elec < FirstIndex || elec >= LastIndex)
    {
      posIndex += numQuad;
      continue;
    }
    // Check to see if the buffer is full
    if (rowIndex + numQuad > NLrowBufferRows)
    {
      break;
    }
    int iw = job.walker;
    NLAinvList_host[counter].push_back(&(walkers[iw]->cuda_DataSet.data()[AinvOffset]));
    NLnumRatioList_host[counter].push_back(numQuad);
    NLelecList_host[counter].push_back(job.elec-FirstIndex);
    NLratioList_host[counter].push_back(&(NLratios_d[counter].data()[rowIndex]));
    RatioRowList_host[counter].push_back(&(NLrowBuffer_d.data()[rowIndex*RowStride]));
    for (int iq=0; iq < numQuad; iq++)
    {
      posBuffer[counter].push_back(quadPoints[posIndex]);
      ratio_pointers[counter].push_back(&(psi_ratios[posIndex]));
      posIndex++;
    }
    rowIndex += numQuad;
    numJobs++;
  }
  //Next batches (PT)
  while(ijob < jobList.size())
  {
    NLjob &job = jobList[ijob];
    int numQuad = job.numQuadPoints;
    int elec    = job.elec;
    if(rowIndex + numQuad > NLrowBufferRows)
    {
      // Compute orbital rows
      Phi->evaluate (posBuffer[counter], SplineRowList_d);
      // Compute ratios
      NLAinvList_d.asyncCopy(NLAinvList_host[counter]);
      NLnumRatioList_d.asyncCopy(NLnumRatioList_host[counter]);
      NLelecList_d.asyncCopy(NLelecList_host[counter]);
      NLratioList_d.asyncCopy(NLratioList_host[counter]);
      RatioRowList_d.asyncCopy(RatioRowList_host[counter]);
      calc_many_ratios (NLAinvList_d.data(), RatioRowList_d.data(),
                        NLratioList_d.data(), NLnumRatioList_d.data(),
                        NumOrbitals, RowStride, NLelecList_d.data(),
                        numJobs);
      // Write ratios out output vector
      rowIndex=0;
      numJobs=0;
      counter = (counter+1)%2;
    }
    if(hasResults)
    {
      cudaStreamSynchronize(gpu::memoryStream);
      for (int i=0; i<ratio_pointers[counter].size(); i++)
        *(ratio_pointers[counter][i]) *= NLratios_host[i];
      hasResults = false;
      ratio_pointers[counter].clear();
    }
    for (; ijob < jobList.size(); ijob++)
    {
      NLjob &job = jobList[ijob];
      int numQuad = job.numQuadPoints;
      int elec    = job.elec;
      // Check if this electron belongs to this determinant
      if (elec < FirstIndex || elec >= LastIndex)
      {
        posIndex += numQuad;
        continue;
      }
      // Check to see if the buffer is full
      if (rowIndex + numQuad > NLrowBufferRows)
      {
        break;
      }
      int iw = job.walker;
      NLAinvList_host[counter].push_back(&(walkers[iw]->cuda_DataSet.data()[AinvOffset]));
      NLnumRatioList_host[counter].push_back(numQuad);
      NLelecList_host[counter].push_back(job.elec-FirstIndex);
      NLratioList_host[counter].push_back(&(NLratios_d[counter].data()[rowIndex]));
      RatioRowList_host[counter].push_back(&(NLrowBuffer_d.data()[rowIndex*RowStride]));
      for (int iq=0; iq < numQuad; iq++)
      {
        posBuffer[counter].push_back(quadPoints[posIndex]);
        ratio_pointers[counter].push_back(&(psi_ratios[posIndex]));
        posIndex++;
      }
      rowIndex += numQuad;
      numJobs++;
    }
    NLratios_host.asyncCopy(NLratios_d[(counter+1)%2]);
    hasResults = true;
    // Reset counters
    posBuffer[(counter+1)%2].clear();
    //ratio_pointers[(counter+1)%2].clear();
    NLAinvList_host[(counter+1)%2].clear();
    NLnumRatioList_host[(counter+1)%2].clear();
    NLelecList_host[(counter+1)%2].clear();
    NLratioList_host[(counter+1)%2].clear();
    RatioRowList_host[(counter+1)%2].clear();
  }
  if(hasResults)
  {
    cudaStreamSynchronize(gpu::memoryStream);
    for (int i=0; i<ratio_pointers[(counter+1)%2].size(); i++)
      *(ratio_pointers[(counter+1)%2][i]) *= NLratios_host[i];
  }
  if (posBuffer[counter].size())
  {
    // Compute whatever remains in the buffer
    // Compute orbital rows
    Phi->evaluate (posBuffer[counter], SplineRowList_d);
    // Compute ratios
    NLAinvList_d.asyncCopy(NLAinvList_host[counter]);
    NLnumRatioList_d.asyncCopy(NLnumRatioList_host[counter]);
    NLelecList_d.asyncCopy(NLelecList_host[counter]);
    NLratioList_d.asyncCopy(NLratioList_host[counter]);
    RatioRowList_d.asyncCopy(RatioRowList_host[counter]);
    calc_many_ratios (NLAinvList_d.data(), RatioRowList_d.data(),
                      NLratioList_d.data(), NLnumRatioList_d.data(),
                      NumOrbitals, RowStride, NLelecList_d.data(),
                      numJobs);
    // Write ratios out output vector
    NLratios_host = NLratios_d[counter];
    for (int i=0; i<ratio_pointers[counter].size(); i++)
      *(ratio_pointers[counter][i]) *= NLratios_host[i];
  }
  // DEBUG DEBUG DEBUG
  // for (int i=0; i<psi_ratios.size(); i++) {
  //   double diff = psi_ratios[i] - cpu_ratios[i];
  //   if (std::abs(diff) > 1.0e-8)
  // 	fprintf (stderr, "i=%d  GPU=%1.12f  CPU=%1.12f  FirstIndex=%d\n",
  // 		 i, psi_ratios[i], cpu_ratios[i], FirstIndex);
  // }
}



}

