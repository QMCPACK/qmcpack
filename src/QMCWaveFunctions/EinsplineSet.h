//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EINSPLINE_SET_H
#define QMCPLUSPLUS_EINSPLINE_SET_H

#include "Configuration.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "Utilities/TimerManager.h"
#include "spline/einspline_engine.hpp"
#ifdef QMC_CUDA
#include "einspline/multi_bspline_create_cuda.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/AtomicOrbitalCuda.h"
#endif

namespace qmcplusplus
{
class EinsplineSetBuilder;

class EinsplineSet : public SPOSet
{
  friend class EinsplineSetBuilder;

public:
  //////////////////////
  // Type definitions //
  //////////////////////
  using UnitCellType = CrystalLattice<ParticleSet::Scalar_t, OHMMS_DIM>;

  ///////////
  // Flags //
  ///////////
  /// True if all Lattice is diagonal, i.e. 90 degree angles
  bool Orthorhombic;
  /// True if we are using localize orbitals
  bool Localized;
  /// True if we are tiling the primitive cell
  bool Tiling;

  //////////////////////////
  // Lattice and geometry //
  //////////////////////////
  TinyVector<int, 3> TileFactor;
  Tensor<int, OHMMS_DIM> TileMatrix;
  UnitCellType SuperLattice, PrimLattice;
  /// The "Twist" variables are in reduced coords, i.e. from 0 to1.
  /// The "k" variables are in Cartesian coordinates.
  PosType TwistVector, kVector;
  /// This stores which "true" twist vector this clone is using.
  /// "True" indicates the physical twist angle after untiling
  int TwistNum;
  /// metric tensor to handle generic unitcell
  Tensor<RealType, OHMMS_DIM> GGt;

  int NumValenceOrbs;

public:
  UnitCellType GetLattice();
  void resetParameters(const opt_variables_type& active) override {}
  void resetSourceParticleSet(ParticleSet& ions);
  void setOrbitalSetSize(int norbs) override;
  inline std::string Type() { return "EinsplineSet"; }
  EinsplineSet() : TwistNum(0), NumValenceOrbs(0) { className = "EinsplineSet"; }
};

////////////////////////////////////////////////////////////////////
// This is just a template trick to avoid template specialization //
// in EinsplineSetExtended.                                       //
////////////////////////////////////////////////////////////////////
template<typename StorageType, int dim>
struct MultiOrbitalTraits
{};

template<>
struct MultiOrbitalTraits<double, 2>
{
  using SplineType = multi_UBspline_2d_d;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_2d_d_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<std::complex<double>, 2>
{
  using SplineType = multi_UBspline_2d_z;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_2d_z_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<float, 2>
{
  using SplineType = multi_UBspline_2d_s;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_2d_s_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<std::complex<float>, 2>
{
  using SplineType = multi_UBspline_2d_c;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_2d_c_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<double, 3>
{
  using SplineType = multi_UBspline_3d_d;
  using BCType     = BCtype_d;
  using DataType   = double;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_3d_d_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<std::complex<double>, 3>
{
  using SplineType = multi_UBspline_3d_z;
  using BCType     = BCtype_z;
  using DataType   = std::complex<double>;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_3d_z_cuda;
#endif
};


template<>
struct MultiOrbitalTraits<float, 3>
{
  using SplineType = multi_UBspline_3d_s;
  using BCType     = BCtype_s;
  using DataType   = float;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_3d_s_cuda;
#endif
};

template<>
struct MultiOrbitalTraits<std::complex<float>, 3>
{
  using SplineType = multi_UBspline_3d_c;
  using BCType     = BCtype_c;
  using DataType   = std::complex<float>;
#ifdef QMC_CUDA
  using CudaSplineType = multi_UBspline_3d_c_cuda;
#endif
};


#ifdef QMC_CUDA
template<typename StoreType, typename CudaPrec>
struct StorageTypeConverter;
template<>
struct StorageTypeConverter<double, double>
{
  using CudaStorageType = double;
};
template<>
struct StorageTypeConverter<double, float>
{
  using CudaStorageType = float;
};
template<>
struct StorageTypeConverter<std::complex<double>, float>
{
  using CudaStorageType = std::complex<float>;
};
template<>
struct StorageTypeConverter<std::complex<double>, std::complex<double>>
{
  using CudaStorageType = std::complex<double>;
};
template<>
struct StorageTypeConverter<std::complex<double>, double>
{
  using CudaStorageType = std::complex<double>;
};
#endif


////////////////////////////////////////////////////////////////////
// Template class for evaluating multiple extended Bloch orbitals //
// quickly.  Currently uses einspline library.                    //
////////////////////////////////////////////////////////////////////
template<typename StorageType>
class EinsplineSetExtended : public EinsplineSet
{
  friend class EinsplineSetBuilder;

protected:
  //////////////////////
  // Type definitions //
  //////////////////////
  //using UnitCellType = CrystalLattice<RealType,OHMMS_DIM>;
  using SplineType = typename MultiOrbitalTraits<StorageType, OHMMS_DIM>::SplineType;
  using BCType     = typename MultiOrbitalTraits<StorageType, OHMMS_DIM>::BCType;

  using StorageValueVector    = typename OrbitalSetTraits<StorageType>::ValueVector;
  using StorageGradVector     = typename OrbitalSetTraits<StorageType>::GradVector;
  using StorageHessVector     = typename OrbitalSetTraits<StorageType>::HessVector;
  using StorageGradHessVector = typename OrbitalSetTraits<StorageType>::GradHessVector;
  using RealValueVector       = Vector<double>;
  using ComplexValueVector    = Vector<std::complex<double>>;
  using RealGradVector        = Vector<TinyVector<double, OHMMS_DIM>>;
  using ComplexGradVector     = Vector<TinyVector<std::complex<double>, OHMMS_DIM>>;
  using RealHessType          = Tensor<double, OHMMS_DIM>;
  using ComplexHessType       = Tensor<std::complex<double>, OHMMS_DIM>;
  using RealHessVector        = Vector<RealHessType>;
  using RealHessMatrix        = Matrix<RealHessType>;
  using ComplexHessVector     = Vector<ComplexHessType>;
  using ComplexHessMatrix     = Matrix<ComplexHessType>;
  using RealValueMatrix       = Matrix<double>;
  using ComplexValueMatrix    = Matrix<std::complex<double>>;
  using RealGradMatrix        = Matrix<TinyVector<double, OHMMS_DIM>>;
  using ComplexGradMatrix     = Matrix<TinyVector<std::complex<double>, OHMMS_DIM>>;
  using RealGGGType           = TinyVector<RealHessType, 3>;
  using RealGGGVector         = Vector<RealGGGType>;
  using RealGGGMatrix         = Matrix<RealGGGType>;
  using ComplexGGGType        = TinyVector<ComplexHessType, 3>;
  using ComplexGGGVector      = Vector<ComplexGGGType>;
  using ComplexGGGMatrix      = Matrix<ComplexGGGType>;

  /////////////////////////////
  /// Orbital storage object //
  /////////////////////////////
  SplineType* MultiSpline;

  //////////////////////////////////////
  // Radial/Ylm orbitals around atoms //
  //////////////////////////////////////
  std::vector<AtomicOrbital<StorageType>> AtomicOrbitals;

  // First-order derivative w.r.t. the ion positions
  std::vector<TinyVector<SplineType*, OHMMS_DIM>> FirstOrderSplines;
  // Temporary storage for Eispline calls
  StorageValueVector storage_value_vector_, storage_lapl_vector_;
  StorageGradVector storage_grad_vector_;
  StorageHessVector storage_hess_vector_;
  StorageGradHessVector storage_grad_hess_vector_;

  // True if we should unpack this orbital into two copies
  std::vector<bool> MakeTwoCopies;
  // k-points for each orbital
  Vector<TinyVector<double, OHMMS_DIM>> kPoints;

  ///////////////////
  // Phase factors //
  ///////////////////
  Vector<double> phase;
  Vector<std::complex<double>> eikr;
  void computePhaseFactors(const TinyVector<double, OHMMS_DIM>& r);
  // For running at half G-vectors with real orbitals;
  // 0 if the twist is zero, 1 if the twist is G/2.
  TinyVector<int, OHMMS_DIM> HalfG;

  ////////////
  // Timers //
  ////////////
  NewTimer &ValueTimer, &VGLTimer, &VGLMatTimer;
  NewTimer& EinsplineTimer;

#ifdef QMC_CUDA
  // Cuda equivalents of the above
  using CudaStorageType = typename StorageTypeConverter<StorageType, CUDA_PRECISION>::CudaStorageType;
  using CudaSplineType  = typename MultiOrbitalTraits<CudaStorageType, OHMMS_DIM>::CudaSplineType;

  CudaSplineType* CudaMultiSpline;
  gpu::device_vector<CudaStorageType> CudaValueVector, CudaGradLaplVector;
  gpu::device_vector<CudaStorageType*> CudaValuePointers, CudaGradLaplPointers;
  std::vector<cudaIpcMemHandle_t> spline_rank_handles;
  std::vector<CudaStorageType*> spline_rank_pointers;
  std::vector<cudaEvent_t> spline_events;
  std::vector<cudaStream_t> spline_streams;
  int abort_counter  = 0;
  bool split_splines = false;
  void resize_cuda(int numWalkers);
  void get_split_spline_pointers();
  // Cuda equivalent
  gpu::device_vector<int> CudaMakeTwoCopies;
  gpu::device_vector<int> CudaTwoCopiesIndex;
  // Cuda equivalent
  gpu::device_vector<TinyVector<CUDA_PRECISION, OHMMS_DIM>> CudakPoints, CudakPoints_reduced;
  void applyPhaseFactors(gpu::device_vector<CudaStorageType*>& storageVector, gpu::device_vector<CTS::RealType*>& phi);
  // Data for vectorized evaluations
  gpu::host_vector<CTS::PosType> hostPos, hostPhasePos, NLhostPos;
  gpu::device_vector<CTS::PosType> cudapos, cudaphasepos, NLcudapos;
  gpu::host_vector<CTS::RealType> hostSign, NLhostSign;
  gpu::device_vector<CTS::RealType> cudaSign, NLcudaSign;
  // This stores the inverse of the lattice vector matrix in
  // GPU memory.
  gpu::device_vector<CTS::RealType> Linv_cuda, L_cuda;
  gpu::host_vector<CTS::RealType> L_host, Linv_host;
#endif

public:
  /** create MultiSpline
   * @param xyz_g grid data
   * @param xyz_bc boundary conditions
   */
  template<typename GT, typename BCT>
  void allocate(GT& xyz_g, BCT& xyz_bc, int nv)
  {
    SplineType* dummy = nullptr;
    MultiSpline       = einspline::create(dummy, xyz_g, xyz_bc, nv);
  }

  inline void resizeStorage(int n, int nvals)
  {
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
    storage_value_vector_.resize(n);
    storage_lapl_vector_.resize(n);
    storage_grad_vector_.resize(n);
    storage_hess_vector_.resize(n);
    storage_grad_hess_vector_.resize(n);
    phase.resize(n);
    eikr.resize(n);
    NumValenceOrbs = nvals;
  }

#if !defined(QMC_COMPLEX)
  // Real return values
  void evaluateValue(const ParticleSet& P, int iat, RealValueVector& psi) override;
  void evaluateVGL(const ParticleSet& P,
                   int iat,
                   RealValueVector& psi,
                   RealGradVector& dpsi,
                   RealValueVector& d2psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            RealValueMatrix& psi,
                            RealGradMatrix& dpsi,
                            RealValueMatrix& d2psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            RealValueMatrix& psi,
                            RealGradMatrix& dpsi,
                            RealHessMatrix& grad_grad_psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            RealValueMatrix& psi,
                            RealGradMatrix& dpsi,
                            RealHessMatrix& grad_grad_psi,
                            RealGGGMatrix& grad_grad_grad_logdet) override;

  //    void evaluate (const ParticleSet& P, const PosType& r, std::vector<double> &psi);
  // This is the gradient of the orbitals w.r.t. the ion iat
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          RealGradMatrix& gradphi) override;
  // Evaluate the gradient w.r.t. to ion iat of the gradient and
  // laplacian of the orbitals w.r.t. the electrons
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          RealGradMatrix& dphi,
                          RealHessMatrix& dgrad_phi,
                          RealGradMatrix& dlaplphi) override;
#else
  // Complex return values
  void evaluateValue(const ParticleSet& P, int iat, ComplexValueVector& psi) override;
  void evaluateVGL(const ParticleSet& P,
                   int iat,
                   ComplexValueVector& psi,
                   ComplexGradVector& dpsi,
                   ComplexValueVector& d2psi) override;
  void evaluateVGH(const ParticleSet& P,
                   int iat,
                   ComplexValueVector& psi,
                   ComplexGradVector& dpsi,
                   ComplexHessVector& grad_grad_psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ComplexValueMatrix& psi,
                            ComplexGradMatrix& dpsi,
                            ComplexValueMatrix& d2psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ComplexValueMatrix& psi,
                            ComplexGradMatrix& dpsi,
                            ComplexHessMatrix& grad_grad_psi) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ComplexValueMatrix& psi,
                            ComplexGradMatrix& dpsi,
                            ComplexHessMatrix& grad_grad_psi,
                            ComplexGGGMatrix& grad_grad_grad_logdet) override;
#endif

#ifdef QMC_CUDA
  void finalizeConstruction() override;

  // Vectorized evaluation functions
#if !defined(QMC_COMPLEX)
  void evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::RealType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::RealType*>& phi) override;
  inline void evaluate(std::vector<Walker_t*>& walkers,
                       std::vector<PosType>& newpos,
                       gpu::device_vector<CTS::RealType*>& phi,
                       gpu::device_vector<CTS::RealType*>& grad_lapl,
                       int row_stride) override
  {
    evaluate(walkers, newpos, phi, grad_lapl, row_stride, 0, false);
  }
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::RealType*>& phi,
                gpu::device_vector<CTS::RealType*>& grad_lapl,
                int row_stride,
                int k,
                bool klinear) override;

  void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi) override;
#else
  void evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::ComplexType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::ComplexType*>& phi) override;
  inline void evaluate(std::vector<Walker_t*>& walkers,
                       std::vector<PosType>& newpos,
                       gpu::device_vector<CTS::ComplexType*>& phi,
                       gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                       int row_stride) override
  {
    evaluate(walkers, newpos, phi, grad_lapl, row_stride, 0, false);
  }
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::ComplexType*>& phi,
                gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                int row_stride,
                int k,
                bool klinear) override;
  void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi) override;
#endif
#endif

  void resetParameters(const opt_variables_type& active) override;
  void setOrbitalSetSize(int norbs) override;
  std::string Type();

  void registerTimers();
  PosType get_k(int orb) override { return kPoints[orb]; }


  std::unique_ptr<SPOSet> makeClone() const override;

  EinsplineSetExtended()
      : MultiSpline(NULL),
        ValueTimer(*timer_manager.createTimer("EinsplineSetExtended::ValueOnly")),
        VGLTimer(*timer_manager.createTimer("EinsplineSetExtended::VGL")),
        VGLMatTimer(*timer_manager.createTimer("EinsplineSetExtended::VGLMatrix")),
        EinsplineTimer(*timer_manager.createTimer("libeinspline"))
#ifdef QMC_CUDA
        ,
        CudaMultiSpline(NULL),
        CudaValueVector("EinsplineSetExtended::CudaValueVector"),
        CudaGradLaplVector("EinsplineSetExtended::CudaGradLaplVector"),
        CudaValuePointers("EinsplineSetExtended::CudaValuePointers"),
        CudaGradLaplPointers("EinsplineSetExtended::CudaGradLaplPointers"),
        CudaMakeTwoCopies("EinsplineSetExtended::CudaMakeTwoCopies"),
        CudaTwoCopiesIndex("EinsplineSetExtended::CudaTwoCopiesIndex"),
        CudakPoints("EinsplineSetExtended::CudakPoints"),
        CudakPoints_reduced("EinsplineSetExtended::CudakPoints_reduced"),
        cudapos("EinsplineSetExtended::cudapos"),
        NLcudapos("EinsplineSetExtended::NLcudapos"),
        cudaSign("EinsplineSetExtended::cudaSign"),
        NLcudaSign("EinsplineSetExtended::NLcudaSign"),
        Linv_cuda("EinsplineSetExtended::Linv_cuda"),
        L_cuda("EinsplineSetExtended::L_cuda")
#endif
  {
    className = "EinsplineSetExtended";
    for (int i = 0; i < OHMMS_DIM; i++)
      HalfG[i] = 0;
  }
};

#ifdef QMC_CUDA
template<typename T>
struct AtomicSplineJob
{
  T dist, SplineDelta;
  T rhat[OHMMS_DIM];
  int lMax, YlmIndex;
  T* SplineCoefs;
  T *phi, *grad_lapl;
  T PAD[3];
  //T PAD[(64 - (2*OHMMS_DIM*sizeof(T) + 2*sizeof(int) + 2*sizeof(T*)))/sizeof(T)];
};

template<typename T>
struct AtomicPolyJob
{
  T dist, SplineDelta;
  T rhat[OHMMS_DIM];
  int lMax, PolyOrder, YlmIndex;
  T* PolyCoefs;
  T *phi, *grad_lapl;
  T PAD[2];
  //T PAD[(64 - (2*OHMMS_DIM*sizeof(T) + 2*sizeof(int) + 2*sizeof(T*)))/sizeof(T)];
};


template<typename StorageType>
class EinsplineSetHybrid : public EinsplineSetExtended<StorageType>
{
  friend class EinsplineSetBuilder;

protected:
  int a;
  //////////////////////
  // Type definitions //
  //////////////////////
  using CTS             = CUDAGlobalTypes;
  using Walker_t        = typename EinsplineSetExtended<StorageType>::Walker_t;
  using PosType         = typename EinsplineSetExtended<StorageType>::PosType;
  using CudaStorageType = typename EinsplineSetExtended<StorageType>::CudaStorageType;

  std::vector<gpu::device_vector<CTS::RealType>> AtomicSplineCoefs_GPU, AtomicPolyCoefs_GPU;
  gpu::device_vector<AtomicOrbitalCuda<CTS::RealType>> AtomicOrbitals_GPU;

  // gpu::host_vector<AtomicPolyJob<CTS::RealType> >   AtomicPolyJobs_CPU;
  // gpu::device_vector<AtomicPolyJob<CTS::RealType> >   AtomicPolyJobs_GPU;
  // gpu::host_vector<AtomicSplineJob<CTS::RealType> > AtomicSplineJobs_CPU;
  // gpu::device_vector<AtomicSplineJob<CTS::RealType> > AtomicSplineJobs_GPU;

  gpu::device_vector<HybridJobType> HybridJobs_GPU;
  gpu::device_vector<CTS::RealType> IonPos_GPU;
  gpu::device_vector<CTS::RealType> CutoffRadii_GPU, PolyRadii_GPU;
  gpu::device_vector<HybridData<CTS::RealType>> HybridData_GPU;

  gpu::device_vector<CTS::RealType> Ylm_GPU;
  gpu::device_vector<CTS::RealType*> Ylm_ptr_GPU, dYlm_dtheta_ptr_GPU, dYlm_dphi_ptr_GPU;
  gpu::host_vector<CTS::RealType*> Ylm_ptr_CPU, dYlm_dtheta_ptr_CPU, dYlm_dphi_ptr_CPU;
  gpu::device_vector<CTS::RealType> rhats_GPU;
  gpu::host_vector<CTS::RealType> rhats_CPU;
  gpu::device_vector<int> JobType;

  // Vectors for 3D Bspline evaluation
  gpu::device_vector<CTS::RealType> BsplinePos_GPU;
  gpu::host_vector<CTS::RealType> BsplinePos_CPU;
  gpu::device_vector<CudaStorageType*> BsplineVals_GPU, BsplineGradLapl_GPU;
  gpu::host_vector<CudaStorageType*> BsplineVals_CPU, BsplineGradLapl_CPU;

  // The maximum lMax across all atomic orbitals
  int lMax;
  int numlm, NumOrbitals, Ylm_BS;
  // Stores the maximum number of walkers that can be handled by currently
  // allocated GPU memory.  Must resize if we have more walkers than this.
  int CurrentWalkers;

  //////////////////////////////
  /// Orbital storage objects //
  //////////////////////////////

  ////////////
  // Timers //
  ////////////
  // Data for vectorized evaluations

  void sort_electrons(std::vector<PosType>& pos);

public:
  void finalizeConstruction() override;
  //    void registerTimers();

  // Resize cuda objects
  void resize_cuda(int numwalkers);

  // Vectorized evaluation functions
#if !defined(QMC_COMPLEX)
  void evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::RealType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::RealType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::RealType*>& phi,
                gpu::device_vector<CTS::RealType*>& grad_lapl,
                int row_stride) override;
  void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi) override;
#else
  void evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::ComplexType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::ComplexType*>& phi) override;
  void evaluate(std::vector<Walker_t*>& walkers,
                std::vector<PosType>& newpos,
                gpu::device_vector<CTS::ComplexType*>& phi,
                gpu::device_vector<CTS::ComplexType*>& grad_lapl,
                int row_stride) override;
  void evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi) override;
#endif

  std::string Type();

  std::unique_ptr<SPOSet> makeClone() const override;

  EinsplineSetHybrid();
};

#endif

//  template<typename StorageType>
//  inline void EinsplineSetExtended<StorageType>::computePhaseFactors
//  (TinyVector<RealType,OHMMS_DIM> r)
//  {
//    for (int i=0; i<kPoints.size(); i++) phase[i] = -dot(r, kPoints[i]);
//    eval_e2iphi(phase,eikr);
////#ifdef HAVE_MKL
////    for (int i=0; i<kPoints.size(); i++)
////      phase[i] = -dot(r, kPoints[i]);
////    vzCIS(OrbitalSetSize, phase, (double*)eikr.data());
////#else
////    double s, c;
////    for (int i=0; i<kPoints.size(); i++) {
////      phase[i] = -dot(r, kPoints[i]);
////      qmcplusplus::sincos (phase[i], &s, &c);
////      eikr[i] = std::complex<double>(c,s);
////    }
////#endif
//  }
//


} // namespace qmcplusplus
#endif
