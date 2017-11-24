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
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "QMCWaveFunctions/MuffinTin.h"
#include "Utilities/NewTimer.h"
#include <spline/einspline_engine.hpp>
#ifdef QMC_CUDA
#include <einspline/multi_bspline_create_cuda.h>
#include "QMCWaveFunctions/AtomicOrbitalCuda.h"
#endif

namespace qmcplusplus
{

class EinsplineSetBuilder;

class EinsplineSet : public SPOSetBase
{
  friend class EinsplineSetBuilder;
public:
  //////////////////////
  // Type definitions //
  //////////////////////
  typedef CrystalLattice<ParticleSet::Scalar_t,OHMMS_DIM> UnitCellType;

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
  TinyVector<int,3> TileFactor;
  Tensor<int,OHMMS_DIM> TileMatrix;
  UnitCellType SuperLattice, PrimLattice;
  /// The "Twist" variables are in reduced coords, i.e. from 0 to1.
  /// The "k" variables are in Cartesian coordinates.
  PosType TwistVector, kVector;
  /// This stores which "true" twist vector this clone is using.
  /// "True" indicates the physical twist angle after untiling
  int TwistNum;
  /// metric tensor to handle generic unitcell
  Tensor<RealType,OHMMS_DIM> GGt;

  ///////////////////////////////////////////////
  // Muffin-tin orbitals from LAPW calculation //
  ///////////////////////////////////////////////
  std::vector<MuffinTinClass> MuffinTins;
  int NumValenceOrbs, NumCoreOrbs;

public:
  UnitCellType GetLattice();
  virtual void resetParameters(const opt_variables_type& active)
  {}
  void resetTargetParticleSet(ParticleSet& e);
  void resetSourceParticleSet(ParticleSet& ions);
  void setOrbitalSetSize(int norbs);
  inline std::string Type()
  {
    return "EinsplineSet";
  }
  EinsplineSet() :  TwistNum(0), NumValenceOrbs(0), NumCoreOrbs(0)
  {
    className = "EinsplineSet";
  }
};

////////////////////////////////////////////////////////////////////
// This is just a template trick to avoid template specialization //
// in EinsplineSetExtended.                                       //
////////////////////////////////////////////////////////////////////
template<typename StorageType, int dim>  struct MultiOrbitalTraits {};

template<> struct MultiOrbitalTraits<double,2>
{
  typedef multi_UBspline_2d_d SplineType;
#ifdef QMC_CUDA
  typedef multi_UBspline_2d_d_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<std::complex<double>,2>
{
  typedef multi_UBspline_2d_z SplineType;
#ifdef QMC_CUDA
  typedef multi_UBspline_2d_z_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<float,2>
{
  typedef multi_UBspline_2d_s SplineType;
#ifdef QMC_CUDA
  typedef multi_UBspline_2d_s_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<std::complex<float>,2>
{
  typedef multi_UBspline_2d_c SplineType;
#ifdef QMC_CUDA
  typedef multi_UBspline_2d_c_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<double,3>
{
  typedef multi_UBspline_3d_d SplineType;
  typedef BCtype_d BCType;
  typedef double DataType;
#ifdef QMC_CUDA
  typedef multi_UBspline_3d_d_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<std::complex<double>,3>
{
  typedef multi_UBspline_3d_z SplineType;
  typedef BCtype_z BCType;
  typedef std::complex<double> DataType;
#ifdef QMC_CUDA
  typedef multi_UBspline_3d_z_cuda CudaSplineType;
#endif
};


template<> struct MultiOrbitalTraits<float,3>
{
  typedef multi_UBspline_3d_s SplineType;
  typedef BCtype_s BCType;
  typedef float DataType;
#ifdef QMC_CUDA
  typedef multi_UBspline_3d_s_cuda CudaSplineType;
#endif
};

template<> struct MultiOrbitalTraits<std::complex<float>,3>
{
  typedef multi_UBspline_3d_c SplineType;
  typedef BCtype_c BCType;
  typedef std::complex<float> DataType;
#ifdef QMC_CUDA
  typedef multi_UBspline_3d_c_cuda CudaSplineType;
#endif
};


#ifdef QMC_CUDA
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
  //typedef CrystalLattice<RealType,OHMMS_DIM> UnitCellType;
  typedef typename MultiOrbitalTraits<StorageType,OHMMS_DIM>::SplineType SplineType;
  typedef typename MultiOrbitalTraits<StorageType,OHMMS_DIM>::BCType     BCType;

  typedef typename OrbitalSetTraits<StorageType>::ValueVector_t StorageValueVector_t;
  typedef typename OrbitalSetTraits<StorageType>::GradVector_t  StorageGradVector_t;
  typedef typename OrbitalSetTraits<StorageType>::HessVector_t  StorageHessVector_t;
  typedef typename OrbitalSetTraits<StorageType>::GradHessVector_t  StorageGradHessVector_t;
  typedef Vector<double>                                        RealValueVector_t;
  typedef Vector<std::complex<double> >                              ComplexValueVector_t;
  typedef Vector<TinyVector<double,OHMMS_DIM> >                 RealGradVector_t;
  typedef Vector<TinyVector<std::complex<double>,OHMMS_DIM> >        ComplexGradVector_t;
  typedef Tensor<double,OHMMS_DIM>                              RealHessType;
  typedef Tensor<std::complex<double>,OHMMS_DIM>                     ComplexHessType;
  typedef Vector<RealHessType>                                  RealHessVector_t;
  typedef Matrix<RealHessType>                                  RealHessMatrix_t;
  typedef Vector<ComplexHessType>                               ComplexHessVector_t;
  typedef Matrix<ComplexHessType>                               ComplexHessMatrix_t;
  typedef Matrix<double>                                        RealValueMatrix_t;
  typedef Matrix<std::complex<double> >                              ComplexValueMatrix_t;
  typedef Matrix<TinyVector<double,OHMMS_DIM> >                 RealGradMatrix_t;
  typedef Matrix<TinyVector<std::complex<double>,OHMMS_DIM> >        ComplexGradMatrix_t;
  typedef TinyVector<RealHessType, 3>                           RealGGGType;
  typedef Vector<RealGGGType>                                   RealGGGVector_t;
  typedef Matrix<RealGGGType>                                   RealGGGMatrix_t;
  typedef TinyVector<ComplexHessType, 3>                        ComplexGGGType;
  typedef Vector<ComplexGGGType>                                ComplexGGGVector_t;
  typedef Matrix<ComplexGGGType>                                ComplexGGGMatrix_t;

  /////////////////////////////
  /// Orbital storage object //
  /////////////////////////////
  SplineType *MultiSpline;

  //////////////////////////////////////
  // Radial/Ylm orbitals around atoms //
  //////////////////////////////////////
  std::vector<AtomicOrbital<StorageType> > AtomicOrbitals;

  // First-order derivative w.r.t. the ion positions
  std::vector<TinyVector<SplineType*,OHMMS_DIM> > FirstOrderSplines;
  // Temporary storage for Eispline calls
  StorageValueVector_t StorageValueVector, StorageLaplVector;
  StorageGradVector_t  StorageGradVector;
  StorageHessVector_t  StorageHessVector;
  StorageGradHessVector_t  StorageGradHessVector;
  // Temporary storage used when blending functions
  StorageValueVector_t BlendValueVector, BlendLaplVector;
  StorageGradVector_t BlendGradVector;
  StorageHessVector_t  BlendHessVector;

  // True if we should unpack this orbital into two copies
  std::vector<bool>         MakeTwoCopies;
  // k-points for each orbital
  Vector<TinyVector<double,OHMMS_DIM> > kPoints;

  ///////////////////
  // Phase factors //
  ///////////////////
  Vector<double> phase;
  Vector<std::complex<double> > eikr;
  void computePhaseFactors(const TinyVector<double,OHMMS_DIM>& r);
  // For running at half G-vectors with real orbitals;
  // 0 if the twist is zero, 1 if the twist is G/2.
  TinyVector<int,OHMMS_DIM> HalfG;

  ////////////
  // Timers //
  ////////////
  NewTimer ValueTimer, VGLTimer, VGLMatTimer;
  NewTimer EinsplineTimer;

#ifdef QMC_CUDA
  // Cuda equivalents of the above
  typedef typename StorageTypeConverter<StorageType,CUDA_PRECISION>::CudaStorageType CudaStorageType;
  typedef typename MultiOrbitalTraits<CudaStorageType,OHMMS_DIM>::CudaSplineType CudaSplineType;

  CudaSplineType *CudaMultiSpline;
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
#endif

public:
  /** create MultiSpline
   * @param xyz_g grid data
   * @param xyz_bc boundary conditions
   */
  template<typename GT, typename BCT>
  void allocate(GT& xyz_g, BCT& xyz_bc, int nv)
  {
    SplineType* dummy=0;
    MultiSpline=einspline::create(dummy,xyz_g,xyz_bc,nv);
  }

  inline void resizeStorage(int n, int nvals, int ncores=0)
  {
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
    StorageValueVector.resize(n);
    BlendValueVector.resize(n);
    StorageLaplVector.resize(n);
    BlendLaplVector.resize(n);
    StorageGradVector.resize(n);
    BlendGradVector.resize(n);
    StorageHessVector.resize(n);
    StorageGradHessVector.resize(n);
    phase.resize(n);
    eikr.resize(n);
    NumValenceOrbs = nvals;
    NumCoreOrbs    = ncores;
  }

  // Real return values
  void evaluate(const ParticleSet& P, int iat, RealValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat, RealValueVector_t& psi,
                RealGradVector_t& dpsi, RealValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat, RealValueVector_t& psi,
                RealGradVector_t& dpsi, RealHessVector_t& grad_grad_psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            RealValueMatrix_t& psi, RealGradMatrix_t& dpsi,
                            RealValueMatrix_t& d2psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            RealValueMatrix_t& psi, RealGradMatrix_t& dpsi,
                            RealHessMatrix_t& grad_grad_psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            RealValueMatrix_t& psi, RealGradMatrix_t& dpsi,
                            RealHessMatrix_t& grad_grad_psi,
                            RealGGGMatrix_t& grad_grad_grad_logdet);

  //    void evaluate (const ParticleSet& P, const PosType& r, std::vector<double> &psi);
#if !defined(QMC_COMPLEX)
  // This is the gradient of the orbitals w.r.t. the ion iat
  void evaluateGradSource (const ParticleSet &P, int first, int last,
                           const ParticleSet &source, int iat_src,
                           RealGradMatrix_t &gradphi);
  // Evaluate the gradient w.r.t. to ion iat of the gradient and
  // laplacian of the orbitals w.r.t. the electrons
  void evaluateGradSource (const ParticleSet &P, int first, int last,
                           const ParticleSet &source, int iat_src,
                           RealGradMatrix_t &dphi,
                           RealHessMatrix_t  &dgrad_phi,
                           RealGradMatrix_t &dlaplphi);
  void evaluateGradSource (const ParticleSet &P, int first, int last,
                           const ParticleSet &source, int iat_src,
                           ComplexGradMatrix_t &gradphi);
#endif
  // Complex return values
  void evaluate(const ParticleSet& P, int iat, ComplexValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat, ComplexValueVector_t& psi,
                ComplexGradVector_t& dpsi, ComplexValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat, ComplexValueVector_t& psi,
                ComplexGradVector_t& dpsi, ComplexHessVector_t& grad_grad_psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ComplexValueMatrix_t& psi, ComplexGradMatrix_t& dpsi,
                            ComplexValueMatrix_t& d2psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ComplexValueMatrix_t& psi, ComplexGradMatrix_t& dpsi,
                            ComplexHessMatrix_t& grad_grad_psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ComplexValueMatrix_t& psi, ComplexGradMatrix_t& dpsi,
                            ComplexHessMatrix_t& grad_grad_psi,
                            ComplexGGGMatrix_t& grad_grad_grad_logdet);
#ifdef QMC_CUDA
  void initGPU();

  // Vectorized evaluation functions
  void evaluate (std::vector<Walker_t*> &walkers, int iat,
                 gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, int iat,
                 gpu::device_vector<CudaComplexType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaComplexType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaRealType*> &phi,
                 gpu::device_vector<CudaRealType*> &grad_lapl,
                 int row_stride);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaComplexType*> &phi,
                 gpu::device_vector<CudaComplexType*> &grad_lapl,
                 int row_stride);
  void evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi);
#endif

  void resetParameters(const opt_variables_type& active);
  void resetTargetParticleSet(ParticleSet& e);
  void setOrbitalSetSize(int norbs);
  std::string Type();

  void registerTimers();
  PosType get_k(int orb)
  {
    return kPoints[orb];
  }


  SPOSetBase* makeClone() const;

  EinsplineSetExtended() :
    ValueTimer  ("EinsplineSetExtended::ValueOnly"),
    VGLTimer    ("EinsplineSetExtended::VGL"),
    VGLMatTimer ("EinsplineSetExtended::VGLMatrix"),
    EinsplineTimer("libeinspline"),
    MultiSpline(NULL)
#ifdef QMC_CUDA
    , CudaMultiSpline(NULL),
    cudapos("EinsplineSetExtended::cudapos"),
    NLcudapos("EinsplineSetExtended::NLcudapos"),
    cudaSign("EinsplineSetExtended::cudaSign"),
    NLcudaSign("EinsplineSetExtended::NLcudaSign"),
    Linv_cuda("EinsplineSetExtended::Linv_cuda"),
    L_cuda("EinsplineSetExtended::L_cuda"),
    CudaValueVector("EinsplineSetExtended::CudaValueVector"),
    CudaGradLaplVector("EinsplineSetExtended::CudaGradLaplVector"),
    CudaValuePointers("EinsplineSetExtended::CudaValuePointers"),
    CudaGradLaplPointers("EinsplineSetExtended::CudaGradLaplPointers"),
    CudaMakeTwoCopies("EinsplineSetExtended::CudaMakeTwoCopies"),
    CudaTwoCopiesIndex("EinsplineSetExtended::CudaTwoCopiesIndex"),
    CudakPoints("EinsplineSetExtended::CudakPoints"),
    CudakPoints_reduced("EinsplineSetExtended::CudakPoints_reduced")
#endif
  {
    className = "EinsplineSetExtended";
    TimerManager.addTimer (&ValueTimer);
    TimerManager.addTimer (&VGLTimer);
    TimerManager.addTimer (&VGLMatTimer);
    TimerManager.addTimer (&EinsplineTimer);
    for (int i=0; i<OHMMS_DIM; i++)
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
  typedef typename EinsplineSetExtended<StorageType>::Walker_t     Walker_t;
  typedef typename EinsplineSetExtended<StorageType>::PosType      PosType;
  typedef typename EinsplineSetExtended<StorageType>::CudaRealType CudaRealType;
  typedef typename EinsplineSetExtended<StorageType>::CudaComplexType CudaComplexType;
  typedef typename EinsplineSetExtended<StorageType>::CudaStorageType CudaStorageType;

  std::vector<gpu::device_vector<CudaRealType> > AtomicSplineCoefs_GPU,
         AtomicPolyCoefs_GPU;
  gpu::device_vector<AtomicOrbitalCuda<CudaRealType> > AtomicOrbitals_GPU;

  // gpu::host_vector<AtomicPolyJob<CudaRealType> >   AtomicPolyJobs_CPU;
  // gpu::device_vector<AtomicPolyJob<CudaRealType> >   AtomicPolyJobs_GPU;
  // gpu::host_vector<AtomicSplineJob<CudaRealType> > AtomicSplineJobs_CPU;
  // gpu::device_vector<AtomicSplineJob<CudaRealType> > AtomicSplineJobs_GPU;

  gpu::device_vector<HybridJobType> HybridJobs_GPU;
  gpu::device_vector<CudaRealType>  IonPos_GPU;
  gpu::device_vector<CudaRealType>  CutoffRadii_GPU, PolyRadii_GPU;
  gpu::device_vector<HybridData<CudaRealType> > HybridData_GPU;

  gpu::device_vector<CudaRealType> Ylm_GPU;
  gpu::device_vector<CudaRealType*> Ylm_ptr_GPU, dYlm_dtheta_ptr_GPU, dYlm_dphi_ptr_GPU;
  gpu::host_vector<CudaRealType*> Ylm_ptr_CPU, dYlm_dtheta_ptr_CPU, dYlm_dphi_ptr_CPU;
  gpu::device_vector<CudaRealType> rhats_GPU;
  gpu::host_vector<CudaRealType> rhats_CPU;
  gpu::device_vector<int> JobType;

  // Vectors for 3D Bspline evaluation
  gpu::device_vector<CudaRealType> BsplinePos_GPU;
  gpu::host_vector<CudaRealType> BsplinePos_CPU;
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

  void sort_electrons(std::vector<PosType> &pos);

public:
  void initGPU();
  //    void registerTimers();

  // Resize cuda objects
  void resize_cuda(int numwalkers);

  // Vectorized evaluation functions
  void evaluate (std::vector<Walker_t*> &walkers, int iat,
                 gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, int iat,
                 gpu::device_vector<CudaComplexType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaComplexType*> &phi);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaRealType*> &phi,
                 gpu::device_vector<CudaRealType*> &grad_lapl,
                 int row_stride);
  void evaluate (std::vector<Walker_t*> &walkers, std::vector<PosType> &newpos,
                 gpu::device_vector<CudaComplexType*> &phi,
                 gpu::device_vector<CudaComplexType*> &grad_lapl,
                 int row_stride);
  void evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi);
  void evaluate (std::vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi);

  std::string Type();

  SPOSetBase* makeClone() const;

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
////      sincos (phase[i], &s, &c);
////      eikr[i] = std::complex<double>(c,s);
////    }
////#endif
//  }
//


}
#endif
