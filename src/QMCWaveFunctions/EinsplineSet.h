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
  void resetSourceParticleSet(ParticleSet& ions);
  void setOrbitalSetSize(int norbs) override;
  inline std::string Type() { return "EinsplineSet"; }
  EinsplineSet(const std::string& my_name) : SPOSet(my_name), TwistNum(0), NumValenceOrbs(0) {}

  virtual std::string getClassName() const override { return "EinsplineSet"; }
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
};

template<>
struct MultiOrbitalTraits<std::complex<double>, 2>
{
  using SplineType = multi_UBspline_2d_z;
};

template<>
struct MultiOrbitalTraits<float, 2>
{
  using SplineType = multi_UBspline_2d_s;
};

template<>
struct MultiOrbitalTraits<std::complex<float>, 2>
{
  using SplineType = multi_UBspline_2d_c;
};

template<>
struct MultiOrbitalTraits<double, 3>
{
  using SplineType = multi_UBspline_3d_d;
  using BCType     = BCtype_d;
  using DataType   = double;
};

template<>
struct MultiOrbitalTraits<std::complex<double>, 3>
{
  using SplineType = multi_UBspline_3d_z;
  using BCType     = BCtype_z;
  using DataType   = std::complex<double>;
};


template<>
struct MultiOrbitalTraits<float, 3>
{
  using SplineType = multi_UBspline_3d_s;
  using BCType     = BCtype_s;
  using DataType   = float;
};

template<>
struct MultiOrbitalTraits<std::complex<float>, 3>
{
  using SplineType = multi_UBspline_3d_c;
  using BCType     = BCtype_c;
  using DataType   = std::complex<float>;
};

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
  /** kpoints for each unique orbitals.
   * Note: for historic reason, this sign is opposite to what was used in DFT when orbitals were generated.
   * Changing the sign requires updating all the evaluation code.
   */
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

  void setOrbitalSetSize(int norbs) override;
  std::string Type();

  void registerTimers();
  PosType get_k(int orb) override { return kPoints[orb]; }

  virtual std::string getClassName() const override { return "EinsplineSetExtended"; }
  bool hasIonDerivs() const override { return true; }

  std::unique_ptr<SPOSet> makeClone() const override;

  EinsplineSetExtended(const std::string& my_name)
      : EinsplineSet(my_name),
        MultiSpline(NULL),
        ValueTimer(createGlobalTimer("EinsplineSetExtended::ValueOnly")),
        VGLTimer(createGlobalTimer("EinsplineSetExtended::VGL")),
        VGLMatTimer(createGlobalTimer("EinsplineSetExtended::VGLMatrix")),
        EinsplineTimer(createGlobalTimer("libeinspline"))
  {
    for (int i = 0; i < OHMMS_DIM; i++)
      HalfG[i] = 0;
  }
};

} // namespace qmcplusplus
#endif
