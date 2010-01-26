//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)             //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EINSPLINE_SET_H
#define QMCPLUSPLUS_EINSPLINE_SET_H

//#include <einspline/bspline.h>
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Optimize/VarList.h"
#include "QMCWaveFunctions/EinsplineOrb.h"
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "QMCWaveFunctions/MuffinTin.h"
#include "Utilities/NewTimer.h"
#include <einspline/multi_bspline_structs.h>
#include "Configuration.h"
#include "Numerics/e2iphi.h"

namespace qmcplusplus {

  class EinsplineSetBuilder;

  class EinsplineSet : public SPOSetBase
  {
    friend class EinsplineSetBuilder;
  protected:
    //////////////////////
    // Type definitions //
    //////////////////////
    typedef CrystalLattice<RealType,OHMMS_DIM> UnitCellType;
    
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
    vector<MuffinTinClass> MuffinTins;
    int NumValenceOrbs, NumCoreOrbs;
        
  public:  
    UnitCellType GetLattice();

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, 
		  ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& psi, GradMatrix_t& dpsi, 
		  ValueMatrix_t& d2psi);

    void resetTargetParticleSet(ParticleSet& e);
    void resetSourceParticleSet(ParticleSet& ions);
    void setOrbitalSetSize(int norbs);
    string Type();
    EinsplineSet() :  
      TwistNum(0), NumValenceOrbs(0), NumCoreOrbs(0)
    {
      className = "EinsplineSet";
    }
  };

  class EinsplineSetLocal : public EinsplineSet
  {
    friend class EinsplineSetBuilder;
  protected:
    /////////////////////
    // Orbital storage //
    /////////////////////
    /// Store the orbital objects.  Using template class allows us to
    /// avoid making separate real and complex versions of this class.
    //std::vector<EinsplineOrb<ValueType,OHMMS_DIM>*> Orbitals;
    std::vector<EinsplineOrb<complex<double>,OHMMS_DIM>*> Orbitals;

  public:
    SPOSetBase* makeClone() const;
    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, 
		  ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& psi, GradMatrix_t& dpsi, 
		  ValueMatrix_t& d2psi);

    void resetParameters(const opt_variables_type& active);

    EinsplineSetLocal() 
    {
      className = "EinsplineSetLocal";
    }
  };



  ////////////////////////////////////////////////////////////////////
  // This is just a template trick to avoid template specialization //
  // in EinsplineSetExtended.                                       //
  ////////////////////////////////////////////////////////////////////
  template<typename StorageType, int dim>  struct MultiOrbitalTraits{};

  template<> struct MultiOrbitalTraits<double,2>
  {  typedef multi_UBspline_2d_d SplineType;  };

  template<> struct MultiOrbitalTraits<double,3>
  {  typedef multi_UBspline_3d_d SplineType;  };

  template<> struct MultiOrbitalTraits<complex<double>,2>
  {  typedef multi_UBspline_2d_z SplineType;  };

  template<> struct MultiOrbitalTraits<complex<double>,3>
  {  typedef multi_UBspline_3d_z SplineType;  };



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
    typedef typename OrbitalSetTraits<StorageType>::ValueVector_t StorageValueVector_t;
    typedef typename OrbitalSetTraits<StorageType>::GradVector_t  StorageGradVector_t;
    typedef typename OrbitalSetTraits<StorageType>::HessVector_t  StorageHessVector_t;
    typedef Vector<double>                                        RealValueVector_t;
    typedef Vector<complex<double> >                              ComplexValueVector_t;
    typedef Vector<TinyVector<double,OHMMS_DIM> >                 RealGradVector_t;
    typedef Vector<TinyVector<complex<double>,OHMMS_DIM> >        ComplexGradVector_t;
    typedef Vector<Tensor<double,OHMMS_DIM> >                     RealHessVector_t;
    typedef Matrix<Tensor<double,OHMMS_DIM> >                     RealHessMatrix_t;
    typedef Vector<Tensor<complex<double>,OHMMS_DIM> >            ComplexHessVector_t;
    typedef Matrix<Tensor<complex<double>,OHMMS_DIM> >            ComplexHessMatrix_t;
    typedef Matrix<double>                                        RealValueMatrix_t;
    typedef Matrix<complex<double> >                              ComplexValueMatrix_t;
    typedef Matrix<TinyVector<double,OHMMS_DIM> >                 RealGradMatrix_t;
    typedef Matrix<TinyVector<complex<double>,OHMMS_DIM> >        ComplexGradMatrix_t;
//     typedef typename OrbitalSetTraits<ReturnType >::ValueVector_t ReturnValueVector_t;
//     typedef typename OrbitalSetTraits<ReturnType >::GradVector_t  ReturnGradVector_t;
//     typedef typename OrbitalSetTraits<ReturnType >::HessVector_t  ReturnHessVector_t;

//     typedef typename OrbitalSetTraits<ReturnType >::ValueMatrix_t ReturnValueMatrix_t;
//     typedef typename OrbitalSetTraits<ReturnType >::GradMatrix_t  ReturnGradMatrix_t;
//     typedef typename OrbitalSetTraits<ReturnType >::HessMatrix_t  ReturnHessMatrix_t;
       
    /////////////////////////////
    /// Orbital storage object //
    /////////////////////////////
    SplineType *MultiSpline;

    //////////////////////////////////////
    // Radial/Ylm orbitals around atoms //
    //////////////////////////////////////
    vector<AtomicOrbital<StorageType> > AtomicOrbitals;


    // First-order derivative w.r.t. the ion positions
    vector<TinyVector<SplineType*,OHMMS_DIM> > FirstOrderSplines;
    
    // Temporary storage for Eispline calls
    StorageValueVector_t StorageValueVector, StorageLaplVector;
    StorageGradVector_t  StorageGradVector;
    StorageHessVector_t  StorageHessVector;
    // Temporary storage used when blending functions        
    StorageValueVector_t BlendValueVector, BlendLaplVector;   
    StorageGradVector_t BlendGradVector;
        
    // True if we should unpack this orbital into two copies
    vector<bool>         MakeTwoCopies;
    // k-points for each orbital
    Vector<TinyVector<double,OHMMS_DIM> > kPoints;

    ///////////////////
    // Phase factors //
    ///////////////////
    Vector<double> phase;
    Vector<complex<double> > eikr;
    inline void computePhaseFactors(TinyVector<double,OHMMS_DIM> r);
    // For running at half G-vectors with real orbitals;  
    // 0 if the twist is zero, 1 if the twist is G/2.
    TinyVector<int,OHMMS_DIM> HalfG;

    ////////////
    // Timers //
    ////////////
    NewTimer ValueTimer, VGLTimer, VGLMatTimer;
    NewTimer EinsplineTimer;
  public:
    void registerTimers();

    // Real return values
    void evaluate(const ParticleSet& P, int iat, RealValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, RealValueVector_t& psi, 
		  RealGradVector_t& dpsi, RealValueVector_t& d2psi);
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
		  RealValueMatrix_t& psi, RealGradMatrix_t& dpsi, 
		  RealValueMatrix_t& d2psi);

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
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
		  ComplexValueMatrix_t& psi, ComplexGradMatrix_t& dpsi, 
		  ComplexValueMatrix_t& d2psi);
    
    void resetParameters(const opt_variables_type& active);
    void resetTargetParticleSet(ParticleSet& e);
    void setOrbitalSetSize(int norbs);
    string Type();
    
    SPOSetBase* makeClone() const;
    
    EinsplineSetExtended() : 
      ValueTimer  ("EinsplineSetExtended::ValueOnly"),
      VGLTimer    ("EinsplineSetExtended::VGL"),
      VGLMatTimer ("EinsplineSetExtended::VGLMatrix"),
      EinsplineTimer("libeinspline")
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

  template<typename StorageType>
  inline void EinsplineSetExtended<StorageType>::computePhaseFactors(TinyVector<double,OHMMS_DIM> r)
  {
    for (int i=0; i<kPoints.size(); i++) phase[i] = -dot(r, kPoints[i]);
    eval_e2iphi(phase,eikr);
//#ifdef HAVE_MKL
//    for (int i=0; i<kPoints.size(); i++) 
//      phase[i] = -dot(r, kPoints[i]);
//    vzCIS(OrbitalSetSize, phase, (double*)eikr.data());
//#else
//    double s, c;
//    for (int i=0; i<kPoints.size(); i++) {
//      phase[i] = -dot(r, kPoints[i]);
//      sincos (phase[i], &s, &c);
//      eikr[i] = complex<double>(c,s);
//    }
//#endif
  }
  

  
}
#endif
