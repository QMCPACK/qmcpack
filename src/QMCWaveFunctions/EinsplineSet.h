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
#include <einspline/multi_bspline_structs.h>


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
    UnitCellType SuperLattice, PrimLattice, PrimLatticeInv;
    /// The "Twist" variables are in reduced coords, i.e. from 0 to1.
    /// The "k" variables are in Cartesian coordinates.
    PosType TwistVector, kVector;
    /// This stores which "true" twist vector this clone is using.
    /// "True" indicates the physical twist angle after untiling
    int TwistNum;
    /// metric tensor to handle generic unitcell
    Tensor<RealType,OHMMS_DIM> GGt;
    
    /////////////////////
    // Orbital storage //
    /////////////////////
    /// Store the orbital objects.  Using template class allows us to
    /// avoid making separate real and complex versions of this class.
    //std::vector<EinsplineOrb<ValueType,OHMMS_DIM>*> Orbitals;
    std::vector<EinsplineOrb<complex<double>,OHMMS_DIM>*> Orbitals;
    
    
  public:  
    UnitCellType GetLattice();

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, 
		  ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
		  ValueMatrix_t& d2logdet);
    
    void resetParameters(VarRegistry<RealType>& vlist);
    void resetTargetParticleSet(ParticleSet& e);
    void setOrbitalSetSize(int norbs);
    EinsplineSet() :  TwistNum(0)
    {
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
  template<typename StorageType, typename ReturnType>
  class EinsplineSetExtended : public SPOSetBase
  {
        friend class EinsplineSetBuilder;
  protected:
    //////////////////////
    // Type definitions //
    //////////////////////
    typedef CrystalLattice<RealType,OHMMS_DIM> UnitCellType;
    typedef typename MultiOrbitalTraits<StorageType,OHMMS_DIM>::SplineType SplineType; 
    typedef typename OrbitalSetTraits<StorageType>::ValueVector_t StorageValueVector_t;
    typedef typename OrbitalSetTraits<StorageType>::GradVector_t  StorageGradVector_t;
    typedef typename OrbitalSetTraits<StorageType>::HessVector_t  StorageHessVector_t;
    typedef typename OrbitalSetTraits<ReturnType >::ValueVector_t ReturnValueVector_t;
    typedef typename OrbitalSetTraits<ReturnType >::GradVector_t  ReturnGradVector_t;
    typedef typename OrbitalSetTraits<ReturnType >::HessVector_t  ReturnHessVector_t;

    typedef typename OrbitalSetTraits<ReturnType >::ValueMatrix_t ReturnValueMatrix_t;
    typedef typename OrbitalSetTraits<ReturnType >::GradMatrix_t  ReturnGradMatrix_t;
    typedef typename OrbitalSetTraits<ReturnType >::HessMatrix_t  ReturnHessMatrix_t;
   
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
    UnitCellType SuperLattice, PrimLattice, PrimLatticeInv;
    /// The "Twist" variables are in reduced coords, i.e. from 0 to1.
    /// The "k" variables are in Cartesian coordinates.
    PosType TwistVector, kVector;
    /// This stores which "true" twist vector this clone is using.
    /// "True" indicates the physical twist angle after untiling
    int TwistNum;
    /// metric tensor to handle generic unitcell
    Tensor<RealType,OHMMS_DIM> GGt;
    
    /////////////////////////////
    /// Orbital storage object //
    /////////////////////////////
    SplineType *MultiSpline;
    // Temporary storage for Eispline calls
    StorageValueVector_t StorageValueVector;
    StorageGradVector_t  StorageGradVector;
    StorageHessVector_t  StorageHessVector;
    // k-points for each orbital
    Vector<TinyVector<double,OHMMS_DIM> > kPoints;
   
  public:
    UnitCellType GetLattice();

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, 
		  GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
		  ValueMatrix_t& d2logdet);
    
    void resetParameters(VarRegistry<RealType>& vlist);
    void resetTargetParticleSet(ParticleSet& e);
    void setOrbitalSetSize(int norbs);

    EinsplineSetExtended() : TwistNum(0)
    {
    }
  };

  
}
#endif
