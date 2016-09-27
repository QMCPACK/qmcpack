//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file EinsplineSetTemp.h
 * @brief Refactoring EinsplineSet to handle mixed data sets.
 *
 * Warning: heavy use of templates
 */
#ifndef QMCPLUSPLUS_EINSPLINE_SET_TEMPATE_H
#define QMCPLUSPLUS_EINSPLINE_SET_TEMPATE_H

#include <QMCWaveFunctions/OrbitalSetTraits.h>
#include <QMCWaveFunctions/Spline3D/twist_handler.hpp>
extern "C"
{
#include <einspline/bspline.h>
}
#include <einspline/multi_bspline_structs.h>

namespace qmcplusplus
{
  /** determine if EngT (e.g., einspline engine) handles real data or complex data
   *
   * Default is true and complex cases are specialized
   */
  template<typename EngT>
  struct is_real_bspline
  {
      static const bool value=true;
  };

  ///specialization for multi_UBspline_3d_z
  template<>
  struct is_real_bspline<multi_UBspline_3d_z>
  {
      static const bool value=false;
  };

  /////specialization for multi_UBspline_3d_z
  //template<>
  //struct is_real_bspline<multi_UBspline_3d_c>
  //{
  //    static const bool value=false;
  //};
  //
  /** struct to check if two types are the same
   */
  template<typename T1, typename T2>
    struct type_check
    {
      static const bool value=false;
    };

  template<typename T1>
    struct type_check<T1,T1>
    {
      static const bool value=true;
    };

  /** dummy traits class for bspline engine
   *
   * Should fail to instantiate invalid engines if the trait class is not implemented
   * The specialization provides
   * - DIM, physical dimension
   * - Real_t real data type
   * - Type_t value data type
   * - Spline_t einspline engine type
   */
  template<typename EngT>
  struct bspline_engine_traits
  {};

  /** specialization with multi_UBspline_3d_d
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_d>
  {
    enum {DIM=3};
    typedef double Real_t;
    typedef double Type_t;
    typedef multi_UBspline_3d_d Spline_t;
  };

  ///specialization with multi_UBspline_3d_z
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_z>
  {
    enum {DIM=3};
    typedef double Real_t;
    typedef std::complex<double> Type_t;
    typedef multi_UBspline_3d_z Spline_t;
  };

  /** bitset index for bspline_engine to assist builder
   *
   * (Gamma,HalfG)x(Real,Complex)x(Orthorombic,General)
   */
  enum {BSPLINE_TWIST_BIT, BSPLINE_TYPE_BIT, BSPLINE_CELL_BIT};


  /** generic declaration of bspline_engin
   *
   * Template parameters
   * - EngT engine which evaluates quantities
   * - IsGamma true, if the twist angle is the gamma point
   * - IsOrtho true, the cell is orthorombic
   * - IsReal true, if EngT handles real data
   * IsReal can be deduced from is_real_bsplne<EngT>::value
   * This implementation is to abort when the proper specialization is not found.
   */
  template<typename EngT, bool IsGamma, bool IsOrtho, bool IsReal>
    struct bspline_engine
    {
      enum {DIM=bspline_engine_traits<EngT>::DIM};

      /** specialization should provide this function to evaluate V
       * @param r Carteisan position
       * @param psi starting address of the values
       */
      template<typename ValT>
        inline void evaluate(const TinyVector<ValT,DIM>& r, ValT* restrict psi)
        {
          APP_ABORT("No proper specialization for bspline_engine<EngT,IsGamma,IsOrtho,IsReal>");
        }
      /** specialization should provide this function to evaluate VGL 
       * @param r Carteisan position
       * @param psi starting address of the values
       * @param dpsi starting address of the gradients
       * @param d2psi starting address of the laplacians
       */
      template<typename ValT>
      inline void evaluate_vgl(const TinyVector<ValT,DIM>& r
          , ValT* restrict psi , TinyVector<ValT,DIM>* restrict dpsi , ValT* restrict d2psi)
      { 
        APP_ABORT("No proper specialization for bspline_engine<EngT,IsGamma,IsOrtho,IsReal>");
      }
    };

  /** specialization for complex bspline
   *
   * IsGamma is not used and specialization for IsOrtho is not implemented.
   * When IsReal==false, evaluate and evaluate_vgl are implemented
   * for real and complex arguments.
   */
  template<typename EngT, bool IsGamma, bool IsOrtho>
    struct bspline_engine<EngT,IsGamma,IsOrtho,false>
    : public twist_handler<typename bspline_engine_traits<EngT>::Real_t
      , bspline_engine_traits<EngT>::DIM>
  {
    ///typedef for the data type handled by the engine
    typedef typename bspline_engine_traits<EngT>::Type_t Type_t;
    ///typedef for the engine
    typedef typename bspline_engine_traits<EngT>::Spline_t Spline_t;
    /// Orbital storage objec
    Spline_t *MultiSpline;
    ///temporary values
    typename OrbitalSetTraits<Type_t>::ValueVector_t myPsi;
    ///temporary gradients
    typename OrbitalSetTraits<Type_t>::GradVector_t myGrad;
    ///Hessians are stored
    typename OrbitalSetTraits<Type_t>::HessVector_t myHess;

    enum {DIM=bspline_engine_traits<EngT>::DIM};

    /** value-only evaluation at r
     * @param r current position
     * @param psi starting address of the real values
     */
    template<typename ValT>
      inline void evaluate(const TinyVector<ValT,DIM>& r, ValT* restrict psi)
      {
        std::cout << "bspline_engine<EngT,IsGamma,IsOrtho,false>::evaluate, real" << std::endl;
//        PosType ru(PrimLattice.toUnit(r));
//        for (int i=0; i<DIM; i++) ru[i] -= std::floor (ru[i]);
//        bspline_evaluate(MultiSpline, ru, myPsi.data());
//        Real_t s,c;
//        for (int j=0,i=0; j<N; j++) 
//        {
//          sincos(-dot(r,kPoints[j]),&s,&c);
//          //assign the real value
//          psi[i++] = ValT(myPsi[j].real()*c-myPsi[j].imag()*s);
//          if (MakeTwoCopies[j]) 
//          {
//            //assign the imag value
//            psi[i++] = ValT(myPsi[j].imag()*c+myPsi[j].real()*s);
//          }
//        }
      }

    /** value-only evaluation at r
     * @param r current position
     * @param psi starting address of the complex values
     */
    template<typename ValT>
      inline void evaluate(const TinyVector<ValT,DIM>& r, std::complex<ValT>* restrict psi)
      {
        std::cout << "bspline_engine<EngT,IsGamma,IsOrtho,false>::evaluate, complex" << std::endl;
//        PosType ru(PrimLattice.toUnit(r));
//        for (int i=0; i<DIM; i++) ru[i] -= std::floor (ru[i]);
//        bspline_evaluate(MultiSpline,ru,myPsi.data());
//        Real_t s,c;
//        for (int i=0; i<N; i++) 
//        {
//          sincos(-dot(r,kPoints[j]), &s, &c);
//          psi[i]=std::complex<ValT>(c*myPsi[i].real()-s*myPsi[i].imag(),c*myPsi[i].imag()+s*myPsi.real());
//        }
      }

    template<typename ValT>
      inline void evaluate(const TinyVector<ValT,DIM>& r 
          , ValT* restrict psi , TinyVector<ValT,DIM>* restrict dpsi , ValT* restrict d2psi)
      {
        std::cout << "bspline_engine_vgl<EngT,IsGamma,IsOrtho,false>::evaluate, real" << std::endl;
//        PosType ru(PrimLattice.toUnit(r));
//        for (int i=0; i<DIM; i++) ru[i] -= std::floor (ru[i]);
//        bspline_evaluate(MultiSpline,ru,myPsi.data(),myGrad.data(),myHess.data());
//        Real_t s,c;
//        for (int i=0; i<N; i++) 
//        {
//          sincos(-dot(r,kPoints[j]), &s, &c);
//          Type_t e_mikr(c,s);
//          ///do the right thing
//        }
      }

    template<typename ValT>
      inline void evaluate(const TinyVector<ValT,DIM>& r
          , std::complex<ValT>* restrict psi, TinyVector<std::complex<ValT>,DIM>* restrict dpsi , std::complex<ValT>* restrict d2psi)
      {
        std::cout << "bspline_engine_vgl<EngT,IsGamma,IsOrtho,false>::evaluate, complex" << std::endl;
//        PosType ru(PrimLattice.toUnit(r));
//        for (int i=0; i<DIM; i++) ru[i] -= std::floor (ru[i]);
//        bspline_evaluate(MultiSpline,ru,myPsi.data(),myGrad.data(),myHess.data());
//        Real_t s,c;
//        for (int i=0; i<N; i++) 
//        {
//          sincos(-dot(r,kPoints[j]), &s, &c);
//          Type_t e_mikr(c,s);
//          ///do the right thing
//        }
      }
  };

  /** specialization for EngT for real storage, HalfG and general cell 
   *
   * The data type and dimention for twist_handler are deduced from the trait class
   */
  template<typename EngT,bool IsGamma, bool IsOrtho>
    struct bspline_engine<EngT,IsGamma,IsOrtho,true>
    : public twist_handler<typename bspline_engine_traits<EngT>::Real_t, bspline_engine_traits<EngT>::DIM>
    {
      ///typedef for the data type handled by the engine
      typedef typename bspline_engine_traits<EngT>::Type_t Type_t;
      ///typedef for the engine
      typedef typename bspline_engine_traits<EngT>::Spline_t Spline_t;
      /// Orbital storage objec
      Spline_t *MultiSpline;
      ///temporary values
      typename OrbitalSetTraits<Type_t>::ValueVector_t myPsi;
      ///temporary gradients
      typename OrbitalSetTraits<Type_t>::GradVector_t myGrad;
      ///Hessians are stored
      typename OrbitalSetTraits<Type_t>::HessVector_t myHess;

      enum {DIM=bspline_engine_traits<EngT>::DIM};
      //using PrimLattice;
      //using HalfG;
      ////this is going to be fixed
      //int N;

      template<typename ValT>
      inline void evaluate(const TinyVector<ValT,3>& r, ValT* restrict psi)
      {
        std::cout << "bspline_engine<EngT,IsGamma,IsOrtho,IsReal==true>::evaluate called" << std::endl;
        //PosType ru(PrimLattice.toUnit(r));
        //int asign=0;
        //for (int i=0; i<OHMMS_DIM; i++)
        //{
        //  RealType img = std::floor(ru[i]);
        //  ru[i] -= img;
        //  asign += HalfG[i] * (int)img;
        //}
        //bspline_evaluate(MultiSpline,ru,myPsi.data());
        //if(asign&1)
        //  for(int i=0; i<N; ++i) convert(-myPsi[i],psi[i]);
        //else
        //  for(int i=0; i<N; ++i) convert(myPsi[i],psi[i]);
      }

      /** function to evaluate vgl at r
       * @param r Carteisan position
       * @param psi starting address of the values
       * @param dpsi starting address of the gradients
       * @param d2psi starting address of the laplacians
       */
      template<typename ValT>
      inline void evaluate_vgl(const TinyVector<ValT,DIM>& r
          , ValT* restrict psi , TinyVector<ValT,DIM>* restrict dpsi , ValT* restrict d2psi)
      {
        std::cout << "bspline_engine<EngT,IsGamma,IsOrtho,IsReal==true>::evaluate_vgl called" << std::endl;
        //PosType ru(PrimLattice.toUnit(r));
        //int asign=0;
        //for (int i=0; i<OHMMS_DIM; i++)
        //{
        //  RealType img = std::floor(ru[i]);
        //  ru[i] -= img;
        //  asign += HalfG[i] * (int)img;
        //}
        //bspline_evaluate(MultiSpline,ru,psi,dpsi,Hessian.data());
        //if(asign&1)
        //{
        //  for(int i=0; i<N; ++i) convert(-myVal[i],psi[i]);
        //  for(int i=0; i<N; ++i) convert(-dot(PrimLattice.G,dpsi[i]),dpsi[i]);
        //  for(int i=0; i<N; ++i) convert(-trace(Hessian[i],GGt),d2psi[i]);
        //}
        //else
        //{
        //  for(int i=0; i<N; ++i) convert(myVal[i],psi[i]);
        //  for(int i=0; i<N; ++i) convert(dot(PrimLattice.G,dpsi[i]),dpsi[i]);
        //  for(int i=0; i<N; ++i) convert(trace(Hessian[i],GGt),d2psi[i]);
        //}
      }
    };

  /** specialization for EngT for real storage, gamma and orthogonal cell
   *
   * The data type and dimention for twist_handler are deduced from the trait class
   */
  template<typename EngT>
    struct bspline_engine<EngT,true,true,true>
    : public twist_handler<typename bspline_engine_traits<EngT>::Real_t, bspline_engine_traits<EngT>::DIM>
    {
      ///typedef for the data type handled by the engine
      typedef typename bspline_engine_traits<EngT>::Type_t Type_t;
      ///typedef for the engine
      typedef typename bspline_engine_traits<EngT>::Spline_t Spline_t;
      /// Orbital storage objec
      Spline_t *MultiSpline;
      ///temporary values
      typename OrbitalSetTraits<Type_t>::ValueVector_t myPsi;
      ///temporary gradients
      typename OrbitalSetTraits<Type_t>::GradVector_t myGrad;
      ///temporary laplacian
      typename OrbitalSetTraits<Type_t>::ValueVector_t myLap;
      ///
      enum {DIM=bspline_engine_traits<EngT>::DIM};
      ////this is going to be fixed
      //int N;

      template<typename ValT>
      inline void evaluate(const TinyVector<ValT,DIM>& r, ValT* restrict psi)
      {
        //TinyVector<Type_t,DIM> ru(PrimLattice.toUnit(r));
        if(type_check<Type_t,ValT>::value)
          std::cout << "bspline_engine<EngT,true,true,true>::evaluate no conversion " << std::endl;
          //bspline_evaluate(MultiSpline,ru,psi);
        else
        {
          std::cout << "bspline_engine<EngT,true,true,true>::evaluate conversion different datatype" << std::endl;
          //bspline_evaluate(MultiSpline,ru,myPsi.data());
          //for(int i=0; i<N; ++i) convert(myPsi[i],psi[i]);
        }
      }

      template<typename ValT>
      inline void evaluate_vgl(const TinyVector<ValT,DIM>& r
          , ValT* restrict psi , TinyVector<ValT,DIM>* restrict dpsi , ValT* restrict d2psi)
      {
        //TinyVector<Type_t,DIM> ru(PrimLattice.toUnit(r));
        if(type_check<Type_t,ValT>::value)
          std::cout << "bspline_engine<EngT,true,true,true>::evaluate_vgl no conversion " << std::endl;
          //bspline_evaluate(MultiSpline,ru,psi,dpsi,d2psi);
        else
        {
          std::cout << "bspline_engine<EngT,true,true,true>::evaluate_vgl conversion different datatype" << std::endl;
          //bspline_evaluate(MultiSpline,ru,myPs.data(),myGrad.data(),myLap.data());
          //for(int i=0; i<N; ++i) convert(myPsi[i],psi[i]);
          //for(int i=0; i<N; ++i) convert(myGrad[i],dpsi[i]);
          //for(int i=0; i<N; ++i) convert(myLap[i],d2psi[i]);
        }
      }
    };

//  /** default is set to complex einspline */
//  template<typename EngT, bool IsGamma=false, bool IsOrtho=false, bool IsReal=false>
//    struct EinsplineSetTemp : public SPOSetBase
//    {
//      ///typedef derived from bspline_engine<EngT,ValT>
//      typedef bspline_engine<EngT,IsGamma,IsReal,IsOrtho> spline_engine_type;
//      ///object which handles the spline engine and aux storages
//      spline_engine_type* my_bspline;
//
//      /** constructor */
//      EinsplineSetTemp(spline_engine_type* a=0): my_bspline(a)
//      {
//        className="EinsplineSetTemp";
//      }
//
//      /** use default constructor */
//      SPOSetBase* makeClone() const
//      {
//        return EinsplineSetTemp(*this);
//      }
//
//      void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
//      {
//        my_spline->evaluate(P.R[iat],psi.data());
//      }
//
//      void evaluate(const ParticleSet& P, int iat
//          , ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
//      {
//        my_spline->evaluate_vgl(P.R[iat],psi.data(),dpsi.data(),d2psi.data());
//      }
//
//      void evaluate_notranspose(const ParticleSet& P, int first, int last
//          , ValueMatrix_t& psi, GradMatrix_t& dpsi, ValueMatrix_t& d2psi)
//      {
//        for(int iat=first,i=0;iat<last;++iat,++i)
//          my_spline->evaluate_vgl(P.R[iat],psi[i],dpsi[i],d2psi[i]);
//      }
//
//      void resetParameters(const opt_variables_type& active)
//      {}
//  };
}
#endif
