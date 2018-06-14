//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_CUBICFUNCTORSFORJASTROW_H
#define QMCPLUSPLUS_CUBICFUNCTORSFORJASTROW_H
#include "Numerics/OneDimGridBase.h"
#include "Numerics/CubicSpline.h"
#include "Numerics/CubicBspline.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

/** A numerical functor
 *
 * implements interfaces to be used for Jastrow functions and replaces  CubicBsplineSingle.
 * Template parameters
 * - RT real data type
 * - FNOUT final numerical functor
 *  An example is in CBSOBuilder.h which uses CubicBspline
 *  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
 *  typedef CubicSplineSingle<RealType,SplineEngineType> RadialOrbitalType;
 */
template <typename RT, typename FNOUT>
struct CubicSplineSingle: public OptimizableFunctorBase
{

  ///typedef for the value_type
  typedef RT value_type;
  ///typedef of the source functor
  typedef OptimizableFunctorBase FNIN;
  ///typedef for the grid
  typedef OneDimGridBase<real_type> grid_type;

  // mmorales: until I figure out how to go around this
  int NumParams;
  FNIN *InFunc;
  FNOUT OutFunc;
  int NumGridPoints;
  real_type Rmax;
  real_type GridDelta;
  real_type Y;
  real_type dY;
  real_type d2Y;
  ///constructor
  CubicSplineSingle(): InFunc(0) { }

  CubicSplineSingle(const CubicSplineSingle& old):
    OptimizableFunctorBase(old),
    NumGridPoints(old.NumGridPoints),
    Rmax(old.Rmax),
    GridDelta(old.GridDelta)
  {
    NumParams=0;
    if(old.InFunc)
    {
      initialize(old.InFunc->makeClone(),old.Rmax,old.NumGridPoints);
    }
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new CubicSplineSingle<RT,FNOUT>(*this);
  }

  ///constructor with arguments
  CubicSplineSingle(FNIN* in_, grid_type* agrid): InFunc(in_)
  {
    initialize(in_,agrid);
  }
  ///constructor with arguments
  CubicSplineSingle(FNIN* in_, real_type rc, int npts):InFunc(in_)
  {
    initialize(in_,rc,npts);
  }
  ///set the input, analytic function
  void setInFunc(FNIN* in_)
  {
    InFunc=in_;
  }

  void reportStatus(std::ostream& os)
  {
    //myVars.print(os);
  }

  bool isOptimizable()
  {
    return false;
  }

  /** evaluate everything: value, first, second and third derivatives
   */
  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type &d3udr3)
  {
    std::cerr << "Third derivative not implemented for CubicSplineSingle.\n";
    return OutFunc.splint(r,dudr,d2udr2);
  }


  /** evaluate everything: value, first and second derivaties
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    return OutFunc.splint(r,dudr,d2udr2);
  }

  /** evaluate value only
   */
  inline real_type evaluate(real_type r)
  {
    return OutFunc.splint(r);
  }

  /** evaluate value only
   *
   * Function required for SphericalBasisSet
   */
  inline real_type evaluate(real_type r, real_type rinv)
  {
    return Y=OutFunc.splint(r);
  }

  /** evaluate everything: value, first and second derivaties
   *
   * Function required for SphericalBasisSet
   */
  inline real_type evaluateAll(real_type r, real_type rinv)
  {
    return Y=OutFunc.splint(r,dY,d2Y);
  }

  /** implement the virtual function of OptimizableFunctorBase */
  inline real_type f(real_type r)
  {
    return OutFunc.splint(r);
  }

  /** implement the virtual function of OptimizableFunctorBase  */
  inline real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    OutFunc.splint(r,dudr,d2udr2);
    return dudr;
  }

  bool put(xmlNodePtr cur)
  {
    bool s=false;
    if(InFunc)
      s=InFunc->put(cur);
    return s;
  }

  void checkInVariables(opt_variables_type& active)
  {
    if(InFunc)
      InFunc->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    if(InFunc)
      InFunc->checkOutVariables(active);
  }

  ///reset the input/output function
  void resetParameters(const opt_variables_type& active)
  {
    if(InFunc)
    {
      InFunc->resetParameters(active);
      reset();
    }
  }

  void print(std::ostream& os)
  {
    real_type r=0;
    for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
      os << r << " " << OutFunc.splint(r) << std::endl;
  }

  ///set the input, analytic function
  void initialize(FNIN* in_, grid_type* agrid)
  {
    initialize(in_,agrid->rmax(),agrid->size());
  }

  void initialize(FNIN* in_, real_type rmax, int npts)
  {
    InFunc=in_;
    Rmax=rmax;
    NumGridPoints=npts;
    GridDelta=Rmax/static_cast<real_type>(NumGridPoints-1);
    reset();
  }

  void reset()
  {
    if(InFunc)
    {
      typename FNOUT::container_type datain(NumGridPoints);
      real_type r=0;
      for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
      {
        datain[i] = InFunc->f(r);
      }
      OutFunc.Init(0.0,Rmax,datain,true,InFunc->df(0.0),0.0);
    }
    else
    {
      APP_ABORT("CubicSplineSingle::reset has no input functor");
    }
  }
};


template <typename RT>
struct CubicSplineBasisSet: public OptimizableFunctorBase
{

  ///typedef of the source functor
  typedef OptimizableFunctorBase FNIN;
  ///typedef for the argument
  typedef CubicBspline<RT,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> FNOUT;
  ///typedef for the grid
  typedef OneDimGridBase<real_type> grid_type;

  FNIN *InFunc;
  FNOUT *OutFunc;
  int NumGridPoints;
  real_type Rmax;
  real_type GridDelta;

  ///constructor
  CubicSplineBasisSet(): InFunc(0), OutFunc(0) { }
  ///constructor with arguments
  CubicSplineBasisSet(FNIN* in_, grid_type* agrid)
  {
    initialize(in_,agrid);
  }
  ///set the input, analytic function
  void setInFunc(FNIN* in_)
  {
    InFunc=in_;
  }
  ///set the output numerical function
  void setOutFunc(FNOUT* out_)
  {
    OutFunc=out_;
  }
  ///reset the input/output function
  void resetParameters(const opt_variables_type& active)
  {
    if(!InFunc)
      APP_ABORT("CubicSplineBasisSet::resetParameters failed due to null input function ");
    InFunc->resetParameters(active);
    reset();
  }

  void reset()
  {
    if(!OutFunc)
      OutFunc = new FNOUT;
    typename FNOUT::container_type datain(NumGridPoints);
    real_type r=0;
    for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
      datain[i] = InFunc->f(r);
    OutFunc->Init(0.0,Rmax,datain,true,InFunc->df(0.0),0.0);
  }

  /** evaluate everything: value, first and second derivaties
  */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    return OutFunc->splint(r,dudr,d2udr2);
  }

  /** evaluate value only
  */
  inline real_type evaluate(real_type r)
  {
    return OutFunc->splint(r);
  }

  /** implement the virtual function of OptimizableFunctorBase */
  real_type f(real_type r)
  {
    return OutFunc->splint(r);
  }

  /** implement the virtual function of OptimizableFunctorBase  */
  real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    OutFunc->splint(r,dudr,d2udr2);
    return dudr;
  }

  bool put(xmlNodePtr cur)
  {
    return InFunc->put(cur);
  }

  void print(std::ostream& os)
  {
    real_type r=0;
    for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
      os << r << " " << OutFunc->splint(r) << std::endl;
  }

  ///set the input, analytic function
  void initialize(FNIN* in_, grid_type* agrid)
  {
    Rmax=agrid->rmax();
    NumGridPoints=agrid->size();
    GridDelta=Rmax/static_cast<real_type>(NumGridPoints-1);
    InFunc=in_;
    reset();
  }
};

//  /** A numerical functor using cubic bspline
//   *
//   * This class is replaced by CubicSplineSingle<RT,FNOUT>
//   */
//  template <typename RT>
//    struct CubicBsplineSingle: public OptimizableFunctorBase<RT> {
//
//      ///typedef of the source functor
//      typedef OptimizableFunctorBase<RT> FNIN;
//      ///typedef for the argument
//      typedef typename FNIN::real_type real_type;
//      ///typedef for OptimizableSetType
//      typedef typename FNIN::OptimizableSetType OptimizableSetType;
//      ///typedef for the argument
//      typedef CubicBspline<RT,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> FNOUT;
//      ///typedef for the grid
//      typedef OneDimGridBase<real_type> grid_type;
//
//      FNIN *InFunc;
//      FNOUT *OutFunc;
//      int NumGridPoints;
//      real_type Rmax;
//      real_type GridDelta;
//      real_type Y;
//      real_type dY;
//      real_type d2Y;
//
//
//      ///constructor
//      CubicBsplineSingle(): InFunc(0), OutFunc(0) { }
//      ///constructor with arguments
//      CubicBsplineSingle(FNIN* in_, grid_type* agrid): InFunc(0), OutFunc(0) {
//        initialize(in_,agrid);
//      }
//      ///constructor with arguments
//      CubicBsplineSingle(FNIN* in_, real_type rc, int npts):InFunc(0), OutFunc(0){
//        initialize(in_,rc,npts);
//      }
//      ///set the input, analytic function
//      void setInFunc(FNIN* in_) { InFunc=in_;}
//      ///set the output numerical function
//      void setOutFunc(FNOUT* out_) { OutFunc=out_;}
//
//      /** evaluate everything: value, first and second derivaties
//       */
//      inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//        return OutFunc->splint(r,dudr,d2udr2);
//      }
//
//      /** evaluate value only
//       */
//      inline real_type evaluate(real_type r) {
//        return OutFunc->splint(r);
//      }
//
//      /** evaluate value only
//       *
//       * Function required for SphericalBasisSet
//       */
//      inline real_type evaluate(real_type r, real_type rinv)
//      {
//        return Y=OutFunc->splint(r);
//      }
//
//      /** evaluate everything: value, first and second derivaties
//       *
//       * Function required for SphericalBasisSet
//       */
//      inline real_type evaluateAll(real_type r, real_type rinv)
//      {
//        return Y=OutFunc->splint(r,dY,d2Y);
//      }
//
//      /** implement the virtual function of OptimizableFunctorBase */
//      inline real_type f(real_type r)
//      {
//        return OutFunc->splint(r);
//      }
//
//      /** implement the virtual function of OptimizableFunctorBase  */
//      inline real_type df(real_type r) {
//        real_type dudr,d2udr2;
//        OutFunc->splint(r,dudr,d2udr2);
//        return dudr;
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        if(InFunc)
//          return InFunc->put(cur);
//        else
//          return false;
//      }
//
//      void addOptimizables(OptimizableSetType& vlist)
//      {
//        if(InFunc)
//          InFunc->addOptimizables(vlist);
//      }
//
//      ///reset the input/output function
//      void resetParameters(OptimizableSetType& optVariables)
//      {
//        if(!InFunc)
//        {
//          app_error() << "  CubicSplineJastrow::reset failed due to null input function " << std::endl;
//          OHMMS::Controller->abort();
//        }
//        InFunc->resetParameters(optVariables);
//        resetInternals();
//      }
//
//      void print(std::ostream& os) {
//        real_type r=0;
//        for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
//          os << r << " " << OutFunc->splint(r) << std::endl;
//      }
//
//      ///set the input, analytic function
//      void initialize(FNIN* in_, grid_type* agrid) {
//        initialize(in_,agrid->rmax(),agrid->size());
//      }
//      void initialize(FNIN* in_, real_type rmax, int npts)
//      {
//        InFunc=in_;
//        Rmax=rmax;
//        NumGridPoints=npts;
//        GridDelta=Rmax/static_cast<real_type>(NumGridPoints-1);
//        resetInternals();
//      }
//
//      void resetInternals()
//      {
//        if(!OutFunc) OutFunc = new FNOUT;
//        typename FNOUT::container_type datain(NumGridPoints);
//        real_type r=0;
//        for(int i=0; i<NumGridPoints; i++, r+=GridDelta)
//        {
//          datain[i] = InFunc->f(r);
//        }
//        OutFunc->Init(0.0,Rmax,datain,true,InFunc->df(0.0),0.0);
//      }
//
//    };
//
//  /** A numerical functor
//   *
//   * implements interfaces to be used for Jastrow functions
//   * - OneBodyJastrow<NumericalJastrow>
//   * - TwoBodyJastrow<NumericalJastrow>
//   */
//  template <class RT>
//    struct SplineJastrow: public OptimizableFunctorBase<RT> {
//
//      ///typedef of the source functor
//      typedef OptimizableFunctorBase<RT> FNIN;
//      ///typedef for the argument
//      typedef typename FNIN::real_type real_type;
//      ///typedef of the target functor
//      typedef OneDimCubicSpline<real_type,real_type>  FNOUT;
//
//      real_type Rmax;
//      FNIN *InFunc;
//      FNOUT *OutFunc;
//
//      ///constructor
//      SplineJastrow(): InFunc(0), OutFunc(0) { }
//      ///constructor with arguments
//      SplineJastrow(FNIN* in_, typename FNOUT::grid_type* agrid){
//        initialize(in_,agrid);
//      }
//      ///set the input, analytic function
//      void setInFunc(FNIN* in_) { InFunc=in_;}
//      ///set the output numerical function
//      void setOutFunc(FNOUT* out_) { OutFunc=out_;}
//      ///reset the input/output function
//      inline void reset() {
//        InFunc->reset();
//        //reference to the output functions grid
//        const typename FNOUT::grid_type& grid = OutFunc->grid();
//        //set cutoff function
//        int last=grid.size()-1;
//        for(int i=0; i<grid.size(); i++) {
//          (*OutFunc)(i) = InFunc->f(grid(i));
////	  std::cout << grid(i) << "   " << (*OutFunc)(i) << std::endl;
//        }
//	(*OutFunc)(last)=0.0;
//        //boundary conditions
//        real_type deriv1=InFunc->df(grid(0));
//        real_type deriv2=0.0;
//        OutFunc->spline(0,deriv1,last,deriv2);
//        Rmax=grid(last);
//      }
//
//      /** evaluate everything: value, first and second derivaties
//      */
//      inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//        return OutFunc->splint(r,dudr,d2udr2);
//      }
//
//      /** evaluate value only
//      */
//      inline real_type evaluate(real_type r) {
//        return OutFunc->splint(r);
//      }
//
//      /** implement the virtual function of OptimizableFunctorBase */
//      real_type f(real_type r) {
//        return OutFunc->splint(r);
//      }
//
//      /** implement the virtual function of OptimizableFunctorBase  */
//      real_type df(real_type r) {
//        real_type dudr,d2udr2;
//        OutFunc->splint(r,dudr,d2udr2);
//        return dudr;
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        return InFunc->put(cur);
//      }
//
//      void addOptimizables( VarRegistry<real_type>& vlist)
//      {
//        InFunc->addOptimizables(vlist);
//      }
//      void print(std::ostream& os) {
//        const typename FNOUT::grid_type& grid = OutFunc->grid();
//        for(int i=0; i<grid.size(); i++) {
//          std::cout << grid(i) << " " << (*OutFunc)(i) << std::endl;
//        }
//      }
//
//      ///set the input, analytic function
//      void initialize(FNIN* in_, typename FNOUT::grid_type* agrid) {
//        InFunc=in_;
//        setOutFunc(new FNOUT(agrid));
//        reset();
//      }
//    };
//
}
#endif
