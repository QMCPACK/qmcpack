//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BSPLINE_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <sstream>

namespace qmcplusplus {

  template<class T>
  struct BsplineFunctor: public OptimizableFunctorBase<T> {
    ///typedef of real values
    typedef typename OptimizableFunctorBase<T>::real_type real_type;
    typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
    
    std::vector<double> SplineCoefs;
    std::vector<double> Parameters;
    std::vector<std::string> ParameterNames;
    double Rcut, DeltaR, DeltaRInv;
    static const double A[16], dA[16], d2A[16];
    double CuspValue;
    int NumParams;

    ///constructor
    BsplineFunctor() : NumParams(0), Rcut(0.0)
    {
      reset();
    }

    void resize(int n) 
    { 
      ParameterNames.resize(n);
      Parameters.resize(n);
      SplineCoefs.resize(n+4);
      for (int i=0; i< n; i++) {
	std::stringstream sstr;
	sstr << "P" << i;
	ParameterNames[i] = sstr.str();
      }
    }
    
    void reset() 
    {

    }
    
    inline real_type evaluate(real_type r) {
      r *= DeltaRInv;

      double ipart, t;
      t = modf (r, &ipart);
      int i = (int) ipart;
      
      double tp[4];
      tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;

      return 
	(SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
	 SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
	 SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
	 SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));
		
    }

    inline real_type 
    evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
      r *= DeltaRInv;
      double ipart, t;
      t = modf (r, &ipart);
      int i = (int) ipart;
      
      double tp[4];
      tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;

      dudr = DeltaRInv * 
	(SplineCoefs[i+0]*(dA[ 0]*tp[0] + dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3])+
	 SplineCoefs[i+1]*(dA[ 4]*tp[0] + dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3])+
	 SplineCoefs[i+2]*(dA[ 8]*tp[0] + dA[ 9]*tp[1] + dA[10]*tp[2] + dA[11]*tp[3])+
	 SplineCoefs[i+3]*(dA[12]*tp[0] + dA[13]*tp[1] + dA[14]*tp[2] + dA[15]*tp[3]));
      d2udr2 = DeltaRInv * DeltaRInv *
	(SplineCoefs[i+0]*(d2A[ 0]*tp[0] + d2A[ 1]*tp[1] + d2A[ 2]*tp[2] + d2A[ 3]*tp[3])+
	 SplineCoefs[i+1]*(d2A[ 4]*tp[0] + d2A[ 5]*tp[1] + d2A[ 6]*tp[2] + d2A[ 7]*tp[3])+
	 SplineCoefs[i+2]*(d2A[ 8]*tp[0] + d2A[ 9]*tp[1] + d2A[10]*tp[2] + d2A[11]*tp[3])+
	 SplineCoefs[i+3]*(d2A[12]*tp[0] + d2A[13]*tp[1] + d2A[14]*tp[2] + d2A[15]*tp[3]));
      return 
	(SplineCoefs[i+0]*(A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3])+
	 SplineCoefs[i+1]*(A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3])+
	 SplineCoefs[i+2]*(A[ 8]*tp[0] + A[ 9]*tp[1] + A[10]*tp[2] + A[11]*tp[3])+
	 SplineCoefs[i+3]*(A[12]*tp[0] + A[13]*tp[1] + A[14]*tp[2] + A[15]*tp[3]));

    }

    inline real_type f(real_type r) {
      if(r>Rcut) return 0.0;
      return evaluate (r);
    }
    inline real_type df(real_type r) {
      if(r>Rcut) return 0.0;
      double du, d2u;
      evaluate (r, du, d2u);
      return du;
    }
    
    bool put(xmlNodePtr cur) 
    {
      OhmmsAttributeSet rAttrib;
      rAttrib.add(NumParams,"n");
      rAttrib.add(CuspValue, "cusp");
      rAttrib.add(Rcut, "Rcut");
      rAttrib.put(cur);
      if (NumParams == 0) {
	app_error() << "You must specify a positive number of parameters for the "
		    << "one-body spline jastrow function.\n";
	abort();
      }
      resize(NumParams);
      return true;
    }
    
    void addOptimizables(OptimizableSetType& vlist)
    {
      for (int i=0; i<ParameterNames.size(); i++)
	vlist[ParameterNames[i]] = Parameters[i];
    }
    
    /** reset the internal variables.
     *
     * USE_resetParameters
     */
    void resetParameters(OptimizableSetType& optVariables) 
    {
      for (int i=0; i<ParameterNames.size(); i++) {
	typename OptimizableSetType::iterator it(optVariables.find(ParameterNames[i]));
	if(it != optVariables.end()) 
	  Parameters[i] = it->second;
      }
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: BsplineConstraints.h 1691 2007-02-01 21:51:50Z jnkim $ 
 ***************************************************************************/
