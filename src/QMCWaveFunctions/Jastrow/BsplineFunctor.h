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
#include <cstdio>

namespace qmcplusplus {

  template<class T>
  struct BsplineFunctor: public OptimizableFunctorBase<T> {
    ///typedef of real values
    typedef typename OptimizableFunctorBase<T>::real_type real_type;
    typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
    
    std::string elementType;
    std::vector<double> SplineCoefs;
    std::vector<double> Parameters;
    std::vector<std::string> ParameterNames;
    double Rcut, DeltaR, DeltaRInv;
    static const real_type A[16], dA[16], d2A[16];
    double CuspValue;
    int NumParams;

    ///constructor
    BsplineFunctor() : NumParams(0), Rcut(0.0), CuspValue(0.0)
    {
    }

    void resize(int n) 
    { 
      NumParams = n;
      int numCoefs = NumParams + 4;
      int numKnots = numCoefs - 2;
      DeltaR = Rcut / (double)(numKnots - 1);
      DeltaRInv = 1.0/DeltaR;

      Parameters.resize(n);
      SplineCoefs.resize(numCoefs);
    }
    
    void reset() 
    {
      for (int i=0; i<SplineCoefs.size(); i++)
	SplineCoefs[i] = 0.0;

      // Ensure that cusp conditions is satsified at the origin
      SplineCoefs[1] = Parameters[0];
      SplineCoefs[2] = Parameters[1];
      SplineCoefs[0] = Parameters[1] - 2.0*DeltaR * CuspValue;
      for (int i=2; i<Parameters.size(); i++)
	SplineCoefs[i+1] = Parameters[i];
      FILE *fout = fopen ((elementType + ".dat").c_str(), "w");
      for (double r=0.0; r<Rcut; r+=0.01)
	fprintf (fout, "%1.10e %1.10e\n", r, evaluate(r));
      fclose (fout);
//       cerr << "SplineCoefs = ";
//       for (int i=0; i<SplineCoefs.size(); i++)
// 	cerr << SplineCoefs[i] << " ";
//       cerr << endl;
    }
    
    inline real_type evaluate(real_type r) {
      if (r >= Rcut)
	return 0.0;
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
      if (r >= Rcut) {
	dudr = d2udr2 = 0.0;
	return 0.0;
      }
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
      CuspValue = 0.0;
      NumParams = 0;
      Rcut = 0.0;
      OhmmsAttributeSet rAttrib;
      rAttrib.add(elementType, "elementType");
      rAttrib.add(NumParams,   "size");
      rAttrib.add(CuspValue,   "cusp");
      rAttrib.add(Rcut,        "rcut");
      rAttrib.put(cur);

      if (NumParams == 0) {
	app_error() << "You must specify a positive number of parameters for the "
		    << "one-body spline jastrow function.\n";
	abort();
      }
      resize (NumParams);
      // Now read coefficents
      xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
      bool haveCoefs = false;
      while (xmlCoefs != NULL) {
	if ((char*)xmlCoefs->name == std::string("coefficients")) {
	  string type, id;
	  OhmmsAttributeSet cAttrib;
	  cAttrib.add(id, "id");
	  cAttrib.add(type, "type");
	  cAttrib.put(xmlCoefs);
	  
	  if (type != "Array") {
	    app_error() << "Unknown correlation type """ << type 
			<< """ in BsplineFunctor.\n"
			<< "Resetting to ""Array"".\n";
	    xmlNewProp (xmlCoefs, (const xmlChar*) "type", 
			(const xmlChar*) "Array");
	  }
	  
	  haveCoefs = true;
	  std::stringstream sstr;
	  sstr << (char*)xmlCoefs->xmlChildrenNode->content;
	  for (int i=0; i<NumParams; i++)
	    sstr >> Parameters[i];
	  
	  // Setup parameter names
	  for (int i=0; i< NumParams; i++) {
	    std::stringstream sstr;
	    sstr << id << "_" << i;
	    ParameterNames.push_back(sstr.str());
	  }

	  cerr << "Parameter     Name      Value\n";
	  for (int i=0; i<ParameterNames.size(); i++)
	    cerr << "    " << i << "         " << ParameterNames[i] 
		 << "       " << Parameters[i] << endl;
	}
	xmlCoefs = xmlCoefs->next;
      }

      reset();
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
	if(it != optVariables.end()) {
// 	  cerr << "Resetting " << ParameterNames[i] << " to " 
// 	       << it->second << endl;
	  Parameters[i] = it->second;
	}
	reset();
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
