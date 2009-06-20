//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BSPLINE3D_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE3D_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/LinearFit.h"
#include <cstdio>

namespace qmcplusplus {

  struct PolynomialFunctor3D: public OptimizableFunctorBase {

    typedef real_type value_type;
    int N_eI, N_ee;
    Array<real_type,3> gamma;
    Array<int,3> index;
    vector<real_type> gamaVec;
    int NumConstraints, NumGamma;
    Matrix<real_type> ConstraintMatrix;
    std::vector<real_type> Parameters;
    std::vector<std::string> ParameterNames;
    std::string iSpecies, eSpecies1, eSpecies2;
    int ResetCount;
    // Order of continuity
    const int C;

    ///constructor
    PolynomialFunctor3D() : 
      NumParams_eI(0), NumParams_ee(0), ResetCount(0),
      cutoff_radius(0.0), C(3);
    { }

    OptimizableFunctorBase* makeClone() const 
    {
      return new PolynomialFunctor3D(*this);
    }

    void resize(int neI, int nee) 
    { 
      N_eI = neI;
      N_ee = nee;

      gamma.resize(N_eI+1, N_eI+1, N_ee+1);
      index.resize(N_eI+1, N_eI+1, N_ee+1);
      NumGamma = ((N_eI+1)*(N_eI+2)/2 * (N_ee+1));
      NumConstraints = (2*N_eI+1) + (N_eI+N_ee+1);
      int numParams = numGamma - numConstraints;
      Parameters.resize(numParams);
      GammaVec.resize(numGamma);
      ConstraintMatrix.resize(numGamma, numGamma);
      // Assign indices
      int n=0;
      for (int m=0; m<=N_eI; m++)
	for (int l=m; l<=N_eI; l++)
	  for (int n=0; n<=N_ee; n++) 
	    index(l,m,n) = n++;
      assert (n == numGamma);


      // Fill up contraint matrix
      ConstraintMatrix = 0.0;
      int k;
      // e-e no-cusp constraint
      for (k=0; k<=2*N_eI; k++) {
	for (int m=0; m<=N_eI; m++) {
	  int l = k - m;
	  int i = index(l,m,1);
	  if (l > m)
	    ConstraintMatrix(k,i) = 2.0;
	  else if (l == m)
	    ConstraintMatrix(k,i) = 1.0;
	}
      }
      // e-I no-cusp constraint
      for (int kp=0; kp<=N_eI+N_ee; kp++) {
	ConstraintMatrix(k+kp,index(0,0,kp)) = (real_type) C;
	ConstraintMatrix(k+kp,index(1,0,kp)) = -0.5*cutoff_radius;
	for (int l=1; l<=NeI+N_ee; l++) {
	  int n = kp - l;
	  if (n > 0) {
	    ConstraintMatrix(k+kp,index(l,0,n)) = C;
	    ConstraintMatrix(k+dp,index(l,1,n)) = -0.5*cutoff_radius;
	  }
	}
      }
      
      for (int i=NumConstraints; i<numGamma; i++)
	ConstraintMatrix(i,i) = 1.0;

      // Now, invert constraint matrix
      Invert(ConstraintMatrix, numGamma, numGamma);
    }
    
    void reset() 
    {
      gammaVec = 0.0;
      // Set constrained parameters
      for (int i=0; i<NumConstraints; i++)
	for (int j=0; j<Parameters.size(); j++)
	  gammaVec[i] += ConstraintMatrix(i,j+NumConstraints)*Parameters[i];

      // Set unconstrained parameters
      for (int i=0; i<Parameters.size(); i++)
	gammaVec[NumConstraints+i] = Parameters[i];

      // Set gamma matrix
      int n=0;
      for (int m=0; m<=N_eI; m++)
	for (int l=m; l<=N_eI; l++)
	  for (int n=0; n<=N_ee; n++) 
	    gamma(m,l,n) = gamma(l,m,n) = gammaVec[n++];

      // Now check that constraints have been satisfied
      // e-e constraints
      for (int k=0; k<=2*N_eI; k++) {
	real_type sum=0.0;
	for (int l=0; l<=k; l++) {
	  int m = k - l;
	  sum += gamma(l,m,1);
	}
	if (fabs(sum) > 1.0e-10) {
	  app_error() << "e-e constraint not satisfied in PolynomialFunctor3D:  k=" 
		      << k << "  sum=" << sum << endl;
	  abort();
	}
      }

      // e-I constraints
      sum = 0.0;
      for (int k=0; k<=N_ei+N_ee; k++) {
	for (int m=0; m<=k; m++) {
	  int n = k - m;
	  sum += C*gamma(0,m,n) - 0.5*cutoff_radius*gamma(1,m,n);
	}
	if (fabs(sum) > 1.0e-10) {
	  app_error() << "e-I constraint not satisfied in PolynomialFunctor3D:  k=" 
		      << k << "  sum=" << sum << endl;
	  abort();
	}
      }
    }

    inline real_type evaluate(real_type r_12, 
			      real_type r_1I,
			      real_type r_2I) 
    {
      if (r_1I >= 0.5*cutoff_radius || r_2I >= 0.5*cutoff_radius) 
        return 0.0;

      real_type val = 0.0;
      real_type r2l=1.0;
      for (int l=0; l<N_eI; l++) {
	real_type r2m=1.0; 
	for (int m=0; m<N_eI; m++) {
	  real_type r2n=1.0;
	  for (int n=0; n<N_ee; n++) {
	    val += gamma(l,m,n)*r2l*r2m*r2n;
	    r2n *= r_12;
	  }
	  r2m *= r_2I;
	}
	r2l *= r_1I;
      }
      for (int i=0; i<C; i++)
	val *= (r_1I - 0.5*cutoff_radius)*(r_2I - 0.5*cutoff_radius);
      return val;
    }


    inline real_type evaluate(real_type r_12, real_type r_1I, real_type r_2I,
			      TinyVector<real_type,3> &grad,
			      Tensor<real_type,3> &hess) 
    {
      if (r_12 >= cutoff_radius || r_1I >= 0.5*cutoff_radius ||
	  r_2I >= 0.5*cutoff_radius) {
	grad = 0.0;
	hess = 0.0;
        return 0.0;
      }
      real_type val = 0.0;
      real_type r2l=1.0;
      for (int l=0; l<N_eI; l++) {
	real_type r2m=1.0; 
	for (int m=0; m<N_eI; m++) {
	  real_type r2n=1.0;
	  for (int n=0; n<N_ee; n++) {
	    val += gamma(l,m,n)*r2l*r2m*r2n;
	    r2n *= r_12;
	  }
	  r2m *= r_2I;
	}
	r2l *= r_1I;
      }
      for (int i=0; i<C; i++)
	val *= (r_1I - 0.5*cutoff_radius)*(r_2I - 0.5*cutoff_radius);

      return val;
    }


    inline real_type evaluate(real_type r, real_type rinv) 
    {
      return 0.0;
    }


    inline bool
    evaluateDerivatives (real_type r, vector<TinyVector<real_type,3> >& derivs)
    {
    }

    inline real_type f(real_type r) {
      return 0.0;
    }
    inline real_type df(real_type r) {
      return 0.0;
    }
    
    bool put(xmlNodePtr cur) 
    {
      // ReportEngine PRE("PolynomialFunctor3D","put(xmlNodePtr)");
      // //CuspValue = -1.0e10;
      // NumParams_eI = NumParams_ee = 0;
      // cutoff_radius = 0.0;
      // OhmmsAttributeSet rAttrib;
      // rAttrib.add(NumParams_ee,   "esize");
      // rAttrib.add(NumParams_eI,   "isize");
      // rAttrib.add(cutoff_radius,  "rcut");
      // rAttrib.put(cur);

      // if (NumParams_eI == 0) 
      //   PRE.error("You must specify a positive number for \"isize\"",true);
      // if (NumParams_ee == 0) 
      //   PRE.error("You must specify a positive number for \"esize\"",true);

      // app_log() << " esize = " << NumParams_ee << " parameters " << endl;
      // app_log() << " isize = " << NumParams_eI << " parameters " << endl;
      // app_log() << " rcut = " << cutoff_radius << endl;

      // resize (NumParams_eI, NumParams_ee);


      // // Now read coefficents
      // xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
      // while (xmlCoefs != NULL) 
      // {
      //   string cname((const char*)xmlCoefs->name);
      //   if (cname == "coefficients") 
      //   {
      //     string type("0"), id("0");
      //     OhmmsAttributeSet cAttrib;
      //     cAttrib.add(id, "id");
      //     cAttrib.add(type, "type");
      //     cAttrib.put(xmlCoefs);

      //     if (type != "Array") 
      //     {
      //       PRE.error( "Unknown correlation type " + type + 
      // 		       " in PolynomialFunctor3D." + "Resetting to \"Array\"");
      //       xmlNewProp (xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
      //     }

      // 	  vector<real_type> params;
      // 	  putContent(params, xmlCoefs);
      // 	  if (params.size() == Parameters.size()) 
      // 	    Parameters = params;
      // 	  else {
      // 	    app_error() << "Expected " << Parameters.size() << " parameters,"
      // 			<< " but found only " << params.size()
      // 			<< " in PolynomialFunctor3D.\n";
      // 	    abort();
      // 	  }
	 
	  
      //     // Setup parameter names
      // 	  int index=0;
      //     for (int i=0; i< NumParams_ee; i++) 
      // 	    for (int j=0; j < NumParams_eI; j++)
      // 	      for (int k=0; k<=j; k++) {
      // 		std::stringstream sstr;
      // 		sstr << id << "_" << i << "_" << j << "_" << k;
      // 		myVars.insert(sstr.str(),Parameters[index],true);
      // 		ParamArray(i,j,k) = ParamArray(i,k,j) = Parameters[index];
      // 		index++;
      // 	      }


      // 	  app_log() << "Parameter     Name      Value\n";
      // 	  myVars.print(app_log());
      // 	}
      // 	xmlCoefs = xmlCoefs->next;
      // }
      // reset();
      // print();
      // return true;
    }
    
    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
    }

    void resetParameters(const opt_variables_type& active)
    {
      int iparam = 0;
      for (int i=0; i<NumParams_ee; i++)
	for (int j=0; j<NumParams_eI; j++)
	  for (int k=0; k<=j; k++) {
	    int loc = myVars.where(iparam);
	    if (loc >=0) Parameters[iparam] = myVars[iparam] = active[loc];
	    ParamArray(i,j,k) = Parameters[iparam];
	    ParamArray(i,k,j) = Parameters[iparam];
	    iparam++;
	  }
      reset();
	    
      // for(int i=0; i<Parameters.size(); ++i) {
      //   int loc=myVars.where(i);
      //   if(loc>=0) Parameters[i]=myVars[i]=active[loc];
      // }
      
      if (ResetCount++ == 100) {
	ResetCount = 0;
	print();
      }
      reset();
    }


    void print()
    {
      const int N = 100;
      string fname = iSpecies + ".J3.h5";
      hid_t hid = H5Fcreate (fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      Array<real_type,3> val (N,N,N);
      for (int i=0; i<N; i++) {
	double r_12 = (real_type)i/(real_type)(N-1);
	for (int j=0; j<N; j++) {
	  double r_1I = (real_type)j/(real_type)(N-1) * 0.5*cutoff_radius;
	  for (int k=0; k<N; k++) {
	    double r_2I = (real_type)k/(real_type)(N-1) * 0.5*cutoff_radius;
	    val(i,j,k) = evaluate (r_12*(r_1I+r_2I), r_1I, r_2I);
	  }
	}
      }
      HDFAttribIO<Array<real_type,3> > coefs_attrib (SplineCoefs);
      HDFAttribIO<Array<real_type,3> > param_attrib (ParamArray);
      HDFAttribIO<Array<real_type,3> > val_attrib (val);      

      val_attrib.write (hid, "val");
      coefs_attrib.write (hid, "coefs");
      param_attrib.write (hid, "params");

      H5Fclose(hid);

      // string fname = (elementType != "") ? elementType : pairType;
      // fname = fname + ".dat";
      // //cerr << "Writing " << fname << " file.\n";
      // FILE *fout = fopen (fname.c_str(), "w");
      // for (double r=0.0; r<cutoff_radius; r+=0.001)
      // 	fprintf (fout, "%8.3f %16.10f\n", r, evaluate(r));
      // fclose(fout);
    }
    

    void print(std::ostream& os)
    {
      int n=100;
      real_type d=cutoff_radius/100.,r=0;
      real_type u,du,d2du;
      for(int i=0; i<n; ++i)
      {
        u=evaluate(r,du,d2du);
        os << setw(22) << r << setw(22) << u << setw(22) << du
          << setw(22) << d2du << std::endl;
        r+=d;
      }
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: PolynomialFunctor3D.h 1691 2007-02-01 21:51:50Z jnkim $ 
 ***************************************************************************/
