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
#ifndef QMCPLUSPLUS_POLYNOMIAL3D_FUNCTOR_H
#define QMCPLUSPLUS_POLYNOMIAL3D_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/LinearFit.h"
#include "Numerics/DeterminantOperators.h"
#include <cstdio>
#include <algorithm>

namespace qmcplusplus {

  struct PolynomialFunctor3D: public OptimizableFunctorBase {

    typedef real_type value_type;
    int N_eI, N_ee;
    Array<real_type,3> gamma;
    // Permutation vector, used when we need to pivot
    // columns
    vector<int> GammaPerm;

    Array<int,3> index;
    vector<bool> IndepVar;
    vector<real_type> GammaVec;
    int NumConstraints, NumGamma;
    Matrix<real_type> ConstraintMatrix;
    std::vector<real_type> Parameters;
    std::vector<std::string> ParameterNames;
    std::string iSpecies, eSpecies1, eSpecies2;
    int ResetCount;
    // Order of continuity
    const int C;

    ///constructor
    PolynomialFunctor3D(real_type ee_cusp=0.0, real_type eI_cusp=0.0) : 
      N_eI(0), N_ee(0), ResetCount(0), C(3)
    { 
      if (std::fabs(ee_cusp) > 0.0 || std::fabs(eI_cusp) > 0.0) {
	app_error() << "PolynomialFunctor3D does not support nonzero cusp.\n";
	abort();
      }
	
      cutoff_radius = 0.0;
    }

    OptimizableFunctorBase* makeClone() const 
    {
      return new PolynomialFunctor3D(*this);
    }

    void resize(int neI, int nee) 
    { 
      N_eI = neI;
      N_ee = nee;
      const double L = 0.5 * cutoff_radius;

      gamma.resize(N_eI+1, N_eI+1, N_ee+1);
      index.resize(N_eI+1, N_eI+1, N_ee+1);
      NumGamma = ((N_eI+1)*(N_eI+2)/2 * (N_ee+1));
      NumConstraints = (2*N_eI+1) + (N_eI+N_ee+1);
      int numParams = NumGamma - NumConstraints;
      Parameters.resize(numParams);
      GammaVec.resize(NumGamma);
      ConstraintMatrix.resize(NumConstraints, NumGamma);
      // Assign indices
      int num=0;
      for (int m=0; m<=N_eI; m++)
	for (int l=m; l<=N_eI; l++)
	  for (int n=0; n<=N_ee; n++) 
	    index(l,m,n) = index(m,l,n) = num++;
      assert (num == NumGamma);

      cerr << "NumGamma = " << NumGamma << endl;

      // Fill up contraint matrix
      // For 3 constraints and 2 parameters, we would have
      // |A00 A01 A02 A03 A04|  |g0|   |0 | 
      // |A11 A11 A12 A13 A14|  |g1|   |0 |
      // |A22 A21 A22 A23 A24|  |g2| = |0 |
      // | 0   0   0   1   0 |  |g3|   |p0|
      // | 0   0   0   0   1 |  |g4|   |p1|

      ConstraintMatrix = 0.0;
      // cerr << "ConstraintMatrix.size = " << ConstraintMatrix.size(0) 
      // 	   << " by " << ConstraintMatrix.size(1) << endl;
      // cerr << "index.size() = (" << index.size(0) << ", "
      // 	   << index.size(1) << ", " << index.size(2) << ").\n";
      int k;
      // e-e no-cusp constraint
      for (k=0; k<=2*N_eI; k++) {
	for (int m=0; m<=k; m++) {
	  int l = k - m;
	  int i = index(l,m,1);
	  if (l<=N_eI && m <=N_eI) {
	    if (l > m)
	      ConstraintMatrix(k,i) = 2.0;
	    else if (l == m)
	      ConstraintMatrix(k,i) = 1.0;
	  }
	}
      }
      // e-I no-cusp constraint
      for (int kp=0; kp<=N_eI+N_ee; kp++) {
      	if (kp <= N_ee) {
      	  ConstraintMatrix(k+kp,index(0,0,kp)) = (real_type) C;
      	  ConstraintMatrix(k+kp,index(0,1,kp)) = -L;
      	}
      	for (int l=1; l<=kp; l++) {
      	  int n = kp - l;
      	  if (n >= 0 && n <= N_ee && l <= N_eI) {
      	    ConstraintMatrix(k+kp,index(l,0,n)) = (real_type)C;
      	    ConstraintMatrix(k+kp,index(l,1,n)) = -L;
      	  }
      	}
      }

      fprintf (stderr, "Constraint matrix:\n");
      for (int i=0; i<NumConstraints; i++) {
	for (int j=0; j<NumGamma; j++)
	  fprintf (stderr, "%5.2f ", ConstraintMatrix(i,j));
	fprintf(stderr, "\n");
      }

      
      // Now, row-reduce constraint matrix
      GammaPerm.resize(NumGamma);
      IndepVar.resize(NumGamma, false);
      // Set identity permutation
      for (int i=0; i<NumGamma; i++)
	GammaPerm[i] = i;
      int col=-1; 
      for (int row=0; row<NumConstraints; row++) {
	int max_loc;
	real_type max_abs;
	do {
	  col++;
	  max_loc = row;
	  max_abs = std::fabs(ConstraintMatrix(row,col));
	  for (int ri=row+1; ri<NumConstraints; ri++) {
	    real_type abs_val = std::fabs(ConstraintMatrix(ri,col));
	    if (abs_val > max_abs) {
	      max_loc = ri;
	      max_abs = abs_val;
	    }
	  }
	  if (max_abs < 1.0e-6)
	    IndepVar[col] = true;
	} while (max_abs < 1.0e-6);
	
	ConstraintMatrix.swap_rows(row,max_loc);
	real_type lead_inv = 1.0/ConstraintMatrix(row,col);
	for (int c=0; c<NumGamma; c++)
	  ConstraintMatrix(row,c) *= lead_inv;
	// Now, eliminate column entries
	for (int ri=0; ri<NumConstraints; ri++) {
	  if (ri != row) {
	    real_type val = ConstraintMatrix(ri,col);
	    for (int c=0; c < NumGamma; c++)
	      ConstraintMatrix(ri,c) -= val * ConstraintMatrix(row,c);
	  }
	}
      }
      for (int c=col+1; c<NumGamma; c++)
	IndepVar[c] = true;

      fprintf (stderr, "Reduced Constraint matrix:\n");
      for (int i=0; i<NumConstraints; i++) {
	for (int j=0; j<NumGamma; j++)
	  fprintf (stderr, "%5.2f ", ConstraintMatrix(i,j));
	fprintf(stderr, "\n");
      }
      fprintf (stderr, "Independent vars = \n");
      for (int i=0; i<NumGamma; i++)
	if (IndepVar[i])
	  fprintf (stderr, "%d ", i);
      fprintf (stderr, "\n");

      // fprintf (stderr, "Inverse matrix:\n");
      // // Now, invert constraint matrix
      // Invert(ConstraintMatrix.data(), NumGamma, NumGamma);
      // for (int i=0; i<NumGamma; i++) {
      // 	for (int j=0; j<NumGamma; j++)
      // 	  fprintf (stderr, "%5.2f ", ConstraintMatrix(i,j));
      // 	fprintf(stderr, "\n");
      // }
    }
    
    void reset() 
    {
      const double L = 0.5 * cutoff_radius;
      std::fill(GammaVec.begin(), GammaVec.end(), 0.0);

      // First, set all independent variables
      int var=0;
      for (int i=0; i<NumGamma; i++)
	if (IndepVar[i])
	  GammaVec[i] = Parameters[var++];
      
      // Now, set dependent variables
      var = 0;
      //      cerr << "NumConstraints = " << NumConstraints << endl;
      for (int i=0; i<NumGamma; i++)
	if (!IndepVar[i]) {
	  // fprintf (stderr, "constraintMatrix(%d,%d) = %1.10f\n",
	  // 	   var, i, ConstraintMatrix(var,i));
	  assert (std::fabs(ConstraintMatrix(var,i) -1.0) < 1.0e-10);
	  for (int j=0; j<NumGamma; j++)
	    if (i != j)
	      GammaVec[i] -= ConstraintMatrix(var,j) * GammaVec[j];
	  var++;
	}

      int num=0;
      for (int m=0; m<=N_eI; m++)
	for (int l=m; l<=N_eI; l++)
	  for (int n=0; n<=N_ee; n++) 
	    //	    gamma(m,l,n) = gamma(l,m,n) = unpermuted[num++];
	    gamma(m,l,n) = gamma(l,m,n) = GammaVec[num++];

      // Now check that constraints have been satisfied
      // e-e constraints
      for (int k=0; k<=2*N_eI; k++) {
	real_type sum = 0.0;
	for (int m=0; m<=k; m++) {
	  int l = k - m;
	  int i = index(l,m,1);
	  if (l<=N_eI && m <=N_eI) {
	    if (l > m) 
	      sum += 2.0*GammaVec[i];
	    else if (l == m)
	      sum += GammaVec[i];
	  }
	}
	if (std::fabs(sum) > 1.0e-9) 
	  cerr << "error in k = " << k << "  sum = " << sum << endl;
      }


      for (int k=0; k<=2*N_eI; k++) {
	real_type sum=0.0;
	for (int l=0; l<=k; l++) {
	  int m = k - l;
	  if (m <= N_eI && l <= N_eI) {
	    // fprintf (stderr, "k = %d gamma(%d, %d, 1) = %1.8f\n", k, l, m,
	    // 	     gamma(l,m,1));
	    sum += gamma(l,m,1);
	  }
	}
	if (std::fabs(sum) > 1.0e-10) {
	  app_error() << "e-e constraint not satisfied in PolynomialFunctor3D:  k=" 
		      << k << "  sum=" << sum << endl;
	  abort();
	}
      }

      // e-I constraints
      for (int k=0; k<=N_eI+N_ee; k++) {
	real_type sum = 0.0;
	for (int m=0; m<=k; m++) {
	  int n = k - m;
	  if (m <= N_eI && n <= N_ee) {
	    sum += (real_type)C*gamma(0,m,n) - L*gamma(1,m,n);
	    // fprintf (stderr, 
	    // 	     "k = %d gamma(0,%d,%d) = %1.8f  gamma(1,%d,%d)=%1.8f\n",
	    // 	     k, m, n, gamma(0,m,n), m, n, gamma(1,m,n));
	  }
	}
	if (std::fabs(sum) > 1.0e-10) {
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
      for (int l=0; l<=N_eI; l++) {
	real_type r2m=1.0; 
	for (int m=0; m<=N_eI; m++) {
	  real_type r2n=1.0;
	  for (int n=0; n<=N_ee; n++) {
	    val += gamma(l,m,n)*r2l*r2m*r2n;
	    r2n *= r_12;
	  }
	  r2m *= r_2I;
	}
	r2l *= r_1I;
      }
      // for (int i=0; i<C; i++)
      // 	val *= (r_1I - 0.5*cutoff_radius)*(r_2I - 0.5*cutoff_radius);
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

      // r2l[0] = r2m[0] = 1.0;
      // for (int i=1; i<N_eI; i++) {
      // 	r2l[i] = r2l[i-1] * r_1I;
      // 	r2m[i] = r2l[i-1] * r_2I;
      // }
      // r2n[0] = 1.0;
      // for (int i=1; i<N_ee; i++) 
      // 	r2n[i] = r2n[i-1] * r_12;

      

      real_type r2l(1.0), r2l_1(0.0), r2l_2(0.0), lf(0.0);
      for (int l=0; l<N_eI; l++) {
	real_type r2m(1.0), r2m_1(0.0), r2m_2(0.0), mf(0.0);
	for (int m=0; m<N_eI; m++) {
	  real_type r2n(1.0), r2n_1(0.0), r2n_2(0.0), nf(0.0);
	  for (int n=0; n<N_ee; n++) {
	    real_type g = gamma(l,m,n);
	    val += g*r2l*r2m*r2n;

	    grad[0] += nf * g *r2l   * r2m   * r2n_1;
	    grad[1] += lf * g *r2l_1 * r2m   * r2n  ;
	    grad[2] += mf * g *r2l   * r2m_1 * r2n  ;

	    hess(0,0) += nf*(nf-1.0) * g * r2l   * r2m   * r2n_2  ;
	    hess(0,1) += nf*lf       * g * r2l_1 * r2m   * r2n_1  ;
	    hess(0,2) += nf*mf       * g * r2l   * r2m_1 * r2n_1  ;
	    hess(1,1) += lf*(lf-1.0) * g * r2l_2 * r2m   * r2n    ;
	    hess(1,2) += lf*mf       * g * r2l_1 * r2m_1 * r2n    ;
	    hess(2,2) += mf*(mf-1.0) * g * r2l   * r2m_2 * r2n    ;

	    r2n_2 = r2n_1;
	    r2n_1 = r2n;
	    r2n *= r_12;
	    nf += 1.0;
	  }
	  r2m_2 = r2m_1;
	  r2m_1 = r2m;
	  r2m *= r_2I;
	  mf += 1.0;
	}
	r2l_2 = r2l_1;
	r2l_1 = r2l;
	r2l *= r_1I;
	lf += 1.0;
      }
      hess(1,0) = hess(0,1);
      hess(2,0) = hess(0,2);
      hess(2,1) = hess(1,2);

      // for (int i=0; i<C; i++)
      // 	val *= (r_1I - 0.5*cutoff_radius)*(r_2I - 0.5*cutoff_radius);

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
      ReportEngine PRE("PolynomialFunctor3D","put(xmlNodePtr)");
      
      // //CuspValue = -1.0e10;
      // NumParams_eI = NumParams_ee = 0;
      cutoff_radius = 0.0;
      OhmmsAttributeSet rAttrib;
      rAttrib.add(N_ee,   "esize");
      rAttrib.add(N_eI,   "isize");
      rAttrib.add(cutoff_radius,  "rcut");
      rAttrib.put(cur);

      if (N_eI == 0) 
        PRE.error("You must specify a positive number for \"isize\"",true);
      if (N_ee == 0) 
        PRE.error("You must specify a positive number for \"esize\"",true);

      // app_log() << " esize = " << NumParams_ee << " parameters " << endl;
      // app_log() << " isize = " << NumParams_eI << " parameters " << endl;
      // app_log() << " rcut = " << cutoff_radius << endl;

      resize (N_eI, N_ee);

      // Now read coefficents
      xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
      while (xmlCoefs != NULL) 
      {
        string cname((const char*)xmlCoefs->name);
        if (cname == "coefficients")  {
          string type("0"), id("0");
          OhmmsAttributeSet cAttrib;
          cAttrib.add(id, "id");
          cAttrib.add(type, "type");
          cAttrib.put(xmlCoefs);

          if (type != "Array") {
            PRE.error( "Unknown correlation type " + type + 
      		       " in PolynomialFunctor3D." + "Resetting to \"Array\"");
            xmlNewProp (xmlCoefs, (const xmlChar*) "type", 
			(const xmlChar*) "Array");
          }

      	  vector<real_type> params;
      	  putContent(params, xmlCoefs);
      	  if (params.size() == Parameters.size()) 
      	    Parameters = params;
      	  else {
      	    app_error() << "Expected " << Parameters.size() << " parameters,"
      			<< " but found " << params.size()
      			<< " in PolynomialFunctor3D.\n";
      	    abort();
      	  }
	  
	  
          // Setup parameter names
      	  int index=0;
	  for (int i=0; i<Parameters.size(); i++) {
	    std::stringstream sstr;
	    sstr << id << "_" << i;;
	    myVars.insert(sstr.str(),Parameters[i],true);
	  }
          // for (int i=0; i< N_ee; i++) 
      	  //   for (int j=0; j < N_eI; j++)
      	  //     for (int k=0; k<=j; k++) {
      	  // 	std::stringstream sstr;
      	  // 	sstr << id << "_" << i << "_" << j << "_" << k;
      	  // 	myVars.insert(sstr.str(),Parameters[index],true);
      	  // 	ParamArray(i,j,k) = ParamArray(i,k,j) = Parameters[index];
      	  // 	index++;
      	  //     }

      	  app_log() << "Parameter     Name      Value\n";
      	  myVars.print(app_log());
      	}
      	xmlCoefs = xmlCoefs->next;
      }
      reset();
      print();
      return true;
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
      for (int i=0; i<Parameters.size(); i++) {
	int loc = myVars.where(i);
	if (loc >= 0)
	  Parameters[i] = myVars[i] = active[loc];
      }
	
      // int iparam = 0;
      // for (int i=0; i<N_ee; i++)
      // 	for (int j=0; j<N_eI; j++)
      // 	  for (int k=0; k<=j; k++) {
      // 	    int loc = myVars.where(iparam);
      // 	    if (loc >=0) Parameters[iparam] = myVars[iparam] = active[loc];
      // 	    ParamArray(i,j,k) = Parameters[iparam];
      // 	    ParamArray(i,k,j) = Parameters[iparam];
      // 	    iparam++;
      // 	  }
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
      // HDFAttribIO<Array<real_type,3> > coefs_attrib (SplineCoefs);
      // HDFAttribIO<Array<real_type,3> > param_attrib (ParamArray);
      HDFAttribIO<Array<real_type,3> > val_attrib (val);      

      val_attrib.write (hid, "val");
      // coefs_attrib.write (hid, "coefs");
      // param_attrib.write (hid, "params");

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
