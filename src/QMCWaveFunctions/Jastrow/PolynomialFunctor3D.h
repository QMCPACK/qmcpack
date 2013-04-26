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

namespace qmcplusplus
{

struct PolynomialFunctor3D: public OptimizableFunctorBase
{

  typedef real_type value_type;
  int N_eI, N_ee;
  Array<real_type,3> gamma;
  Array<real_type,3> dval_dgamma;
  Array<TinyVector<real_type,3>,3> dgrad_dgamma;
  Array<Tensor<real_type,3>,3> dhess_dgamma;
  // Permutation vector, used when we need to pivot
  // columns
  vector<int> GammaPerm;

  Array<int,3> index;
  vector<bool> IndepVar;
  vector<real_type> GammaVec, dval_Vec;
  vector<TinyVector<real_type,3> > dgrad_Vec;
  vector<Tensor<real_type,3> > dhess_Vec;
  int NumConstraints, NumGamma;
  Matrix<real_type> ConstraintMatrix;
  std::vector<real_type> Parameters, d_valsFD;
  std::vector<TinyVector<real_type,3> > d_gradsFD;
  std::vector<Tensor<real_type,3> > d_hessFD;
  std::vector<std::string> ParameterNames;
  std::string iSpecies, eSpecies1, eSpecies2;
  int ResetCount;
  real_type scale;
  // Order of continuity
  const int C;
  bool notOpt;

  ///constructor
  PolynomialFunctor3D(real_type ee_cusp=0.0, real_type eI_cusp=0.0) :
    N_eI(0), N_ee(0), ResetCount(0), C(3), scale(1.0), notOpt(false)
  {
    if (std::fabs(ee_cusp) > 0.0 || std::fabs(eI_cusp) > 0.0)
    {
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
    dval_dgamma.resize(N_eI+1, N_eI+1, N_ee+1);
    dgrad_dgamma.resize(N_eI+1, N_eI+1, N_ee+1);
    dhess_dgamma.resize(N_eI+1, N_eI+1, N_ee+1);
    index.resize(N_eI+1, N_eI+1, N_ee+1);
    NumGamma = ((N_eI+1)*(N_eI+2)/2 * (N_ee+1));
    NumConstraints = (2*N_eI+1) + (N_eI+N_ee+1);
    int numParams = NumGamma - NumConstraints;
    Parameters.resize(numParams);
    d_valsFD.resize(numParams);
    d_gradsFD.resize(numParams);
    d_hessFD.resize(numParams);
    GammaVec.resize(NumGamma);
    dval_Vec.resize(NumGamma);
    dgrad_Vec.resize(NumGamma);
    dhess_Vec.resize(NumGamma);
    ConstraintMatrix.resize(NumConstraints, NumGamma);
    // Assign indices
    int num=0;
    for (int m=0; m<=N_eI; m++)
      for (int l=m; l<=N_eI; l++)
        for (int n=0; n<=N_ee; n++)
          index(l,m,n) = index(m,l,n) = num++;
    assert (num == NumGamma);
//       cerr << "NumGamma = " << NumGamma << endl;
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
    for (k=0; k<=2*N_eI; k++)
    {
      for (int m=0; m<=k; m++)
      {
        int l = k - m;
        if (l<=N_eI && m <=N_eI)
        {
          int i = index(l,m,1);
          if (l > m)
            ConstraintMatrix(k,i) = 2.0;
          else
            if (l == m)
              ConstraintMatrix(k,i) = 1.0;
        }
      }
    }
    // e-I no-cusp constraint
    for (int kp=0; kp<=N_eI+N_ee; kp++)
    {
      if (kp <= N_ee)
      {
        ConstraintMatrix(k+kp,index(0,0,kp)) = (real_type) C;
        ConstraintMatrix(k+kp,index(0,1,kp)) = -L;
      }
      for (int l=1; l<=kp; l++)
      {
        int n = kp - l;
        if (n >= 0 && n <= N_ee && l <= N_eI)
        {
          ConstraintMatrix(k+kp,index(l,0,n)) = (real_type)C;
          ConstraintMatrix(k+kp,index(l,1,n)) = -L;
        }
      }
    }
//    {
//       fprintf (stderr, "Constraint matrix:\n");
//       for (int i=0; i<NumConstraints; i++) {
// 	for (int j=0; j<NumGamma; j++)
// 	  fprintf (stderr, "%5.2f ", ConstraintMatrix(i,j));
// 	fprintf(stderr, "\n");
//       }
//     }
    // Now, row-reduce constraint matrix
    GammaPerm.resize(NumGamma);
    IndepVar.resize(NumGamma, false);
    // Set identity permutation
    for (int i=0; i<NumGamma; i++)
      GammaPerm[i] = i;
    int col=-1;
    for (int row=0; row<NumConstraints; row++)
    {
      int max_loc;
      real_type max_abs;
      do
      {
        col++;
        max_loc = row;
        max_abs = std::fabs(ConstraintMatrix(row,col));
        for (int ri=row+1; ri<NumConstraints; ri++)
        {
          real_type abs_val = std::fabs(ConstraintMatrix(ri,col));
          if (abs_val > max_abs)
          {
            max_loc = ri;
            max_abs = abs_val;
          }
        }
        if (max_abs < 1.0e-6)
          IndepVar[col] = true;
      }
      while (max_abs < 1.0e-6);
      ConstraintMatrix.swap_rows(row,max_loc);
      real_type lead_inv = 1.0/ConstraintMatrix(row,col);
      for (int c=0; c<NumGamma; c++)
        ConstraintMatrix(row,c) *= lead_inv;
      // Now, eliminate column entries
      for (int ri=0; ri<NumConstraints; ri++)
      {
        if (ri != row)
        {
          real_type val = ConstraintMatrix(ri,col);
          for (int c=0; c < NumGamma; c++)
            ConstraintMatrix(ri,c) -= val * ConstraintMatrix(row,c);
        }
      }
    }
    for (int c=col+1; c<NumGamma; c++)
      IndepVar[c] = true;
//       fprintf (stderr, "Reduced Constraint matrix:\n");
//       for (int i=0; i<NumConstraints; i++) {
// 	for (int j=0; j<NumGamma; j++)
// 	  fprintf (stderr, "%5.2f ", ConstraintMatrix(i,j));
// 	fprintf(stderr, "\n");
//       }
//       fprintf (stderr, "Independent vars = \n");
//       for (int i=0; i<NumGamma; i++)
// 	if (IndepVar[i])
// 	  fprintf (stderr, "%d ", i);
//       fprintf (stderr, "\n");
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
    resize(N_eI, N_ee);
    reset_gamma();
  }

  void reset_gamma()
  {
    // fprintf (stderr, "Paramters:\n");
    // for (int i=0; i<Parameters.size(); i++)
    // 	fprintf (stderr, " %16.10e\n", Parameters[i]);
    const double L = 0.5 * cutoff_radius;
    std::fill(GammaVec.begin(), GammaVec.end(), 0.0);
    // First, set all independent variables
    int var=0;
    for (int i=0; i<NumGamma; i++)
      if (IndepVar[i])
        GammaVec[i] = scale*Parameters[var++];
    assert (var == Parameters.size());
    // Now, set dependent variables
    var = 0;
    //      cerr << "NumConstraints = " << NumConstraints << endl;
    for (int i=0; i<NumGamma; i++)
      if (!IndepVar[i])
      {
        // fprintf (stderr, "constraintMatrix(%d,%d) = %1.10f\n",
        // 	   var, i, ConstraintMatrix(var,i));
        assert (std::fabs(ConstraintMatrix(var,i) -1.0) < 1.0e-6);
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
    for (int k=0; k<=2*N_eI; k++)
    {
      real_type sum = 0.0;
      for (int m=0; m<=k; m++)
      {
        int l = k - m;
        if (l<=N_eI && m <=N_eI)
        {
          int i = index(l,m,1);
          if (l > m)
            sum += 2.0*GammaVec[i];
          else
            if (l == m)
              sum += GammaVec[i];
        }
      }
      if (std::fabs(sum) > 1.0e-9)
        cerr << "error in k = " << k << "  sum = " << sum << endl;
    }
    for (int k=0; k<=2*N_eI; k++)
    {
      real_type sum=0.0;
      for (int l=0; l<=k; l++)
      {
        int m = k - l;
        if (m <= N_eI && l <= N_eI)
        {
          // fprintf (stderr, "k = %d gamma(%d, %d, 1) = %1.8f\n", k, l, m,
          // 	     gamma(l,m,1));
          sum += gamma(l,m,1);
        }
      }
      if (std::fabs(sum) > 1.0e-6)
      {
        app_error() << "e-e constraint not satisfied in PolynomialFunctor3D:  k="
                    << k << "  sum=" << sum << endl;
        abort();
      }
    }
    // e-I constraints
    for (int k=0; k<=N_eI+N_ee; k++)
    {
      real_type sum = 0.0;
      for (int m=0; m<=k; m++)
      {
        int n = k - m;
        if (m <= N_eI && n <= N_ee)
        {
          sum += (real_type)C*gamma(0,m,n) - L*gamma(1,m,n);
          // fprintf (stderr,
          // 	     "k = %d gamma(0,%d,%d) = %1.8f  gamma(1,%d,%d)=%1.8f\n",
          // 	     k, m, n, gamma(0,m,n), m, n, gamma(1,m,n));
        }
      }
      if (std::fabs(sum) > 1.0e-6)
      {
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
    const real_type L = 0.5*cutoff_radius;
    if (r_1I >= L || r_2I >= L)
      return 0.0;
    real_type val = 0.0;
    real_type r2l=1.0;
    for (int l=0; l<=N_eI; l++)
    {
      real_type r2m=1.0;
      for (int m=0; m<=N_eI; m++)
      {
        real_type r2n=1.0;
        for (int n=0; n<=N_ee; n++)
        {
          val += gamma(l,m,n)*r2l*r2m*r2n;
          r2n *= r_12;
        }
        r2m *= r_2I;
      }
      r2l *= r_1I;
    }
    for (int i=0; i<C; i++)
      val *= (r_1I - L)*(r_2I - L);
    return val;
  }


  inline real_type evaluate(real_type r_12, real_type r_1I, real_type r_2I,
                            TinyVector<real_type,3> &grad,
                            Tensor<real_type,3> &hess)
  {
    const real_type L = 0.5*cutoff_radius;
    if (r_1I >= L || r_2I >= L)
    {
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
    grad = 0.0;
    hess = 0.0;
    real_type r2l(1.0), r2l_1(0.0), r2l_2(0.0), lf(0.0);
    for (int l=0; l<=N_eI; l++)
    {
      real_type r2m(1.0), r2m_1(0.0), r2m_2(0.0), mf(0.0);
      for (int m=0; m<=N_eI; m++)
      {
        real_type r2n(1.0), r2n_1(0.0), r2n_2(0.0), nf(0.0);
        for (int n=0; n<=N_ee; n++)
        {
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
    for (int i=0; i<C; i++)
    {
      hess(0,0)=(r_1I - L)*(r_2I - L)*hess(0,0);
      hess(0,1)=(r_1I - L)*(r_2I - L)*hess(0,1)+ (r_2I - L)*grad[0];
      hess(0,2)=(r_1I - L)*(r_2I - L)*hess(0,2)+ (r_1I - L)*grad[0];
      hess(1,1)=(r_1I - L)*(r_2I - L)*hess(1,1)+ 2.0*(r_2I - L)*grad[1];
      hess(1,2)=(r_1I - L)*(r_2I - L)*hess(1,2)+ (r_1I - L)*grad[1] + (r_2I - L)*grad[2] +  val;
      hess(2,2)=(r_1I - L)*(r_2I - L)*hess(2,2)+ 2.0*(r_1I - L)*grad[2];
      grad[0] = (r_1I - L)*(r_2I - L)*grad[0];
      grad[1] = (r_1I - L)*(r_2I - L)*grad[1] + (r_2I - L) * val;
      grad[2] = (r_1I - L)*(r_2I - L)*grad[2] + (r_1I - L) * val;
      val *= (r_1I - L)*(r_2I - L);
    }
    hess(1,0) = hess(0,1);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);
    return val;
  }


  inline real_type evaluate(real_type r_12, real_type r_1I, real_type r_2I,
                            TinyVector<real_type,3> &grad,
                            Tensor<real_type,3> &hess,
                            TinyVector<Tensor<real_type,3>,3> &d3)
  {
    grad = 0.0;
    hess = 0.0;
    d3 = grad;
    const real_type L = 0.5*cutoff_radius;
    if (r_1I >= L || r_2I >= L)
      return 0.0;
    real_type val = 0.0;
    real_type r2l(1.0), r2l_1(0.0), r2l_2(0.0), r2l_3, lf(0.0);
    for (int l=0; l<=N_eI; l++)
    {
      real_type r2m(1.0), r2m_1(0.0), r2m_2(0.0), r2m_3, mf(0.0);
      for (int m=0; m<=N_eI; m++)
      {
        real_type r2n(1.0), r2n_1(0.0), r2n_2(0.0), r2n_3, nf(0.0);
        for (int n=0; n<=N_ee; n++)
        {
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
          d3[0](0,0) += nf*(nf-1.0)*(nf-2.0) * g * r2l   * r2m   * r2n_3 ;
          d3[0](0,1) += nf*(nf-1.0)*lf       * g * r2l_1 * r2m   * r2n_2 ;
          d3[0](0,2) += nf*(nf-1.0)*mf       * g * r2l   * r2m_1 * r2n_2 ;
          d3[0](1,1) += nf*lf*(lf-1.0)       * g * r2l_2 * r2m   * r2n_1 ;
          d3[0](1,2) += nf*lf*mf             * g * r2l_1 * r2m_1 * r2n_1 ;
          d3[0](2,2) += nf*mf*(mf-1.0)       * g * r2l   * r2m_2 * r2n_1 ;
          d3[1](1,1) += lf*(lf-1.0)*(lf-2.0) * g * r2l_3 * r2m   * r2n   ;
          d3[1](1,2) += lf*(lf-1.0)*mf       * g * r2l_2 * r2m_1 * r2n   ;
          d3[1](2,2) += lf*mf*(mf-1.0)       * g * r2l_1 * r2m_2 * r2n   ;
          d3[2](2,2) += mf*(mf-1.0)*(mf-2.0) * g * r2l   * r2m_3 * r2n   ;
          r2n_3 = r2n_2;
          r2n_2 = r2n_1;
          r2n_1 = r2n;
          r2n *= r_12;
          nf += 1.0;
        }
        r2m_3 = r2m_2;
        r2m_2 = r2m_1;
        r2m_1 = r2m;
        r2m *= r_2I;
        mf += 1.0;
      }
      r2l_3 = r2l_2;
      r2l_2 = r2l_1;
      r2l_1 = r2l;
      r2l *= r_1I;
      lf += 1.0;
    }
    for (int i=0; i<C; i++)
    {
      d3[0](0,0) = (r_1I - L)*(r_2I - L)*d3[0](0,0);
      d3[0](0,1) = (r_1I - L)*(r_2I - L)*d3[0](0,1) + (r_2I - L)*hess(0,0);
      d3[0](0,2) = (r_1I - L)*(r_2I - L)*d3[0](0,2) + (r_1I - L)*hess(0,0);
      d3[0](1,1) = (r_1I - L)*(r_2I - L)*d3[0](1,1) + 2.0*(r_2I - L)*hess(0,1);
      d3[0](1,2) = (r_1I - L)*(r_2I - L)*d3[0](1,2) + (r_1I - L)*hess(0,1) + (r_2I - L)*hess(0,2) + grad[0];
      d3[0](2,2) = (r_1I - L)*(r_2I - L)*d3[0](2,2) + 2.0*(r_1I-L)*hess(0,2);
      d3[1](1,1) = (r_1I - L)*(r_2I - L)*d3[1](1,1) + 3.0*(r_2I-L)*hess(1,1);
      d3[1](1,2) = (r_1I - L)*(r_2I - L)*d3[1](1,2) +
                   2.0*(r_2I - L)*hess(1,2) + 2.0*grad[1] + (r_1I - L)*hess(1,1);
      d3[1](2,2) = (r_1I - L)*(r_2I - L)*d3[1](2,2) +
                   2.0*(r_1I - L)*hess(1,2) + 2.0*grad[2] + (r_2I - L)*hess(2,2);
      d3[2](2,2) = (r_1I - L)*(r_2I - L)*d3[2](2,2) + 3.0*(r_1I - L)*hess(2,2);
      hess(0,0)=(r_1I - L)*(r_2I - L)*hess(0,0);
      hess(0,1)=(r_1I - L)*(r_2I - L)*hess(0,1)+ (r_2I - L)*grad[0];
      hess(0,2)=(r_1I - L)*(r_2I - L)*hess(0,2)+ (r_1I - L)*grad[0];
      hess(1,1)=(r_1I - L)*(r_2I - L)*hess(1,1)+ 2.0*(r_2I - L)*grad[1];
      hess(1,2)=(r_1I - L)*(r_2I - L)*hess(1,2)+ (r_1I - L)*grad[1] + (r_2I - L)*grad[2] +  val;
      hess(2,2)=(r_1I - L)*(r_2I - L)*hess(2,2)+ 2.0*(r_1I - L)*grad[2];
      grad[0] = (r_1I - L)*(r_2I - L)*grad[0];
      grad[1] = (r_1I - L)*(r_2I - L)*grad[1] + (r_2I - L) * val;
      grad[2] = (r_1I - L)*(r_2I - L)*grad[2] + (r_1I - L) * val;
      val *= (r_1I - L)*(r_2I - L);
    }
    hess(1,0) = hess(0,1);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);
    d3[0](1,0) = d3[0](0,1);
    d3[0](2,0) = d3[0](0,2);
    d3[0](2,1) = d3[0](1,2);
    d3[1](0,0) = d3[0](1,1);
    d3[1](0,1) = d3[0](0,1);
    d3[1](0,2) = d3[0](1,2);
    d3[1](1,0) = d3[0](0,1);
    d3[1](2,0) = d3[0](1,2);
    d3[1](2,1) = d3[1](1,2);
    d3[2](0,0) = d3[0](0,2);
    d3[2](0,1) = d3[0](1,2);
    d3[2](0,2) = d3[0](2,2);
    d3[2](1,0) = d3[0](1,2);
    d3[2](1,1) = d3[1](1,2);
    d3[2](1,2) = d3[1](2,2);
    d3[2](2,0) = d3[0](2,2);
    d3[2](2,1) = d3[1](2,2);
    return val;
  }


  inline real_type evaluate(real_type r, real_type rinv)
  {
    return 0.0;
  }


  inline bool
  evaluateDerivativesFD (real_type r_12, real_type r_1I, real_type r_2I,
                         vector<double> &d_vals,
                         vector<TinyVector<real_type,3> >& d_grads,
                         vector<Tensor<real_type,3> > &d_hess)
  {
    const real_type eps = 1.0e-6;
    assert (d_vals.size() == Parameters.size());
    assert (d_grads.size() == Parameters.size());
    assert (d_hess.size() == Parameters.size());
    for (int ip=0; ip < Parameters.size(); ip++)
    {
      real_type v_plus, v_minus;
      TinyVector<real_type,3> g_plus, g_minus;
      Tensor<real_type,3> h_plus, h_minus;
      real_type save_p = Parameters[ip];
      Parameters[ip]  = save_p + eps;
      reset_gamma();
      v_plus = evaluate (r_12, r_1I, r_2I, g_plus, h_plus);
      Parameters[ip]  = save_p - eps;
      reset_gamma();
      v_minus = evaluate (r_12, r_1I, r_2I, g_minus, h_minus);
      Parameters[ip]  = save_p;
      reset_gamma();
      real_type dp_inv = 0.5/eps;
      d_vals[ip]  = dp_inv * (v_plus - v_minus);
      d_grads[ip] = dp_inv * (g_plus - g_minus);
      d_hess[ip]  = dp_inv * (h_plus - h_minus);
    }
    return true;
  }


  inline bool
  evaluateDerivatives (real_type r_12, real_type r_1I, real_type r_2I,
                       vector<double> &d_vals,
                       vector<TinyVector<real_type,3> >& d_grads,
                       vector<Tensor<real_type,3> > &d_hess)
  {
    const real_type L = 0.5*cutoff_radius;
    if (r_1I >= L || r_2I >= L)
      return false;
    real_type val = 0.0;
    TinyVector<real_type,3> grad;
    Tensor<real_type,3> hess;
    real_type r2l(1.0), r2l_1(0.0), r2l_2(0.0), r2l_3, lf(0.0);
    for (int l=0; l<=N_eI; l++)
    {
      real_type r2m(1.0), r2m_1(0.0), r2m_2(0.0), r2m_3, mf(0.0);
      for (int m=0; m<=N_eI; m++)
      {
        real_type r2n(1.0), r2n_1(0.0), r2n_2(0.0), r2n_3, nf(0.0);
        for (int n=0; n<=N_ee; n++)
        {
          real_type g = gamma(l,m,n);
          val += g * r2l*r2m*r2n;
          grad[0] += g * nf *r2l   * r2m   * r2n_1;
          grad[1] += g * lf *r2l_1 * r2m   * r2n  ;
          grad[2] += g * mf *r2l   * r2m_1 * r2n  ;
          hess(0,0) += g * nf*(nf-1.0) * r2l   * r2m   * r2n_2  ;
          hess(0,1) += g * nf*lf       * r2l_1 * r2m   * r2n_1  ;
          hess(0,2) += g * nf*mf       * r2l   * r2m_1 * r2n_1  ;
          hess(1,1) += g * lf*(lf-1.0) * r2l_2 * r2m   * r2n    ;
          hess(1,2) += g * lf*mf       * r2l_1 * r2m_1 * r2n    ;
          hess(2,2) += g * mf*(mf-1.0) * r2l   * r2m_2 * r2n    ;
          dval_dgamma(l,m,n) = r2l*r2m*r2n;
          dgrad_dgamma(l,m,n)[0] = nf *r2l   * r2m   * r2n_1;
          dgrad_dgamma(l,m,n)[1] = lf *r2l_1 * r2m   * r2n  ;
          dgrad_dgamma(l,m,n)[2] = mf *r2l   * r2m_1 * r2n  ;
          dhess_dgamma(l,m,n)(0,0) = nf*(nf-1.0) * r2l   * r2m   * r2n_2  ;
          dhess_dgamma(l,m,n)(0,1) = nf*lf       * r2l_1 * r2m   * r2n_1  ;
          dhess_dgamma(l,m,n)(0,2) = nf*mf       * r2l   * r2m_1 * r2n_1  ;
          dhess_dgamma(l,m,n)(1,1) = lf*(lf-1.0) * r2l_2 * r2m   * r2n    ;
          dhess_dgamma(l,m,n)(1,2) = lf*mf       * r2l_1 * r2m_1 * r2n    ;
          dhess_dgamma(l,m,n)(2,2) = mf*(mf-1.0) * r2l   * r2m_2 * r2n    ;
          r2n_3 = r2n_2;
          r2n_2 = r2n_1;
          r2n_1 = r2n;
          r2n *= r_12;
          nf += 1.0;
        }
        r2m_3 = r2m_2;
        r2m_2 = r2m_1;
        r2m_1 = r2m;
        r2m *= r_2I;
        mf += 1.0;
      }
      r2l_3 = r2l_2;
      r2l_2 = r2l_1;
      r2l_1 = r2l;
      r2l *= r_1I;
      lf += 1.0;
    }
    for (int i=0; i<C; i++)
    {
      hess(0,0)=(r_1I - L)*(r_2I - L)*hess(0,0);
      hess(0,1)=(r_1I - L)*(r_2I - L)*hess(0,1)+ (r_2I - L)*grad[0];
      hess(0,2)=(r_1I - L)*(r_2I - L)*hess(0,2)+ (r_1I - L)*grad[0];
      hess(1,1)=(r_1I - L)*(r_2I - L)*hess(1,1)+ 2.0*(r_2I - L)*grad[1];
      hess(1,2)=(r_1I - L)*(r_2I - L)*hess(1,2)+ (r_1I - L)*grad[1] + (r_2I - L)*grad[2] +  val;
      hess(2,2)=(r_1I - L)*(r_2I - L)*hess(2,2)+ 2.0*(r_1I - L)*grad[2];
      grad[0] = (r_1I - L)*(r_2I - L)*grad[0];
      grad[1] = (r_1I - L)*(r_2I - L)*grad[1] + (r_2I - L) * val;
      grad[2] = (r_1I - L)*(r_2I - L)*grad[2] + (r_1I - L) * val;
      val *= (r_1I - L)*(r_2I - L);
      for (int l=0; l<=N_eI; l++)
        for (int m=0; m<=N_eI; m++)
          for (int n=0; n<=N_ee; n++)
          {
            dhess_dgamma(l,m,n)(0,0)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(0,0);
            dhess_dgamma(l,m,n)(0,1)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(0,1)+ (r_2I - L)*dgrad_dgamma(l,m,n)[0];
            dhess_dgamma(l,m,n)(0,2)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(0,2)+ (r_1I - L)*dgrad_dgamma(l,m,n)[0];
            dhess_dgamma(l,m,n)(1,1)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(1,1)+ 2.0*(r_2I - L)*dgrad_dgamma(l,m,n)[1];
            dhess_dgamma(l,m,n)(1,2)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(1,2)+ (r_1I - L)*dgrad_dgamma(l,m,n)[1]
                                     + (r_2I - L)*dgrad_dgamma(l,m,n)[2] +  dval_dgamma(l,m,n);
            dhess_dgamma(l,m,n)(2,2)=(r_1I - L)*(r_2I - L)*dhess_dgamma(l,m,n)(2,2)+ 2.0*(r_1I - L)*dgrad_dgamma(l,m,n)[2];
            dgrad_dgamma(l,m,n)[0] = (r_1I - L)*(r_2I - L)*dgrad_dgamma(l,m,n)[0];
            dgrad_dgamma(l,m,n)[1] = (r_1I - L)*(r_2I - L)*dgrad_dgamma(l,m,n)[1] + (r_2I - L) * dval_dgamma(l,m,n);
            dgrad_dgamma(l,m,n)[2] = (r_1I - L)*(r_2I - L)*dgrad_dgamma(l,m,n)[2] + (r_1I - L) * dval_dgamma(l,m,n);
            dval_dgamma(l,m,n) *= (r_1I - L)*(r_2I - L);
          }
    }
    // Now, pack into vectors
    int num=0;
    for (int m=0; m<=N_eI; m++)
      for (int l=m; l<=N_eI; l++)
        for (int n=0; n<=N_ee; n++)
        {
          // Don't double-count diagonal terms
          real_type nodiag = (l==m) ? 0.0 : 1.0;
          dval_Vec[num] = scale*(dval_dgamma(l,m,n) + nodiag*dval_dgamma(m,l,n));
          for (int i=0; i<3; i++)
          {
            dgrad_Vec[num][i] = scale*(dgrad_dgamma(l,m,n)[i] + nodiag*dgrad_dgamma(m,l,n)[i]);
            for (int j=i; j<3; j++)
            {
              dhess_Vec[num](i,j) = scale*(dhess_dgamma(l,m,n)(i,j) + nodiag*dhess_dgamma(m,l,n)(i,j));
              dhess_Vec[num](j,i) = scale*(dhess_dgamma(l,m,n)(i,j) + nodiag*dhess_dgamma(m,l,n)(i,j));
            }
          }
          num++;
        }
    assert (num == dval_Vec.size());
    // for (int i=0; i<dval_Vec.size(); i++)
    // 	fprintf (stderr, "dval_Vec[%d] = %12.6e\n", i, dval_Vec[i]);
    ///////////////////////////////////////////
    // Now, compensate for constraint matrix //
    ///////////////////////////////////////////
    std::fill (d_vals.begin(), d_vals.end(), 0.0);
    int var = 0;
    for (int i=0; i<NumGamma; i++)
      if (IndepVar[i])
      {
        d_vals[var]  = dval_Vec[i];
        d_grads[var] = dgrad_Vec[i];
        d_hess[var]  = dhess_Vec[i];
        var++;
      }
    int constraint = 0;
    for (int i=0; i<NumGamma; i++)
    {
      if (!IndepVar[i])
      {
        int indep_var = 0;
        for (int j=0; j<NumGamma; j++)
          if (IndepVar[j])
          {
            d_vals[indep_var]  -= ConstraintMatrix(constraint,j) * dval_Vec[i];
            d_grads[indep_var] -= ConstraintMatrix(constraint,j) * dgrad_Vec[i];
            d_hess[indep_var]  -= ConstraintMatrix(constraint,j) * dhess_Vec[i];
            indep_var++;
          }
          else
            if (i != j)
              assert (std::fabs(ConstraintMatrix(constraint,j)) < 1.0e-10);
        constraint++;
      }
    }
    return true;
#ifdef DEBUG_DERIVS
    evaluateDerivativesFD(r_12, r_1I, r_2I, d_valsFD, d_gradsFD, d_hessFD);
    fprintf (stderr, "Param   Analytic   Finite diffference\n");
    for (int ip=0; ip<Parameters.size(); ip++)
      fprintf (stderr, "  %3d  %12.6e  %12.6e\n", ip, d_vals[ip], d_valsFD[ip]);
    fprintf (stderr, "Param   Analytic   Finite diffference\n");
    for (int ip=0; ip<Parameters.size(); ip++)
      fprintf (stderr, "  %3d  %12.6e %12.6e   %12.6e %12.6e   %12.6e %12.6e\n", ip,
               d_grads[ip][0], d_gradsFD[ip][0],
               d_grads[ip][1], d_gradsFD[ip][1],
               d_grads[ip][2], d_gradsFD[ip][2] );
    fprintf (stderr, "Param   Analytic   Finite diffference\n");
    for (int ip=0; ip<Parameters.size(); ip++)
      for (int dim=0; dim<3; dim++)
        fprintf (stderr, "  %3d  %12.6e %12.6e   %12.6e %12.6e   %12.6e %12.6e\n", ip,
                 d_hess[ip](0,dim), d_hessFD[ip](0,dim),
                 d_hess[ip](1,dim), d_hessFD[ip](1,dim),
                 d_hess[ip](2,dim), d_hessFD[ip](2,dim) );
#endif
  }

  inline real_type f(real_type r)
  {
    return 0.0;
  }
  inline real_type df(real_type r)
  {
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
      if (cname == "coefficients")
      {
        string type("0"), id("0"), opt("yes");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.add(opt, "optimize");
        cAttrib.put(xmlCoefs);
        notOpt = (opt=="no");
        if (type != "Array")
        {
          PRE.error( "Unknown correlation type " + type +
                     " in PolynomialFunctor3D." + "Resetting to \"Array\"");
          xmlNewProp (xmlCoefs, (const xmlChar*) "type",
                      (const xmlChar*) "Array");
        }
        vector<real_type> params;
        putContent(params, xmlCoefs);
        if (params.size() == Parameters.size())
          Parameters = params;
        else
        {
          app_log() << "Expected " << Parameters.size() << " parameters,"
                    << " but found " << params.size()
                    << " in PolynomialFunctor3D.\n";
          if (params.size()!=0)
            abort(); //you think you know what they should be but don't.
        }
        // Setup parameter names
        int index=0;
        for (int i=0; i<Parameters.size(); i++)
        {
          std::stringstream sstr;
          sstr << id << "_" << i;;
          if(!notOpt)
            myVars.insert(sstr.str(),Parameters[i],optimize::LOGLINEAR_P,true);
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
        if(!notOpt)
        {
          app_log() << "Parameter     Name      Value\n";
          myVars.print(app_log());
        }
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset_gamma();
    //print();
    return true;
  }

  void resetParameters(const opt_variables_type& active)
  {
    if (notOpt)
      return;
    for(int i=0; i<Parameters.size(); ++i)
    {
      int loc=myVars.where(i);
      if(loc>=0)
        Parameters[i]=myVars[i]=active[loc];
    }
    if (ResetCount++ == 100)
    {
      ResetCount = 0;
      //print();
    }
    reset_gamma();
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  void print()
  {
    const int N = 100;
    string fname = iSpecies + ".J3.h5";
    hid_t hid = H5Fcreate (fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    Array<real_type,3> val (N,N,N);
    for (int i=0; i<N; i++)
    {
      double r_12 = (real_type)i/(real_type)(N-1);
      for (int j=0; j<N; j++)
      {
        double r_1I = (real_type)j/(real_type)(N-1) * 0.5*cutoff_radius;
        for (int k=0; k<N; k++)
        {
          double r_2I = (real_type)k/(real_type)(N-1) * 0.5*cutoff_radius;
          double rmin = std::fabs(r_1I - r_2I);
          double rmax = std::fabs(r_1I + r_2I);
          double r = rmin + r_12*(rmax-rmin);
          val(i,j,k) = evaluate (r, r_1I, r_2I);
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

  inline int getNumParameters()
  {
    return Parameters.size();
  }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: PolynomialFunctor3D.h 1691 2007-02-01 21:51:50Z jnkim $
 ***************************************************************************/
