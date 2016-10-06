//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef NONLINEAR_FIT_H
#define NONLINEAR_FIT_H

#include "Blitz.h"
#include "blitz/tinyvec.h"
#include "blitz/array.h"
#include "MatrixOps.h"

using namespace blitz;

// Template class for fitting a function with M nonlinear parameters.
template<int M, typename ModelType>
class NonlinearFitClass
{
private:
  ModelType Model;
  TinyVector<double,M> Params, Beta, dParams;
  TinyMatrix<double,M,M> Alpha;
  Array<double,2> AlphaInv;
  double Chi2 (const Array<double,1> &x, const Array<double,1> &y,
	       const Array<double,1> &sigma, TinyVector<double,M> params);
  void CalcAlphaBeta(const Array<double,1> &x, 
		     const Array<double,1> &y,
		     const Array<double,1> &sigma,
		     TinyVector<double,M> params);
  void Solve ();

public:
  void Fit(const Array<double,1> &x, const Array<double,1> &y,
	   const Array<double,1> &sigma, TinyVector<double,M> &params);

  inline Array<double,2>& GetCovariance()
  { return AlphaInv; }

  NonlinearFitClass(ModelType model) :
    Model(model)
  {
    AlphaInv.resize(M,M);
  };
};


template<int M, typename ModelType> double
NonlinearFitClass<M,ModelType>::Chi2(const Array<double,1> &x, 
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> params)
{
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);
  Model.SetParams(params);

  double chi2 = 0;
  for (int i=0; i<N; i++) {
    double val = Model (x(i));
//     std::cerr << "val = " << val << std::endl; 
//     std::cerr << "yi  = " << y(i) << std::endl;
//     std::cerr << "sigma = " << sigma(i) << std::endl;
    chi2 += (val-y(i))*(val-y(i))/(sigma(i)*sigma(i));
  }
  return chi2;
}

template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::CalcAlphaBeta (const Array<double,1> &x, 
					       const Array<double,1> &y,
					       const Array<double,1> &sigma,
					       TinyVector<double,M> params)
{
  Model.SetParams(params);
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);
  
  for (int k=0; k<M; k++) {
    Beta[k] = 0.0;
    for (int l=0; l<M; l++)
      Alpha(k,l) = 0.0;
  }

  for (int i=0; i<N; i++) {
    double val = Model(x(i));
    TinyVector<double,M> grad = Model.Grad(x(i));
    for (int k=0; k<M; k++)
      Beta[k] += grad[k]*(y(i)-val)/(sigma(i)*sigma(i));
  }
  for (int i=0; i<N; i++) {
    double val = Model(x(i));
    TinyVector<double,M> grad = Model.Grad(x(i));
    for (int k=0; k<M; k++)
      for (int l=0; l<M; l++)
	Alpha(k,l) += grad[k]*grad[l]/(sigma(i)*sigma(i));
  }
  // std::cerr << "Alpha = " << Alpha << std::endl;
}


template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::Solve()
{
  for (int i=0; i<M; i++) {
    dParams[i] = 0.0;
    for (int j=0; j<M; j++)
      AlphaInv(i,j) = Alpha(i,j);
  }

  GJInverse(AlphaInv);
  TinyMatrix<double,M,M> id;
  id = 0.0;
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      for (int k=0; k<M; k++) 
	id(i,j) += AlphaInv(i,k)*Alpha(k,j);

  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      dParams[i] += AlphaInv(i,j)*Beta(j);

  // Do simple Gaussian elimination
//   double dInv = 1.0/Alpha(0,0);
//   for (int col=0; col<M; col++)
//     Alpha(0,col) *= dInv;
//   Beta(0) *= dInv;
//   for (int row=1; row<M; row++) {
//     for (int col=0; col<row; col++) {
//       Beta(row) -= Alpha(row,col)*Beta(col);
//       for (int k=0;
//     double diag = Alpha(row,row);
//     double diagInv = 1.0/diag;
//     for (int col=0; col<M; col++)
//     Beta(row) *= diagInv;
//     for (int col=row+1; col<


}


template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::Fit (const Array<double,1> &x,
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> &params)
{
  double chiNow = Chi2 (x, y, sigma, params);
  double lambda = 0.01;
  
  bool done = false;
  int iter = 1;
  int numSmallDecrease = 0;
  while (!done) {
//     std::cerr << "Iteration " << iter << ":  Chi2 = " << chiNow << std::endl;
//    std::cerr << "params = " << params << std::endl;
//     std::cerr << "lambda = " << lambda << std::endl;

    CalcAlphaBeta (x, y, sigma, params);
    for (int i=0; i<M; i++)
      Alpha(i,i) *= (1.0+lambda);
    Solve();
//     done = true;
//     for (int i=0; i<M; i++)
//       if (std::abs(Beta[i]) > 1.0e-6)
// 	done = false;
    if (!done) {
      double chiNew = Chi2 (x, y, sigma, params + dParams);
      if (chiNew <= chiNow) {
	lambda *= 0.5;
	params = params + dParams;
	if ((chiNow - chiNew) < 0.001)
	  numSmallDecrease++;
	if (numSmallDecrease > 4)
	  done = true;
	chiNow = chiNew;
      }
      else {
	lambda *= 2.0;
      }
    }
    iter++;
  }
  CalcAlphaBeta (x, y, sigma, params);
  for (int i=0; i<M; i++) 
    for (int j=0; j<M; j++)
      AlphaInv(i,j) = Alpha(i,j);
  GJInverse(AlphaInv);
//   std::cerr << "Covariace matrix:\n";
//   for (int i=0; i<AlphaInv.rows(); i++) {
//     for (int j=0; j<AlphaInv.cols(); j++)
//       fprintf (stderr, "%12.4e ", AlphaInv(i,j));
//     fprintf (stderr, "\n");
//   }
  Model.SetParams (params);
  //  std::cerr << "Chi2 = " << chiNow << std::endl;
}

#endif
