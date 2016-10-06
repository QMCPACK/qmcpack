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
    
    



#include "Numerics/Spline3D/CubicSpline.h"
#include <iostream>



void CubicSpline::Update()
{
  blitz::Array<double,1> mu(n_x);
  F(0)[1] = 1.5 * ( F(1)[0] - F(0)[0] ) * invh;
  F(n_x-1)[1] = 1.5 * ( F(n_x-1)[0] - F(n_x-2)[0] ) * invh;
  mu(0) = 0.5;
  for(int i = 1; i < n_x-1; i++)
  {
    double lambda = 0.25;
    mu(i) = 0.5 - lambda;
    F(i)[1] = 3.0* invh * ( lambda * ( F(i)[0] - F(i-1)[0] ) +
                            mu(i) * ( F(i+1)[0] - F(i)[0] ) );
    double ci = 1.0 - lambda * mu(i-1);
    F(i)[1] -= lambda * F(i-1)[1];
    mu(i) /= ci;
    F(i)[1] /= ci;
  }
  int i = n_x - 1;
  double lambda = 0.5;
  mu(i) = 0.5 - lambda;
  F(i)[1] = 3.0 * invh * lambda * ( F(i)[0] - F(i-1)[0] );
  double ci = 1.0 - lambda * mu(i-1);
  F(i)[1] -= lambda * F(i-1)[1];
  mu(i) /= ci;
  F(i)[1] /= ci;
  for( i = n_x - 2; i >= 0; i-- )
    F(i)[1] -= mu(i) * F(i+1)[1];
  for(int j = 0; j < n_x; j++)
    std::cout << F(j)[0] << '\t' << F(j)[1] << std::endl;
  exit(-1);
  UpToDate = true;
  return;
}

void CubicSpline::Update(double dfi, double dff)
{
  blitz::Array<double,1> mu(n_x);
  F(0)[1] = dfi;
  F(n_x-1)[1] = dff;
  mu(0) = 0.5;
  for(int i = 1; i < n_x-1; i++)
  {
    double lambda = 0.25;
    mu(i) = 0.5 - lambda;
    F(i)[1] = 3.0* invh * ( lambda * ( F(i)[0] - F(i-1)[0] ) +
                            mu(i) * ( F(i+1)[0] - F(i)[0] ) );
    double ci = 1.0 - lambda * mu(i-1);
    F(i)[1] -= lambda * F(i-1)[1];
    mu(i) /= ci;
    F(i)[1] /= ci;
  }
  for( int i = n_x - 2; i > 0; i-- )
    F(i)[1] -= mu(i) * F(i+1)[1];
  UpToDate = true;
  return;
}






/*
void CubicSpline::Update(){

  blitz::Array<double,1> a(n_x),b(n_x),c(n_x),r(n_x);

  /// initilaise the coefficient matrix
  a(0) = 0.0; b(0) = 2.0; b(n_x-1) = 2; c(n_x-1) = 0.0;
  for(int i = 1; i < n_x-1; i++){
    a(i) = 1.0; b(i) = 4.0; c(i) = 1.0;
  }

  /// initialise the right hand side
  double fac = 3* invh ;
  r(0) = fac*(F(1)[0] - F(0)[0]) - 0.5*d2i;
  for(int i = 1; i < n_x-1; i++) r(i) = fac*(F(i+1)[0]-F(i-1)[0]);
  r(n_x-1) = fac*(F(n_x-1)[0] - F(n_x-2)[0]) + 0.5*d2f;


  blitz::Array<double,1> gamma(n_x);

  /// decomposition and forward substitution
  double beta=b(0);
  F(0)[1] = r(0)/beta;
  for(int j = 1; j < n_x; j++){
    gamma(j) = c(j-1)/beta;
    beta = b(j) - a(j)*gamma(j);
    F(j)[1] = (r(j)-a(j)*F(j-1)[1])/beta;
  }
  /// backsubstitution
  for(int j = n_x-2; j >= 0; j--)
    F(j)[1] -= gamma(j+1)*F(j+1)[1];


//   for(int i = 0; i < n_x-1; i++)
//     std::cout << F(i)[0] << '\t' << F(i)[1] << std::endl;
//   exit(-1);

  UpToDate = true;

  return;

}

*/
