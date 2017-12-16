//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef CLEBSCH_GORDAN_H2
#define CLEBSCH_GORDAN_H2

#include <blitz/array.h>

// -*-  -*-
/*! \author Jordan Vincent
 *  \author Curry Taylor
 *  \note  The original Prim was written in F90 by Tim Wilkens.
 */


/**class Clebsch_Gordan
 *\brief Calculates the Clebsch-Gordan coefficients
 */
class Clebsch_Gordan
{
public:

  Clebsch_Gordan(const int lmax);
  ///destructor
  ~Clebsch_Gordan();

  ///maximum angular momentum
  int Lmax;
  ///array to store the Clebsch-Gordan coefficients
  blitz::Array<double,5> cg;
  ///returns \f$c_g(l1,l2,l3,m1,m2) \f$
  inline double operator()(int l1, int l2, int l3, int m1, int m2) const
  {
    return cg(l1,l2,l3,m1+Lmax,m2+Lmax);
  }

private:
  /// default constructor not implemented
  Clebsch_Gordan() { }

  void build_coefficients();
};

#endif



