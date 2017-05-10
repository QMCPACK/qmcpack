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
    
    



#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "Numerics/Clebsch_Gordan.h"

/**
 *@param lmax the maximum angular momentum
 *@brief Constructs all the Clebsch-Gordan coefficients
 for \f$ 0 \le l_1 \le l_{max}, 0 \le l_2 \le l_1 \f$.
 *
 This routine was adopted from a FORTRAN 77 algorithm taken
 from  Rose's 'Elementary Theory of Angular Momentum', p. 39,
 Wigner's formula.
 Coeffienents known to be zero (because of either the
 L or M selection rules) are not computed, and these should not
 be sought.
 *@note The indexing of the vector is as follows:
 \f[ index = i + jI + kIJ + lIJK + mIJKL \f]
 where \f$ (I,J,K,L,M) \f$ are the maximum extent
 of the Clebsch-Gordon coefficients.
 *@note States with \f$ l > 6 \f$ are not allowed and will
 result in an error.
 *
 */

Clebsch_Gordan::Clebsch_Gordan(const int lmax):
  Lmax(lmax), L1max(lmax+1), L2max(2*lmax+1)
{
  CG_coeff.resize(L1max*L1max*L2max*L2max*L2max);
  for(int i=0; i<CG_coeff.size(); i++)
    CG_coeff[i] = 0.0;
  double si[33], fa[33], sum, prefactor;
  int lmin, i, l1, l2, l3, m1, m2, m3, nmin, nmax;
  si[0] = 1.0;
  fa[0] = 1.0;
  for(i=1; i<=32; i++)
  {
    si[i] = -si[i-1];
    fa[i] = static_cast<double>(i) * fa[i-1];
  }
  for(l1=0; l1<=Lmax; l1++)
  {
    for(l2=0; l2<=l1; l2++)
    {
      for(m1=-l1; m1<=l1; m1++)
      {
        for(m2=-l2; m2<=l2; m2++)
        {
          m3 = m1 + m2;
          lmin = std::abs(l1-l2);
          if(lmin < std::abs(m3))
          {
            lmin = std::abs(m3);
          }
          for(l3=lmin; l3<=l1+l2; l3++)
          {
            prefactor = 2.0*l3 + 1.0;
            prefactor *= fa[l3+l1-l2] / fa[l1+l2+l3+1];
            prefactor *= fa[l3-l1+l2] / fa[l1-m1];
            prefactor *= fa[l1+l2-l3] / fa[l1+m1];
            prefactor *= fa[l3+m3] / fa[l2-m2];
            prefactor *= fa[l3-m3] / fa[l2+m2];
            prefactor = sqrt(prefactor);
            nmax = l3 - l1 + l2;
            if(l3+m3 < nmax)
            {
              nmax = l3+m3;
            }
            nmin = 0;
            if(l1-l2-m3 < nmin)
            {
              nmin = -(l1-l2-m3);
            }
            sum = 0;
            for(i=nmin; i<=nmax; i++)
            {
              sum += (si[i+l2+m2]/fa[i]) * fa[l2+l3+m1-i] * fa[l1-m1+i]
                     / fa[l3-l1+l2-i] / fa[l3+m3-i] / fa[i+l1-l2-m3];
            }
            /* apply the relationship
               \langle l_1 m_1 l_2 m_2 | l_3 m_3 \rangle
               = (-1)^{l_1+l_2+l_3} \langle l_2 m_2 l_1 m_1 | l_3 m_3 */
            CG_coeff[index(l1,l2,l3,m1,m2)] = prefactor*sum;
            CG_coeff[index(l2,l1,l3,m2,m1)] = si[l1+l2+l3]*prefactor*sum;
          }
        }
      }
    }
  }
}



