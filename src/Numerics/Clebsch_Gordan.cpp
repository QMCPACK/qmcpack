// -*- C++ -*-
/*! \author Jordan Vincent
 *  \author Curry Taylor
 *  \note  The original Prim was written in F90 by Tim Wilkens.
 */

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "Numerics/Clebsch_Gordan.h"

/*!\fn Clebsch_Gordan::Clebsch_Gordan(const int lmax)
 * \param lmax the maximum angular momentum
 * \brief Constructs all the Clebsch-Gordan coefficients
 for \f$ 0 \le l_1 \le l_{max}, 0 \le l_2 \le l_1 \f$. 
 This routine was adopted from a FORTRAN 77 algorithm taken 
 from  Rose's 'Elementary Theory of Angular Momentum', p. 39, 
 Wigner's formula. The coefficients listed are only those 
 for which \f$ l_1 \ge l_2 \f$, otherwise the relationship
 \f[ \langle l_1 m_1 l_2 m_2 | l_3 m_3 \rangle
 = (-1)^{l_1+l_2+l_3} \langle l_2 m_2 l_1 m_1 | l_3 m_3 
 \rangle \f] applies.
 Coeffienents known to be zero (because of either the
 L or M selection rules) are not computed, and these should not
 be sought. 
 * \note The indexing of the array is as follows:
 \f[ cg(l1,l2,l3,m1+Lmax,m2+Lmax) \f]
 this is due to the need to index arrays from 0 to N-1 
 in C++.  
 * \note As the main routine is statically defined (hard-coded),
 \f$ l \f$ values greater than 6 are not allowed and will cause 
 error.
 *
 */

Clebsch_Gordan::Clebsch_Gordan(const int lmax):Lmax(lmax){
 
  int l1 = Lmax+1;
  int l2 = 2*Lmax+1;

  cg.resize(l1,l1,l2,l2,l2);

  // std::cout << cg.shape() << std::endl;

  cg = 0.0;

  build_coefficients();
}

Clebsch_Gordan::~Clebsch_Gordan() { cg.free(); }

/*!
 * \fn void Clebsch_Gordan::build_coefficients()
 *
 * \brief Calculates the Clebsch-Gordan coefficients and stores 
 them in a 5-dimensional Blitz++ array.
 *
 */

void Clebsch_Gordan::build_coefficients() {

  double si[33], fa[33], sum, prefactor;
  int lmin, i, l1, l2, l3, m1, m2, m3, nmin, nmax;

  si[0] = 1.0;
  fa[0] = 1.0;

  for(i=1; i<=32; i++) {
    si[i] = -si[i-1];
    fa[i] = (double)i * fa[i-1];
  }

  for(l1=0; l1<=Lmax; l1++) {
    for(l2=0; l2<=l1; l2++) {
      for(m1=-l1; m1<=l1; m1++) {
        for(m2=-l2; m2<=l2; m2++) {
	  m3 = m1 + m2;
	  lmin = abs(l1-l2);
		  
	  if(lmin < abs(m3)) { lmin = abs(m3); }

	  for(l3=lmin; l3<=l1+l2; l3++) {
	    prefactor = 2.0*l3 + 1.0;
	    prefactor *= fa[l3+l1-l2] / fa[l1+l2+l3+1];
	    prefactor *= fa[l3-l1+l2] / fa[l1-m1];
	    prefactor *= fa[l1+l2-l3] / fa[l1+m1];
	    prefactor *= fa[l3+m3] / fa[l2-m2];
	    prefactor *= fa[l3-m3] / fa[l2+m2];
	    prefactor = sqrt(prefactor);

	    nmax = l3 - l1 + l2;
	    if(l3+m3 < nmax) { nmax = l3+m3; }
			
	    nmin = 0;
	    if(l1-l2-m3 < nmin) { nmin = -(l1-l2-m3); }
	
	    sum = 0;
	    for(i=nmin; i<=nmax; i++) {
	      sum += (si[i+l2+m2]/fa[i]) * fa[l2+l3+m1-i] * fa[l1-m1+i]
		/ fa[l3-l1+l2-i] / fa[l3+m3-i] / fa[i+l1-l2-m3];
	    }

            cg(l1,l2,l3,m1+Lmax,m2+Lmax) = prefactor*sum;
	    cg(l2,l1,l3,m2+Lmax,m1+Lmax) = si[l1+l2+l3]*prefactor*sum;
	  }
	}
      }
    }
  }
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

