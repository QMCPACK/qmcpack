//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMSHF_SJPSEUDOPOTENTIAL_H
#define OHMMSHF_SJPSEUDOPOTENTIAL_H
#include "SQD/SphericalPotential/RadialPotential.h"
#include "Optimize/VarList.h"
namespace ohmmshf {

  /**
     @ingroup RadialPotential
     @class SJPseudoPotential
     @brief implements the Starkloff-Joannopoulos pseudopotential
     *
     \f[
     V_{SJ}(r) = -\frac{Z_{Eff}}{r}\frac{1-e^{-\lambda r}}
     {1+e^{-\lambda (r-r_c)}} 
     \f]   
     *
     See Th. Starkloff and J.D. Joannopoulos, Phys. Rev. B, \textbf{16},
     5212, (1977).
  */
  struct SJPseudoPotential: public RadialPotentialBase {
    value_type Zeff, SJ_lambda, rc; 
    SJPseudoPotential(VarRegistry<value_type>&,value_type,
		      value_type,value_type);
    SJPseudoPotential(value_type,value_type,value_type);
    value_type evaluate(const BasisSetType& psi,
			RadialOrbitalSet_t&V, int norb);
    ///must always return 1
    inline int getNumOfNodes(int n, int l);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
