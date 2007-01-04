//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_WALKERSETREF_H
#define QMCPLUSPLUS_WALKERSETREF_H

#include "Particle/ParticleSet.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus {

  /** class to assist vectorized operations on Walkers
   */
  struct WalkerSetRef: public QMCTraits {

    /** enum for the particle and walker indices
     *@brief Introduced to handle fast access to Particles or Walkers
     *indices.
     */
    enum {Particles =0, Walkers };

    /**
     *<ul>
     *<li>N[Particles] = number of particles
     *<li>N[Walkers] = number of walkers
     */
    TinyVector<int,2> N;
    typedef Matrix<PosType>   WalkerPosition_t;
    typedef Matrix<GradType>  WalkerGradient_t;
    typedef Matrix<ValueType> WalkerLaplacian_t;

    const ParticleSet& PtclRef;
    WalkerPosition_t  R;
    WalkerGradient_t  G;
    WalkerLaplacian_t L;

    ///default constructor
    WalkerSetRef(const ParticleSet& p): PtclRef(p) {N = 0;}

    inline int tag() const { return PtclRef.tag();}
    inline int walkers() const { return N[Walkers];}
    inline int particles() const { return N[Particles];}

    /** resize the containers */
    inline void resize(int nw, int nptcl) {
      R.resize(nw,nptcl);
      G.resize(nw,nptcl);
      L.resize(nw,nptcl);
      N[Particles] = nptcl;
      N[Walkers] = nw;
    }

    /**@param first an iterator of the first walker to work on
     *@param last an iterator of the last walker to work on
     *@param tauinv timestep
     *@brief Make Metropolis move to the walkers and save in a temporary array.
     */
    template<class IT>
    inline void sample(IT first, IT last, RealType tauinv) {
      makeGaussRandom(R);
      R *= tauinv;
      ///reset the number of active walkers to zero to count them
      int iw = 0;
      while(first != last) {
	const ParticleSet::ParticlePos_t& r = (*first)->R;
	const ParticleSet::ParticlePos_t& drift = (*first)->Drift;
	for(int jat=0; jat<N[Particles]; jat++) {
	  R(iw,jat) += r[jat]+drift[jat];
	}
	///increment the index and iterator
	iw++; first++;
      }
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
