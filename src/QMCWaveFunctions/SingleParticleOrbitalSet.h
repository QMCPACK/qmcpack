//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef OHMMS_SINGLEPARTICLEORBITALSET_H
#define OHMMS_SINGLEPARTICLEORBITALSET_H
#include <vector>
/**class a set of single-particle orbitals. 
 *@brief This class provides necessary interfaces for SlaterDeterminant<SPOSet>
 * and the example can be found in examples/simpledet.cpp
 * The template parameter is any analytic function with a member function
 * which can provde the value, gradient and laplacian at a point.
 * value_type OT::evaluate(const point_type&, gradient_type&, value_type& ) 
 * Example classes can be found in Numerics/CosineFunction.h
 */

namespace ohmmsqmc {

template<class OT>
struct SingleParticleOrbitalSet {

  ///the type of single-particle orbtials 
  typedef OT                      SPOrbital_t;
  typedef typename OT::value_type value_type;

  vector<SPOrbital_t*> Phi;

  ///constructor
  SingleParticleOrbitalSet(){ }

  /**add a single-particle orbital */
  int add(SPOrbital_t* afunction ) {
    Phi.push_back(afunction);
    return Phi.size()-1;
  }

  void reset() {
    for(int i=0; i<Phi.size(); i++) Phi[i]->reset();
  }

  void resizeByWalkers(int ) { }

  inline int size() const { return Phi.size();}

  template<class PTCL, class VM, class GM>
  inline void evaluate(const PTCL& P, int first, int last,
		       VM& logdet, GM& dlogdet, VM& d2logdet) {
    int n = last-first;
    int iat = first;
    for(int i=0; i<n; i++,iat++) {
      Phi[0]->set_point(P.R[iat]);
      for(int j=0; j<n; j++) {
	logdet(j,i)= Phi[j]->evaluate(P.R[iat], dlogdet(i,j),d2logdet(i,j));
      }
    }
  }

  template<class WREF, class VM, class GM>
  inline void 
  evaluate(const WREF& W, int first, int last,
	   vector<VM>& logdet, vector<GM>& dlogdet, vector<VM>& d2logdet) {
   int n = last-first;
    for(int i=0; i<n; i++){
      for(int iw=0; iw<W.walkers(); iw++) {
	int jat = first;
	for(int j=0; j<n; j++,jat++) {
	  logdet[iw](i,j) 
	    = Phi[i]->evaluate(W.R(iw,jat),dlogdet[iw](i,j), d2logdet[iw](i,j));
	}
     }
   }
  }

//   template<class PT, class GT>
//   inline void evaluate(const PT* R, 	
// 		       value_type* logdet, GT* dlogdet,value_type* d2logdet,
// 		       int first, int last, int iw, int nw) {
//     int ij = 0;
//     int n = last-first;
//     int iat = nw*first+iw;
//     for(int i=0; i<n; i++,iat+=nw) {
//       for(int j=0; j<n; j++,ij++) {
// 	logdet[j*n+i] = Phi[j]->evaluate(R[iat], dlogdet[ij],d2logdet[ij]);
//       }
//     }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
