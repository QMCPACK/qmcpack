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
#ifndef OHMMS_QMC_ASYMMETRICDISTANCETABLEDATAIMPL_H
#define OHMMS_QMC_ASYMMETROCDISTANCETABLEDATAIMPL_H

namespace ohmmsqmc {
  
  /**@ingroup nnlist
   * @brief A derived classe from DistacneTableData, specialized for dense-asymmetric case
   *
   * AsymmetricDTD stands for Asymmetric Distance Table Data with
   * distinct source and target sets.
   * The template parameter BC provides BC::apply member function
   * to evaluate Cartesian distances.
   */
  template<class BC>
  struct AsymmetricDTD: public DistanceTableData {

    const ParticleSet& Target;

    AsymmetricDTD(const ParticleSet& source, 
		  const ParticleSet& target): DistanceTableData(source,target),
					       Target(target){
      create(1);
    }

    void create(int walkers){
      int nw = (walkers>0)? walkers:1;
      reset(Origin.getTotalNum(),Target.getTotalNum(),nw);
    }

    /*!\fn  void reset(int n1, int n2, int nactive){
     *\param n1 the number of sources
     *\param n2 the number of targets
     *\param nactive the number of copies of the targets
     *\brief Resize the internal data and assign indices
     */
     inline void reset(int n1, int n2, int nactive){
      if( n1!=N[SourceIndex] || n2 != N[VisitorIndex] 
	  || nactive != N[WalkerIndex]) {
	N[SourceIndex] = n1;
	N[VisitorIndex] = n2;
	int m = n1*n2;
	if(m) { 
	  M.resize(n1+1);
	  J.resize(m);
	  PairID.resize(m);
	  resize(m,nactive);
	  M[0] = 0; int ij = 0;
	  for(int i=0; i<n1; i++) {
	    for(int j=0; j<n2; j++, ij++) {
	      J[ij] = j;
	      PairID[ij] = Origin.GroupID[i];
	    }
	    M[i+1] = M[i]+n2;
	  }
	  npairs_m = n1*n2;
	}
      }
    }

    ///evaluate the Distance Table using a set of Particle Positions
    inline void evaluate(const WalkerSetRef& W) {
      int copies = W.walkers();
      int visitors = W.particles();
      int ns = Origin.getTotalNum();

      reset(ns,visitors,copies);
      for(int iw=0; iw<copies; iw++) {
	int nn=0;
	for(int i=0; i<ns; i++) {
	  PosType r0 = Origin.R(i);
	  for(int j=0; j<visitors; j++,nn++) {
	    PosType drij = W.R(iw,j)-r0;
	    RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
#ifdef USE_FASTWALKER
	    r2_m(nn,iw) = sep;
	    rinv2_m(nn,iw) = 1.0/sep;
	    dr2_m(nn,iw) = drij;
#else
	    r2_m(iw,nn) = sep;
	    rinv2_m(iw,nn) = 1.0/sep;
	    dr2_m(iw,nn) = drij;
#endif
	  }
	}
      }
    }

    ///not so useful inline but who knows
    inline void evaluate(const ParticleSet& P){
      //reset(Origin.getTotalNum(),P.getTotalNum(),1);
      int nn=0;
      for(int i=0; i<N[SourceIndex]; i++) {
	PosType r0(Origin.R[i]);
	for(int j=0; j<N[VisitorIndex]; j++,nn++) {
	  PosType drij(P.R[j]-r0);
	  RealType sep2(BC::apply(Origin.Lattice,drij));
	  RealType sep(sqrt(sep2));
	  r_m[nn]    = sep;
	  //rr_m[nn]   = sep2;
	  rinv_m[nn] = 1.0/sep;
	  dr_m[nn]   = drij;
	}
      }
    }

    ///evaluate the temporary pair relations
    inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat) {
      activePtcl=jat;
      for(int iat=0, loc=jat; iat<N[SourceIndex]; iat++,loc+=N[VisitorIndex]) {
	PosType drij(rnew-Origin.R[iat]);
	RealType sep2(BC::apply(Origin.Lattice,drij));
	RealType sep(sqrt(sep2));
	Temp[iat].r1=sep;
	//Temp[iat].rr1=sep2;
	Temp[iat].rinv1=1.0/sep;
	Temp[iat].dr1=drij;
	Temp[iat].r0=r_m[loc];
	Temp[iat].rinv0=rinv_m[loc];
	Temp[iat].dr0=dr_m[loc];
      }
    }

    ///evaluate the temporary pair relations
    inline void moveOnSphere(const ParticleSet& P, const PosType& displ, IndexType jat) {
      //if(activePtcl == jat) {
      //  for(int iat=0; iat<N[SourceIndex]; iat++) {
      //    //PosType drij = rnew-Origin.R[iat];
      //    //RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
      //    PosType drij(displ+Temp[iat].dr0);
      //    RealType sep2 = BC::apply(Origin.Lattice,drij);
      //    RealType sep = sqrt(sep2);
      //    Temp[iat].r1=sep;
      //    //Temp[iat].rr1=sep2;
      //    Temp[iat].rinv1=1.0/sep;
      //    Temp[iat].dr1=drij;
      //  }
      //} else {
      //  move(P,P.R[jat]+displ,jat);
      //}
      move(P,P.R[jat]+displ,jat);
    }

    inline void update(IndexType jat) {
      for(int iat=0,loc=jat; iat<N[SourceIndex]; iat++, loc += N[VisitorIndex]) {
	r_m[loc]=Temp[iat].r1;
	//rr_m[loc]=Temp[iat].rr1;
	rinv_m[loc]=Temp[iat].rinv1;
	dr_m[loc]=Temp[iat].dr1;
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
