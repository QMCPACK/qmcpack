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
#ifndef OHMMS_QMC_SYMMETRICDISTANCETABLEDATAIMPL_H
#define OHMMS_QMC_SYMMETRICDISTANCETABLEDATAIMPL_H

namespace ohmmsqmc {

  /**@ingroup nnlist
   * @brief A derived classe from DistacneTableData, specialized for dense-symmetric case
   *
   * SymmetricDTD stands for Symmetric Distance Table Data.
   * The source and target sets are identical and the sum is over 
   * distict pairs, i.e., \f$\sum_{i}\sum_{j<i}\f$.
   * The template parameter BC provides BC::apply member function
   * to evaluate Cartesian distances.
   */
  template<class BC>
  struct SymmetricDTD: public DistanceTableData{
 
    //blitz::Array<IndexType,2> IJ;
    ///constructor using source and target arrays
    SymmetricDTD(const ParticleSet& source, const ParticleSet& target):
      DistanceTableData(source,target){ }

    void create(int walkers){
      int nw = (walkers>0)? walkers:1;
      reset(Origin.getTotalNum(),nw);
    }

    inline void reset(int m, int nactive) {
      if(m != N[SourceIndex] || nactive != N[WalkerIndex]) {
	N[SourceIndex]=m;
	N[VisitorIndex]=m;
	int nn = m*(m-1)/2;
	M.resize(m+1);
	J.resize(nn);
	IJ.resize(m*m);
	PairID.resize(m*m);	
	resize(nn,nactive);
	M[0] = 0;
	int nsp = Origin.groups();

	int ij = 0;
	for(int i=0; i<m; i++) {
	  for(int j=i+1; j<m; j++, ij++) {
	    J[ij] = j;
	    PairID[ij] = Origin.GroupID[j]+nsp*Origin.GroupID[i];
	    IJ[i*m+j] = ij;
	    IJ[j*m+i] = ij;
	    //@warning: using a simple pair scheme
	    // Upper packed-storage scheme
	    //PairID[ij] = Origin.GroupID[j]+(2*nsp-Origin.GroupID[i]-1)*(Origin.GroupID[i])/2;
	  }
	  M[i+1] = ij;
	}
        npairs_m = ij;
      }
    }


    ///evaluate the Distance Table using a set of Particle Positions
    inline void evaluate(const WalkerSetRef& W) {

      int copies = W.walkers();
      int visitors = W.particles();

      for(int iw=0; iw<copies; iw++) {
	int nn=0;
	for(int i=0; i<visitors-1; i++) {
	  PosType a0 = W.R(iw,i);
	  for(int j=i+1; j<visitors; j++, nn++) {
	    PosType drij = W.R(iw,j)-a0;
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
      int n = P.getTotalNum();
      //reset(n,1);
      int ij = 0;
      for(int i=0; i<n; i++) {
	for(int j=i+1; j<n; j++, ij++) {
	  PosType drij = P.R[j]-P.R[i];
	  RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
	  r_m[ij] = sep;
	  rinv_m[ij] = 1.0/sep;      
	  dr_m[ij] = drij;
	}
      }
    }

    ///evaluate the temporary pair relations
    inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat) {

      activePtcl=jat;
      Temp[jat].reset();
      for(int iat=0; iat<jat; iat++) {
	int loc = IJ[iat*N[SourceIndex]+jat];
	PosType drij = rnew - P.R[iat];
	RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
	Temp[iat].r1=sep;
	Temp[iat].rinv1=1.0/sep;
	Temp[iat].dr1=drij;
	Temp[iat].r0=r_m[loc];
	Temp[iat].rinv0=rinv_m[loc];
	Temp[iat].dr0=-1.0*dr_m[loc];
      }
      for(int iat=jat+1,nn=jat; iat< N[SourceIndex]; iat++) {
	int loc = IJ[iat*N[SourceIndex]+jat];
	PosType drij = rnew - P.R[iat];
	RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
	Temp[iat].r1=sep;
	Temp[iat].rinv1=1.0/sep;
	Temp[iat].dr1=drij;
	Temp[iat].r0=r_m[loc];
	Temp[iat].rinv0=rinv_m[loc];
	Temp[iat].dr0=dr_m[loc];
      }
    }

    ///update the stripe for jat-th particle
    inline void update(IndexType jat) {

      int nn=jat;
      for(int iat=0;iat<jat; iat++,nn+=N[SourceIndex]) {
	int loc =IJ[nn];
	r_m[loc] = Temp[iat].r1;
	rinv_m[loc]= Temp[iat].rinv1;
	dr_m[loc]= -1.0*Temp[iat].dr1;
      }

      for(int nn=M[jat]; nn<M[jat+1]; nn++) {
	int iat = J[nn];
	r_m[nn] = Temp[iat].r1;
	rinv_m[nn]= Temp[iat].rinv1;
        dr_m[nn]= Temp[iat].dr1;
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
