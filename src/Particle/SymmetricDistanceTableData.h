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
#ifndef QMCPLUSPLUS_SYMMETRICDISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_SYMMETRICDISTANCETABLEDATAIMPL_H

namespace qmcplusplus {

  /**@ingroup nnlist
   * @brief A derived classe from DistacneTableData, specialized for dense-symmetric case
   *
   * SymmetricDTD stands for Symmetric Distance Table Data.
   * The source and target sets are identical and the sum is over 
   * distict pairs, i.e., \f$\sum_{i}\sum_{j<i}\f$.
   * The template parameter BC provides BC::apply member function
   * to evaluate Cartesian distances.
   */
  template<typename T, unsigned D, int SC>
  struct SymmetricDTD
  : public DTD_BConds<T,D,SC>, public DistanceTableData
  {
 

    //blitz::Array<IndexType,2> IJ;
    ///constructor using source and target arrays
    SymmetricDTD(const ParticleSet& source, const ParticleSet& target)
      : DTD_BConds<T,D,SC>(source.Lattice), DistanceTableData(source,target)
    { 
        create(1);
    }

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
   // inline void evaluate(const WalkerSetRef& W) {
   //   int copies = W.walkers();
   //   int visitors = W.particles();
   //   for(int iw=0; iw<copies; iw++) {
   //     int nn=0;
   //     for(int i=0; i<visitors-1; i++) {
   //       PosType a0 = W.R(iw,i);
   //       for(int j=i+1; j<visitors; j++, nn++) {
   //         PosType drij = W.R(iw,j)-a0;
   //         RealType sep = std::sqrt(BC::apply(Origin.Lattice,drij));
   //#ifdef USE_FASTWALKER
   //         r2_m(nn,iw) = sep;
   //         rinv2_m(nn,iw) = 1.0/sep;
   //         dr2_m(nn,iw) = drij;
   //#else
   //         r2_m(iw,nn) = sep;
   //         rinv2_m(iw,nn) = 1.0/sep;
   //         dr2_m(iw,nn) = drij;
   //#endif
   //       }
   //     }
   //   }
   // }

    inline void evaluate(const ParticleSet& P)
    {
      const int n = N[SourceIndex];
      for(int i=0,ij=0; i<n; i++) 
        for(int j=i+1; j<n; j++, ij++) 
          dr_m[ij]=P.R[j]-P.R[i];
      //old with static type
      //BC::apply(Origin.Lattice,dr_m,r_m,rinv_m);
      DTD_BConds<T,D,SC>::apply_bc(dr_m,r_m,rinv_m);
    }

    ///evaluate the temporary pair relations
    inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat) {
      activePtcl=jat;
      for(int iat=0; iat<N[SourceIndex]; ++iat)
      {
        PosType drij(rnew - P.R[iat]);
        Temp[iat].dr1_nobox=drij;
        RealType sep=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
        Temp[iat].r1=sep;
        Temp[iat].rinv1=1.0/sep;
        Temp[iat].dr1=drij;
      }

      //for(int iat=0; iat<jat; iat++) {
      //  int loc = IJ[iat*N[SourceIndex]+jat];
      //  PosType drij(rnew - P.R[iat]);
      //  //old with static type
      //  //RealType sep=std::sqrt(BC::apply(Origin.Lattice,drij));
      //  RealType sep=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      //  Temp[iat].r1=sep;
      //  Temp[iat].rinv1=1.0/sep;
      //  Temp[iat].dr1=drij;
      //  //Temp[iat].r0=r_m[loc];
      //  //Temp[iat].rinv0=rinv_m[loc];
      //  //Temp[iat].dr0=-1.0*dr_m[loc];
      //}
      //Temp[jat].reset();
      //for(int iat=jat+1,nn=jat; iat< N[SourceIndex]; iat++) {
      //  int loc = IJ[iat*N[SourceIndex]+jat];
      //  PosType drij(rnew - P.R[iat]);
      //  //old with static type
      //  //RealType sep=std::sqrt(BC::apply(Origin.Lattice,drij));
      //  RealType sep=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      //  Temp[iat].r1=sep;
      //  Temp[iat].rinv1=1.0/sep;
      //  Temp[iat].dr1=drij;
      //  //Temp[iat].r0=r_m[loc];
      //  //Temp[iat].rinv0=rinv_m[loc];
      //  //Temp[iat].dr0=dr_m[loc];
      //}
    }

    ///evaluate the temporary pair relations
    inline void moveby(const ParticleSet& P, const PosType& displ, IndexType iat) 
    {
      activePtcl=iat;
      for(int jat=0; jat<iat; ++jat) 
        temp_dr[jat]=-1.0*(displ+dr_m[IJ[jat*N[SourceIndex]+iat]]);
      temp_dr[iat]=0.0;
      for(int jat=iat+1; jat< N[SourceIndex]; ++jat) 
        temp_dr[jat]=dr_m[IJ[jat*N[SourceIndex]+iat]]-displ;
      DTD_BConds<T,D,SC>::apply_bc(temp_dr,temp_r);
      //BC::get_min_distanceX(P.Lattice,Temp);
      //for(int jat=0; jat<iat; ++jat) 
      //{
      //  Temp[jat].dr1=-1.0*(displ+dr_m[IJ[jat*N[SourceIndex]+iat]]);
      //  Temp[jat].r1=std::sqrt(BC::get_min_distance(Origin.Lattice,Temp[jat].dr1,Rmax2));
      //  //Temp[jat].rinv1=1.0/Temp[jat].r1;
      //}
      //Temp[iat].reset();
      //for(int jat=iat+1; jat< N[SourceIndex]; ++jat) 
      //{
      //  Temp[jat].dr1=dr_m[IJ[jat*N[SourceIndex]+iat]]-displ;
      //  Temp[jat].r1=std::sqrt(BC::get_min_distance(Origin.Lattice,Temp[jat].dr1,Rmax2));
      //  //Temp[jat].rinv1=1.0/Temp[jat].r1;
      //}
    }

    ///evaluate the temporary pair relations
    inline void moveOnSphere(const ParticleSet& P, const PosType& displ, IndexType iat) 
    {
      activePtcl=iat;
      for(int jat=0; jat<iat; ++jat) {
	PosType& drij=Temp[jat].dr1=-1.0*(displ+dr_m[IJ[jat*N[SourceIndex]+iat]]);
	Temp[jat].r1=std::sqrt(dot(drij,drij));
	//Temp[jat].rinv1=1.0/Temp[jat].r1;
      }
      Temp[iat].reset();
      for(int jat=iat+1; jat< N[SourceIndex]; ++jat) {
	PosType& drij=Temp[jat].dr1=dr_m[IJ[jat*N[SourceIndex]+iat]]-displ;
	Temp[jat].r1=std::sqrt(dot(drij,drij));
	//Temp[jat].rinv1=1.0/Temp[jat].r1;
      }
    }


    ///update the stripe for jat-th particle
    inline void update(IndexType jat) {
      int nn=jat;
      for(int iat=0;iat<jat; iat++,nn+=N[SourceIndex]) {
	int loc =IJ[nn];
	r_m[loc] = Temp[iat].r1;
	rinv_m[loc]= 1.0/Temp[iat].r1;
	//rinv_m[loc]= Temp[iat].rinv1;
	dr_m[loc]= Temp[iat].dr1;
      }

      for(int nn=M[jat]; nn<M[jat+1]; nn++) {
	int iat = J[nn];
	r_m[nn] = Temp[iat].r1;
	rinv_m[nn]= 1.0/Temp[iat].r1;
	//rinv_m[nn]= Temp[iat].rinv1;
        dr_m[nn]= -1.0*Temp[iat].dr1;
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
