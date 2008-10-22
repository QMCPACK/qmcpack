//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"

namespace qmcplusplus {

  NonLocalECPComponent::NonLocalECPComponent(): 
    lmax(0), nchannel(0), nknot(0), Rmax(-1), myRNG(&Random)
    { }

  NonLocalECPComponent::~NonLocalECPComponent() {
    for(int ip=0; ip<nlpp_m.size(); ip++) delete nlpp_m[ip];
  }

  NonLocalECPComponent* NonLocalECPComponent::makeClone()
  {
    NonLocalECPComponent* myclone=new NonLocalECPComponent(*this);
    for(int i=0; i<nlpp_m.size(); ++i)
      myclone->nlpp_m[i]=nlpp_m[i]->makeClone();
    return myclone;
  }
  
  void NonLocalECPComponent::add(int l, RadialPotentialType* pp) {
    angpp_m.push_back(l);
    wgt_angpp_m.push_back(static_cast<RealType>(2*l+1));
    nlpp_m.push_back(pp);
  }

  void NonLocalECPComponent::resize_warrays(int n,int m,int l){
    app_log() << "  NonLocalECPComponent::resize_warrays " << endl;
     psiratio.resize(n);
     vrad.resize(m);
     wvec.resize(m);
     Amat.resize(n*m);
     lpol.resize(l+1,1.0);
     rrotsgrid_m.resize(n);
     nchannel=nlpp_m.size();
     nknot=sgridxyz_m.size();
     //This is just to check
     //for(int nl=1; nl<nlpp_m.size(); nl++) nlpp_m[nl]->setGridManager(false);
     if(lmax) {
       Lfactor1.resize(lmax);
       Lfactor2.resize(lmax);
       for(int nl=0; nl<lmax; nl++) {
         Lfactor1[nl]=static_cast<RealType>(2*nl+1);
         Lfactor2[nl]=1.0e0/static_cast<RealType>(nl+1);
       }
     }
  }

  void NonLocalECPComponent::print(std::ostream& os)
  {
     os << "    Maximum angular mementum = "<<  lmax << endl;
     os << "    Number of non-local channels = " << nchannel << endl;
     for(int l=0; l <nchannel; l++) 
       os << "       l(" << l << ")=" << angpp_m[l] <<endl;
     os << "    Cutoff radius = " << Rmax << endl;
     os << "    Spherical grids and weights: " << endl;
     for(int ik=0;ik<nknot; ik++)
       os << "       " << sgridxyz_m[ik] << setw(20) << sgridweight_m[ik] << endl;
  }
  /** evaluate the non-local potential of the iat-th ionic center
   * @param W electron configuration
   * @param iat ionic index
   * @param psi trial wavefunction
   * @param return the non-local component
   *
   * Currently, we assume that the ratio-only evaluation does not change the state
   * of the trial wavefunction and do not call psi.rejectMove(ieL).
   */
  NonLocalECPComponent::RealType 
  NonLocalECPComponent::evaluate(ParticleSet& W, int iat, TrialWaveFunction& psi) {
    RealType esum=0.0;
    for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++){

      register RealType r(myTable->r(nn));
      if(r>Rmax) continue;

      register RealType rinv(myTable->rinv(nn));
      register PosType  dr(myTable->dr(nn));

      // Compute ratio of wave functions
      for (int j=0; j < nknot ; j++){ 
        PosType deltar(r*rrotsgrid_m[j]-dr);
        W.makeMoveOnSphere(iel,deltar); 
        //W.makeMove(iel,deltar); 
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j];
        W.rejectMove(iel);
        //psi.rejectMove(iel);
      }
      // Compute radial potential
      //int k;
      //RealType rfrac;
      //nlpp_m[0]->locate(r,k,rfrac);
      //for(int ip=0;ip< nchannel; ip++){
      //  vrad[ip]=nlpp_m[ip]->f(k,rfrac)*wgt_angpp_m[ip];
      //}
      for(int ip=0;ip< nchannel; ip++){
        vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];
      }

      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
        RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
        // Forming the Legendre polynomials
        lpol[0]=1.0;
        RealType lpolprev=0.0;
        for (int l=0 ; l< lmax ; l++){
          //Not a big difference
          //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
          //lpol[l+1]/=(l+1);
          lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
          lpol[l+1]*=Lfactor2[l]; 
          lpolprev=lpol[l];
        }
        for(int l=0; l <nchannel; l++,jl++) Amat[jl]=lpol[ angpp_m[l] ]; 
      } 

      if(nchannel==1) {
        esum += vrad[0]*BLAS::dot(nknot, &Amat[0],&psiratio[0]);
      } else {
        BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
        esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
      }
      ////////////////////////////////////
      //Original implmentation by S. C.
      //const char TRANS('T');
      //const int ione=1;
      //const RealType one=1.0;
      //const RealType zero=0.0;
      //dgemv(TRANS,nknot,nchannel,one,&Amat[0],nknot,&psiratio[0],ione,zero,&wvec[0],ione);
      //esum += ddot(nchannel,&vrad[0],ione,&wvec[0],ione);
      ////////////////////////////////////
      //iel++;
    }   /* end loop over electron */
    return esum;
  }

  NonLocalECPComponent::RealType 
  NonLocalECPComponent::evaluate(ParticleSet& W, TrialWaveFunction& psi,int iat, vector<NonLocalData>& Txy) {
    RealType esum=0.0;

    //int iel=0;
    for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++){

      register RealType r(myTable->r(nn));
      if(r>Rmax) continue;

      register RealType rinv(myTable->rinv(nn));
      register PosType  dr(myTable->dr(nn));

      int txyCounter=Txy.size();
      // Compute ratio of wave functions
      for (int j=0; j < nknot ; j++){ 
        PosType deltar(r*rrotsgrid_m[j]-dr);
        PosType newpos(W.makeMoveOnSphere(iel,deltar)); 
        //PosType newpos(W.makeMove(iel,deltar)); 
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j];
        W.rejectMove(iel);
        //psi.rejectMove(iel);
        //first, add a new NonLocalData with ratio
        Txy.push_back(NonLocalData(iel,psiratio[j],deltar));
      }
      // Compute radial potential
      for(int ip=0;ip< nchannel; ip++){
        vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];
      }

      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
        RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
        // Forming the Legendre polynomials
        lpol[0]=1.0;
        RealType lpolprev=0.0;
        for (int l=0 ; l< lmax ; l++){
          //Not a big difference
          //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
          //lpol[l+1]/=(l+1);
          lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
          lpol[l+1]*=Lfactor2[l]; 
          lpolprev=lpol[l];
        }

        //for(int l=0; l <nchannel; l++,jl++) Amat[jl]=lpol[ angpp_m[l] ]; 
        RealType lsum=0;
        for(int l=0; l <nchannel; l++) lsum += vrad[l]*lpol[ angpp_m[l] ]; 
        esum += Txy[txyCounter++].Weight *= lsum;
      } 
     //BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
     //esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
    }   /* end loop over electron */
    return esum;

  }

  ///Randomly rotate sgrid_m
  void NonLocalECPComponent::randomize_grid(ParticleSet::ParticlePos_t& sphere, bool randomize)
  {
    if(randomize) {
      //const RealType twopi(6.28318530718);
      //RealType phi(twopi*Random()),psi(twopi*Random()),cth(Random()-0.5),
      RealType phi(TWOPI*((*myRNG)())), psi(TWOPI*((*myRNG)())), cth(((*myRNG)())-0.5);
      RealType sph(std::sin(phi)),cph(std::cos(phi)),
      sth(std::sqrt(1.0-cth*cth)),sps(std::sin(psi)),
      cps(std::cos(psi));
      TensorType rmat( cph*cth*cps-sph*sps, sph*cth*cps+cph*sps,-sth*cps,
          -cph*cth*sps-sph*cps,-sph*cth*sps+cph*cps, sth*sps,
          cph*sth,             sph*sth,             cth     );
      SpherGridType::iterator it(sgridxyz_m.begin());
      SpherGridType::iterator it_end(sgridxyz_m.end());
      SpherGridType::iterator jt(rrotsgrid_m.begin());
      int ic=0;
      while(it != it_end) {*jt = dot(rmat,*it); ++it; ++jt;}
      //copy the radomized grid to sphere
      std::copy(rrotsgrid_m.begin(), rrotsgrid_m.end(), sphere.begin());
    } else {
      //copy sphere to the radomized grid
      std::copy(sphere.begin(), sphere.end(), rrotsgrid_m.begin());
    }
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
