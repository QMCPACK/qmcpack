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
     psigrad.resize(n);
     psigrad_source.resize(n);
     vrad.resize(m);
     dvrad.resize(m);
     wvec.resize(m);
     Amat.resize(n*m);
     dAmat.resize(n*m);
     lpol.resize(l+1,1.0);
     dlpol.resize(l+1,0.0);
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
	//        W.makeMoveOnSphere(iel,deltar); 
        W.makeMove(iel,deltar); 
#if defined(QMC_COMPLEX)
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j]*std::cos(psi.getPhaseDiff());
#else
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j];
#endif
        W.rejectMove(iel);
        psi.resetPhaseDiff();
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
  NonLocalECPComponent::evaluate(ParticleSet& W, int iat, 
				 TrialWaveFunction& psi,
				 PosType &force_iat) {
#if defined(QMC_COMPLEX)
    return 0.0;
#else
    RealType esum=0.0;
    force_iat = PosType();
    for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++){

      register RealType r(myTable->r(nn));
      if(r>Rmax) continue;

      register RealType rinv(myTable->rinv(nn));
      register PosType  dr(myTable->dr(nn));

      // Compute ratio of wave functions
      // psi.evalGrad(W,iat);
      for (int j=0; j < nknot ; j++){ 
        PosType deltar(r*rrotsgrid_m[j]-dr);
        //W.makeMoveOnSphere(iel,deltar); 
        W.makeMove(iel,deltar) ;
	RealType ratio1 = psiratio[j] = psi.ratio(W,iel)*sgridweight_m[j];
        RealType ratio2 = psi.ratioGrad(W,iel,psigrad[j]) * sgridweight_m[j];
 	if (std::fabs(ratio2 - ratio1) > 1.0e-8)
 	  fprintf (stderr, "ratio1 = %10.8f  ratio2 = %10.8f\n", 
		   ratio1, ratio2);
	psigrad[j] *= psiratio[j];
        W.rejectMove(iel);
        psi.rejectMove(iel);
      }
      // Compute radial potential
      for(int ip=0;ip< nchannel; ip++){
	double dummy;
        vrad[ip]=nlpp_m[ip]->splint(r,dvrad[ip],dummy)*wgt_angpp_m[ip];
	dvrad[ip] *= wgt_angpp_m[ip];
      }

      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
        RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
        // Forming the Legendre polynomials
        lpol[0]=1.0;
	dlpol[0] = 0.0;
        RealType lpolprev=0.0;
	RealType dlpolprev=0.0;
        for (int l=0 ; l< lmax ; l++){
          //Not a big difference
          //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
          //lpol[l+1]/=(l+1);
          lpol[l+1]  = Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
	  dlpol[l+1] = Lfactor1[l]*lpol[l] + Lfactor1[l]*zz*dlpol[l] - l*dlpolprev;
          lpol[l+1]  *= Lfactor2[l]; 
	  dlpol[l+1] *= Lfactor2[l];
          lpolprev  =  lpol[l];
	  dlpolprev = dlpol[l];
        }
        for(int l=0; l <nchannel; l++,jl++){
	  Amat[jl]  =  lpol[ angpp_m[l] ]; 
	  dAmat[jl] = dlpol[ angpp_m[l] ];
	} 
      }
      
      // Force calculation
      for (int j=0,jl=0; j<nknot; j++) {
	for (int l=0; l<nchannel; l++,jl++) {
	  // Term 1:  from dV/dr
	  force_iat += Amat[jl] * psiratio[j] * dvrad[l] * rinv * dr;
	  // Term 2:  from P_l(zz)
	  force_iat -= dAmat[jl] * psiratio[j] * vrad[l] *
	    (-rinv * rrotsgrid_m[j] + dot(dr,rrotsgrid_m[j])*rinv*rinv*rinv * dr);
	  // Term 3:  from grad psi
	  force_iat -= Amat[jl] * vrad[l] *
	    (-dot(psigrad[j],rrotsgrid_m[j])*rinv*dr + psigrad[j]);
	}
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
#endif
  }

  NonLocalECPComponent::RealType 
  NonLocalECPComponent::evaluate(ParticleSet& W, ParticleSet &ions, int iat, 
				 TrialWaveFunction& psi,
				 PosType &force_iat, PosType &pulay_iat) 
  {
#if defined(QMC_COMPLEX)
    return 0.0;
#else
    RealType esum=0.0;
    force_iat = PosType();
    pulay_iat = PosType();
    PosType psi_alpha = psi.evalGradSource(W, ions, iat);
    for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++){

      register RealType r(myTable->r(nn));
      if(r>Rmax) continue;

      register RealType rinv(myTable->rinv(nn));
      register PosType  dr(myTable->dr(nn));

      // Compute ratio of wave functions
      // psi.evalGrad(W,iat);
      for (int j=0; j < nknot ; j++){ 
        PosType deltar(r*rrotsgrid_m[j]-dr);
        //W.makeMoveOnSphere(iel,deltar); 
        W.makeMove(iel,deltar) ;
	RealType ratio1 = psiratio[j] = psi.ratio(W,iel)*sgridweight_m[j];
        RealType ratio2 = psi.ratioGrad(W,iel,psigrad[j]) * sgridweight_m[j];

	W.acceptMove (iel);
	psi.acceptMove(W,iel);
 	if (std::fabs(ratio2 - ratio1) > 1.0e-8)
 	  fprintf (stderr, "ratio1 = %10.8f  ratio2 = %10.8f\n", 
		   ratio1, ratio2);
	psigrad[j] *= psiratio[j];
	psigrad_source[j] = psiratio[j] * (psi.evalGradSource(W, ions, iat) - psi_alpha);
	//psigrad_source[j] = (psi.evalGradSource(W, ions, iat) - psiratio[j] * psi_alpha);
					   
        W.makeMove(iel,-1.0*deltar) ;
	psi.ratio(W,iel);
	GradType junk;
	psi.ratioGrad(W,iel,junk);


	W.acceptMove (iel);
	psi.acceptMove(W,iel);


        // W.rejectMove(iel);
        // psi.rejectMove(iel);
      }
      // Compute radial potential
      for(int ip=0;ip< nchannel; ip++){
	double dummy;
        vrad[ip]=nlpp_m[ip]->splint(r,dvrad[ip],dummy)*wgt_angpp_m[ip];
	dvrad[ip] *= wgt_angpp_m[ip];
      }

      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
        RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
        // Forming the Legendre polynomials
        lpol[0]=1.0;
	dlpol[0] = 0.0;
        RealType lpolprev=0.0;
	RealType dlpolprev=0.0;
        for (int l=0 ; l< lmax ; l++){
          //Not a big difference
          //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
          //lpol[l+1]/=(l+1);
          lpol[l+1]  = Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
	  dlpol[l+1] = Lfactor1[l]*lpol[l] + Lfactor1[l]*zz*dlpol[l] - l*dlpolprev;
          lpol[l+1]  *= Lfactor2[l]; 
	  dlpol[l+1] *= Lfactor2[l];
          lpolprev  =  lpol[l];
	  dlpolprev = dlpol[l];
        }
        for(int l=0; l <nchannel; l++,jl++){
	  Amat[jl]  =  lpol[ angpp_m[l] ]; 
	  dAmat[jl] = dlpol[ angpp_m[l] ];
	} 
      }
      
      // Force calculation
      for (int j=0,jl=0; j<nknot; j++) {
	for (int l=0; l<nchannel; l++,jl++) {
	  // Term 1:  from dV/dr
	  force_iat += Amat[jl] * psiratio[j] * dvrad[l] * rinv * dr;
	  // Term 2:  from P_l(zz)
	  force_iat -= dAmat[jl] * psiratio[j] * vrad[l] *
	    (-rinv * rrotsgrid_m[j] + dot(dr,rrotsgrid_m[j])*rinv*rinv*rinv * dr);
	  // Term 3:  from grad psi
	  force_iat -= Amat[jl] * vrad[l] *
	    (-dot(psigrad[j],rrotsgrid_m[j])*rinv*dr + psigrad[j]);
	  pulay_iat -= Amat[jl] * psigrad_source[j] * vrad[l];
	}
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
#endif
  }



//   NonLocalECPComponent::RealType 
//   NonLocalECPComponent::evaluate(ParticleSet& W, int iat, TrialWaveFunction& psi,
// 				 PosType &force_iat, PosType &pulay) {
//     RealType esum=0.0;
//     force_iat = PosType();
//     pulay     = PosType();
//     int NumElecs = W.G.size();
//     if (WarpNorm.size() < NumElecs) {
//       WarpNorm.resize(NumElecs);
//       Gnew.resize(nknot, NumElecs);
//       dG.resize(NumElecs);
//       dL.resize(NumElecs);
//       Gion.resize(nknot);
//     }

//     // Compute warp norm for each electron
//     for (int kel=0; kel<NumElecs; kel++)
//       WarpNorm[kel] = 0.0;
//     for (int ion=0; ion<myTable->centers(); ion++) 
//       for (int mm=myTable->M[ion],kel=0; mm<myTable->M[ion+1]; mm++,kel++)
// 	WarpNorm[kel] += WarpFunction(myTable->r(mm));
//     for (int kel=0; kel<NumElecs; kel++)
//       WarpNorm[kel] = 1.0/WarpNorm[kel];

//     for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++){

//       register RealType r(myTable->r(nn));
//       if(r>Rmax) continue;

//       register RealType rinv(myTable->rinv(nn));
//       register PosType  dr(myTable->dr(nn));

//       // Compute ratio of wave functions
//       // psi.evalGrad(W,iat);
//       for (int j=0; j < nknot ; j++){ 
//         PosType deltar(r*rrotsgrid_m[j]-dr);
//         //W.makeMoveOnSphere(iel,deltar); 
//         W.makeMove(iel,deltar) ;
//         //psi.ratioGrad(W,iel,psigrad[j]) * sgridweight_m[j];
// 	RealType ratio1 = psi.ratio(W,iel)                       * sgridweight_m[j];
// 	RealType ratio2 = psiratio[j] = psi.ratio (W,iel,dG, dL) * sgridweight_m[j];

//  	if (std::fabs(ratio2 - ratio1) > 1.0e-8)
//  	  fprintf (stderr, "ratio1 = %10.8f  ratio2 = %10.8f\n", 
// 		   ratio1, ratio2);

// 	for (int i=0; i < NumElecs; i++)
// 	  Gnew(j,i) = psiratio[j] * (W.G[i] + dG[i]);
// 	psigrad[j] = Gnew(j,iel);
	
// 	// Now compute gradient of psi w.r.t. ion iat
// 	Gion[j] = PosType();
// 	for (int mm=myTable->M[iat],kel=0; mm<myTable->M[iat+1]; mm++,kel++)
// 	  // Gion[j] -= WarpNorm[kel] * WarpFunction(myTable->r(mm)) *
// 	  //   Gnew(j,kel);
// 	  Gion[j] -= WarpNorm[kel] * WarpFunction(myTable->r(mm)) * dG[kel] * psiratio[j];

// 	Gion[j] = -1.0*dG[iel] * psiratio[j];

//         W.rejectMove(iel);
//         psi.rejectMove(iel);
//       }
//       // Compute radial potential
//       for(int ip=0;ip< nchannel; ip++){
// 	double dummy;
//         vrad[ip]=nlpp_m[ip]->splint(r,dvrad[ip],dummy)*wgt_angpp_m[ip];
// 	dvrad[ip] *= wgt_angpp_m[ip];
//       }

//       // Compute spherical harmonics on grid
//       for (int j=0, jl=0; j<nknot ; j++){ 
//         RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
//         // Forming the Legendre polynomials
//         lpol[0]=1.0;
// 	dlpol[0] = 0.0;
//         RealType lpolprev=0.0;
// 	RealType dlpolprev=0.0;
//         for (int l=0 ; l< lmax ; l++){
//           //Not a big difference
//           //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
//           //lpol[l+1]/=(l+1);
//           lpol[l+1]  = Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
// 	  dlpol[l+1] = Lfactor1[l]*lpol[l] + Lfactor1[l]*zz*dlpol[l] - l*dlpolprev;
//           lpol[l+1]  *= Lfactor2[l]; 
// 	  dlpol[l+1] *= Lfactor2[l];
//           lpolprev  =  lpol[l];
// 	  dlpolprev = dlpol[l];
//         }
//         for(int l=0; l <nchannel; l++,jl++){
// 	  Amat[jl]  =  lpol[ angpp_m[l] ]; 
// 	  dAmat[jl] = dlpol[ angpp_m[l] ];
// 	} 
//       }
      
//       // Force calculation
//       for (int j=0,jl=0; j<nknot; j++) {
// 	for (int l=0; l<nchannel; l++,jl++) {
// 	  // Term 1:  from dV/dr
// 	  force_iat += Amat[jl] * psiratio[j] * dvrad[l] * rinv * dr;
// 	  // Term 2:  from P_l(zz)
// 	  force_iat -= dAmat[jl] * psiratio[j] * vrad[l] *
// 	    (-rinv * rrotsgrid_m[j] + dot(dr,rrotsgrid_m[j])*rinv*rinv*rinv * dr);
// 	  // Term 3:  from grad psi
// 	  force_iat -= Amat[jl] * vrad[l] *
// 	    (-dot(psigrad[j],rrotsgrid_m[j])*rinv*dr + psigrad[j]);

// 	  // Pulay force, term 1
// 	  // HACK HACK HACK
// 	  pulay -= Amat[jl] * vrad[l] * Gion[j];
// 	}
//       }


//       if(nchannel==1) {
//         esum += vrad[0]*BLAS::dot(nknot, &Amat[0],&psiratio[0]);
//       } else {
//         BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
//         esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
//       }
//     }   /* end loop over electron */
    
//     // Now that we have esum=(W \psi)/psi, we can compute the second
//     // part of the Pulay correction
//     // Compute grad psi_T w.r.t. ion iat

//     // Compute warp norm for each electron
//     for (int kel=0; kel<NumElecs; kel++)
//       WarpNorm[kel] = 0.0;
    
//     for (int ion=0; ion<myTable->centers(); ion++) 
//       for (int mm=myTable->M[ion],kel=0; mm<myTable->M[ion+1]; mm++,kel++)
// 	WarpNorm[kel] += WarpFunction(myTable->r(mm));
//     for (int kel=0; kel<NumElecs; kel++)
//       WarpNorm[kel] = 1.0/WarpNorm[kel];
    
//     // Now, compute gradient
//     PosType Giat = PosType();
//     for (int mm=myTable->M[iat],kel=0; mm<myTable->M[iat+1]; mm++,kel++)
//       Giat -= WarpNorm[kel] * WarpFunction(myTable->r(mm)) *W.G[kel];
//     //pulay += esum*Giat;
    


//     return esum;
//   }




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
        //W.makeMoveOnSphere(iel,deltar); 
        PosType newpos(W.makeMove(iel,deltar));
#if defined(QMC_COMPLEX)
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j]*std::cos(psi.getPhaseDiff());
#else
        psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j];
#endif
        W.rejectMove(iel);
        psi.resetPhaseDiff();
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

  NonLocalECPComponent::RealType 
  NonLocalECPComponent::evaluate(ParticleSet& W, TrialWaveFunction& psi,int iat, 
				 vector<NonLocalData>& Txy, PosType &force_iat) 
  {
#if defined(QMC_COMPLEX)
    return 0.0;
#else
    RealType esum=0.0;
    force_iat = PosType();

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
        //W.makeMoveOnSphere(iel,deltar); 
        PosType newpos(W.makeMove(iel,deltar)); 
        //psiratio[j]=psi.ratio(W,iel)*sgridweight_m[j];
	psiratio[j]=psi.ratioGrad(W,iel,psigrad[j]) * sgridweight_m[j];
	psigrad[j] *= psiratio[j];
        W.rejectMove(iel);
        psi.rejectMove(iel);
        //first, add a new NonLocalData with ratio
        Txy.push_back(NonLocalData(iel,psiratio[j],deltar));
      }
      // Compute radial potential
      for(int ip=0;ip< nchannel; ip++){
	double dummy;
        vrad[ip]=nlpp_m[ip]->splint(r,dvrad[ip],dummy)*wgt_angpp_m[ip];
	dvrad[ip] *= wgt_angpp_m[ip];
      }

      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
        RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
	// Forming the Legendre polynomials
        lpol[0]=1.0;
	dlpol[0] = 0.0;
        RealType lpolprev=0.0;
	RealType dlpolprev=0.0;
        for (int l=0 ; l< lmax ; l++){
          lpol[l+1]  = Lfactor1[l]*zz*lpol[l]-l*lpolprev; 
	  dlpol[l+1] = Lfactor1[l]*lpol[l] + Lfactor1[l]*zz*dlpol[l] - l*dlpolprev;
          lpol[l+1]  *= Lfactor2[l]; 
	  dlpol[l+1] *= Lfactor2[l];
          lpolprev  =  lpol[l];
	  dlpolprev = dlpol[l];
        }
	for(int l=0; l <nchannel; l++,jl++){
	  Amat[jl]  =  lpol[ angpp_m[l] ]; 
	  dAmat[jl] = dlpol[ angpp_m[l] ];
	}

        //for(int l=0; l <nchannel; l++,jl++) Amat[jl]=lpol[ angpp_m[l] ]; 
        RealType lsum=0;
        for(int l=0; l <nchannel; l++) lsum += vrad[l]*lpol[ angpp_m[l] ]; 
        esum += Txy[txyCounter++].Weight *= lsum;
      } 

      // Force calculation
      for (int j=0,jl=0; j<nknot; j++) {
	for (int l=0; l<nchannel; l++,jl++) {
	  // Term 1:  from dV/dr
	  force_iat += Amat[jl] * psiratio[j] * dvrad[l] * rinv * dr;
	  // Term 2:  from P_l(zz)
	  force_iat -= dAmat[jl] * psiratio[j] * vrad[l] *
	    (-rinv * rrotsgrid_m[j] + dot(dr,rrotsgrid_m[j])*rinv*rinv*rinv * dr);
	  // Term 3:  from grad psi
	  force_iat -= Amat[jl] * vrad[l] *
	    (-dot(psigrad[j],rrotsgrid_m[j])*rinv*dr + psigrad[j]);
	}
      }

     //BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
     //esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
    }   /* end loop over electron */
    return esum;
#endif
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

  template<typename T>
  void NonLocalECPComponent::randomize_grid(vector<T> &sphere)
  {
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
    //copy the randomized grid to sphere
    for (int i=0; i<rrotsgrid_m.size(); i++)
      for (int j=0; j<OHMMS_DIM; j++)
	sphere[OHMMS_DIM*i+j] = rrotsgrid_m[i][j];
  }

  template void NonLocalECPComponent::randomize_grid(vector<float> &sphere);
  template void NonLocalECPComponent::randomize_grid(vector<double> &sphere);


}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
