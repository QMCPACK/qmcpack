//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
// 
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_EWALDSUM_H
#define OHMMS_EWALDSUM_H

#include <math.h>
#include <string>
#include <vector>
using namespace std;

template<class PL_t, class PosArray_t, class ScalarArray_t>
class EwaldSum {

public:
  typedef typename PL_t::Scalar_t Scalar_t;
  typedef typename PL_t::SingleParticlePos_t SingleParticlePos_t;
  typedef typename PL_t::Tensor_t Tensor_t;

  Scalar_t Eps;   //!< Eps = alpha^2
  Scalar_t ConvE, ConvF;//!< conversion factors for the energy and force
  Scalar_t Tr_over_Tg; // !< ratio to determine an optimal eps
  EwaldSum() { 
    Eps = -1;
    ConvE = rydberg*bohr;     // default energy conversion
    ConvF = rydberg*bohr*bohr;// default force conversion
    Tr_over_Tg = 1.0/9.0; // PIV with icc, need to be modified on different machines
    init();
  }
  ~EwaldSum() { }

  void init();//!< intializing aux values depending on Eps
  void makecells(const PL_t& lattice);
  Scalar_t energy(const PL_t& lattice, 
		  const PosArray_t& r, const ScalarArray_t& q);

  Scalar_t forcestress(const PL_t& lattice, const PosArray_t& r, const ScalarArray_t& q,
		       PosArray_t& f, Tensor_t& stress);

  void print(ostream& os) {
    os << "Eps                 = " << Eps << endl;
    os << "Reciprocal     sum  = " << gamrec << endl;
    os << "Direct         sum  = " << gamdir << endl;
    os << "Point self-senergy  = " << s3 << endl;
    os << "Charged system term = " << s4 << endl;
    os << "Maximum K components= " << maxg1 << " " << maxg2 << " " << maxg3 << endl;
    os << "Maximum R components= " << maxx1 << " " << maxx2 << " " << maxx3 << endl;

  }

private:

  Scalar_t acclog;//!< tolerance
  Scalar_t g2max; //!< |G_max|^2 for a given Eps
  Scalar_t x2max; //!< |R_n|^2   for a given Eps
  Scalar_t sqeps; //!< alpha = sqrt(eps) 
  Scalar_t gamrec;//!< sum over reciprocal terms
  Scalar_t gamdir;//!< sum over direct terms
  Scalar_t s3, s4;//!< point self-energy; charged system term
  int maxg1, maxg2, maxg3,  maxx1,  maxx2, maxx3;
  PosArray_t R;//!< cartensian vectors within a supercell
};


template<class PL_t, class PosArray_t, class ScalarArray_t>
void
EwaldSum<PL_t, PosArray_t, ScalarArray_t>::init(){

  const Scalar_t accur =1.0e-10;
  acclog = log(accur);
  g2max = 4.0*Eps*fabs(acclog);
  x2max = fabs(acclog)/Eps;                                                  
  sqeps = sqrt(Eps);                                                       

  R.InUnit = false;
}

//use wrapAroundBox to convert any configuration to cartesian in a supercell
#include "ParticleBase/ParticleUtility.h"

template<class PL_t, class PosArray_t, class ScalarArray_t>
typename EwaldSum<PL_t, PosArray_t, ScalarArray_t>::Scalar_t
EwaldSum<PL_t,PosArray_t,ScalarArray_t>::energy(const PL_t& lat, 
						const PosArray_t& rin, 
						const ScalarArray_t& q) {

  const Scalar_t pi = M_PI;
  const Scalar_t twopi = 2*M_PI;
  const Scalar_t spi = sqrt(M_PI);
  int nat = rin.size();
  if(Eps < 0) {
    Eps = pi*(pow(static_cast<Scalar_t>(nat)*0.5/(lat.Volume*lat.Volume),1.0/3.0));
    init();
  }


  maxg1 = static_cast<int>(sqrt( g2max/dot(lat.b(0),lat.b(0)))) + 1;
  maxg2 = static_cast<int>(sqrt( g2max/dot(lat.b(1),lat.b(1)))) + 1;
  maxg3 = static_cast<int>(sqrt( g2max/dot(lat.b(2),lat.b(2)))) + 1;
  maxx1 = static_cast<int>(sqrt( x2max/dot(lat.a(0),lat.a(0)))) + 1;
  maxx2 = static_cast<int>(sqrt( x2max/dot(lat.a(1),lat.a(1)))) + 1;
  maxx3 = static_cast<int>(sqrt( x2max/dot(lat.a(2),lat.a(2)))) + 1;

  R.resize(nat);//make sure the size is consistent

  //convert the current positions to Cartesian coordinates in a super cell
  wrapAroundBox(lat, rin, R);

  gamrec = 0.0;
  for(int ig1=-maxg1; ig1<= maxg1; ig1++) {    
    for(int ig2=-maxg2; ig2<=maxg2; ig2++) {
      for(int ig3=-maxg3; ig3<=maxg3; ig3++) {

	if(ig1 == 0 && ig2 == 0 && ig3 == 0) continue;	// exclude G = 0;

	// this G
	SingleParticlePos_t tau = 
	  static_cast<Scalar_t>(ig1)*lat.b(0) +
	  static_cast<Scalar_t>(ig2)*lat.b(1) + 
	  static_cast<Scalar_t>(ig3)*lat.b(2);
	tau *= twopi; // multiply 2*pi
                                                                              
	Scalar_t tau2 = dot(tau,tau);

	Scalar_t t2e  = tau2/(4.0*Eps);// |2\pi G|^2/(4*Eps)
                                                  
	if ( -t2e < acclog) continue;	//  contribution neglegible 
	Scalar_t expfac = exp( - t2e)/t2e; // 4\eps correction later

	//gamma-ewald: sum over kappa
	Scalar_t sumr = 0.0;
	Scalar_t sumi = 0.0;
	for(int iat=0; iat<nat; iat++) {
	  Scalar_t kdotr = dot(R[iat], tau);
	  sumr += q[iat]*cos(kdotr);
	  sumi += q[iat]*sin(kdotr);
	}

	Scalar_t term = expfac*(sumr*sumr + sumi*sumi); 
	gamrec += term; // unit = 1/(1/A)^2 = A^2 
      }//loop-ig3
    }//loop-ig2
  }//loop-ig1


  //summation over the direct lattice
  gamdir = 0.0;
  for(int ix1=-maxx1; ix1<=maxx1; ix1++) {
    for(int ix2=-maxx2; ix2<=maxx2; ix2++) {
      for(int ix3=-maxx3; ix3<=maxx3; ix3++) {
	SingleParticlePos_t xlp =
	  static_cast<Scalar_t>(ix1)*lat.a(0) +
	  static_cast<Scalar_t>(ix2)*lat.a(1) + 
	  static_cast<Scalar_t>(ix3)*lat.a(2);
	for(int iat=0; iat<nat; iat++) {
	  for(int jat=0; jat<nat; jat++) {
	    SingleParticlePos_t xxx = R[iat]-R[jat]-xlp;
	    Scalar_t xxx2 = dot(xxx,xxx);
	    
	    // (l-prim,kappa) = (0,kappa0) is excluded
	    if (xxx2 < 1.0e-8) continue;
	    Scalar_t arg = sqrt(xxx2) * sqeps;
                                                 
	    // neglegible contribution
	    if ( -arg < acclog) continue;

	    // h(x) = erfc(x)/x, sqrt(eps) is corrected later
	    Scalar_t hx = erfc(arg)/arg; 

	    //.....gamma-ewald: 
	    gamdir += q[iat]*q[jat]*hx;
	  }//jat
	}//iat
      }//ix3
    }//ix2
  }//ix1

  //the third and fourth sums
  s3 = 0.0;
  s4 = 0.0;
  for(int iat=0; iat<nat; iat++) {
    s3 += q[iat]*q[iat];
    s4 += q[iat];
  }

  Scalar_t vol = lat.Volume;
  /*!\latex {
    reciprocal term:
    ${2\pi\over V} {1\over 4\eps} V_G= {\pi\over 2 V \eps}V_G$
    direct term:${\sqrt{\eps}\over 2} V_R$ 
    }
  unit = (au)^2/A
  */
  //Scalar_t esum =  pi/(2.0*vol*Eps)*gamrec + sqeps/2.0*gamdir 
  //  - sqeps/spi*s3 - pi/(2.0*vol*Eps)*s4*s4;

  Scalar_t esum =  pi/(vol*Eps)*gamrec + sqeps*gamdir 
    - 2.0*sqeps/spi*s3 - pi/(vol*Eps)*s4*s4;
  return esum*ConvE;
}


template<class PL_t, class PosArray_t, class ScalarArray_t>
typename EwaldSum<PL_t, PosArray_t, ScalarArray_t>::Scalar_t
EwaldSum<PL_t,PosArray_t,ScalarArray_t>
::forcestress(const PL_t& lat, const PosArray_t& rin, const ScalarArray_t& q,
	      PosArray_t& f0,  EwaldSum<PL_t,PosArray_t,ScalarArray_t>::Tensor_t& stress) {

  const Scalar_t pi = M_PI;
  const Scalar_t twopi = 2*M_PI;
  const Scalar_t spi = sqrt(M_PI);

  int nat = rin.size();
  if(Eps < 0) {
    // Eps = \pi *(t_R/T_G N/V^2)^{1/3} 
    // t_R/t_G ~ 1/9 optimal on PIV with icc
    Eps = pi*(pow(static_cast<Scalar_t>(nat)*Tr_over_Tg/(lat.Volume*lat.Volume),1.0/3.0));
    init();
  }
  maxg1 = static_cast<int>(sqrt( g2max/dot(lat.b(0),lat.b(0)))) + 1;
  maxg2 = static_cast<int>(sqrt( g2max/dot(lat.b(1),lat.b(1)))) + 1;
  maxg3 = static_cast<int>(sqrt( g2max/dot(lat.b(2),lat.b(2)))) + 1;
  maxx1 = static_cast<int>(sqrt( x2max/dot(lat.a(0),lat.a(0)))) + 1;
  maxx2 = static_cast<int>(sqrt( x2max/dot(lat.a(1),lat.a(1)))) + 1;
  maxx3 = static_cast<int>(sqrt( x2max/dot(lat.a(2),lat.a(2)))) + 1;

  R.resize(nat);//make sure the size is consistent

  //convert the current positions to Cartesian coordinates in a super cell
  wrapAroundBox(lat, rin, R);

  Scalar_t vol = lat.Volume;

  //!< pi/(V*eps) factor for rec force
  Scalar_t rec_fac = 2.0*ConvF*pi/(vol*Eps);
  Tensor_t stress_rec, stress_dir;//stress for rec and dir terms

  gamrec = 0.0;
  for(int ig1=-maxg1; ig1<= maxg1; ig1++) {    
    for(int ig2=-maxg2; ig2<=maxg2; ig2++) {
      for(int ig3=-maxg3; ig3<=maxg3; ig3++) {
	if(ig1 == 0 && ig2 == 0 && ig3 == 0) continue;

	// this G
	SingleParticlePos_t tau = 
	  static_cast<Scalar_t>(ig1)*lat.b(0) +
	  static_cast<Scalar_t>(ig2)*lat.b(1) + 
	  static_cast<Scalar_t>(ig3)*lat.b(2);
	tau *= twopi; // multiply 2*pi
                                                                               
	Scalar_t tau2 = dot(tau,tau);   // G^2
	Scalar_t t2e  = tau2/(4.0*Eps); //G^2/4/Eps

	if ( -t2e < acclog) continue;	//contribution neglegible 
	Scalar_t expfac = exp( - t2e)/t2e;

	//gamma-ewald: sum over kappa
	Scalar_t sumr = 0.0;
	Scalar_t sumi = 0.0;
	for(int iat=0; iat<nat; iat++) {
	  Scalar_t kdotr = dot(R[iat], tau);
	  sumr += q[iat]*cos(kdotr);
	  sumi += q[iat]*sin(kdotr);
	}

	Scalar_t term   = expfac*(sumr*sumr + sumi*sumi);
        gamrec += term;                                                    

	//stress:
	Scalar_t factor = term * 2.0*(t2e + 1.0)/tau2;
	for(int idir=0; idir<3; idir++) {
	  for(int jdir=0; jdir<3; jdir++) {
	    Scalar_t sum = tau[idir]*tau[jdir]*factor;
	    if(idir == jdir) sum -= term;
	    stress_rec(idir,jdir) += sum;
	  }
	}

	//forces: loop over the atoms kappa-0   
	for(int iat=0; iat<nat; iat++) {     

	  Scalar_t sum = 0.0;

	  // summation over kappa          
	  for(int jat=0; jat<nat; jat++) {
	    SingleParticlePos_t dr =R[iat]-R[jat];
	    sum += q[jat]*sin(dot(tau,dr));
	  }
	  sum *= expfac*rec_fac; 
	  f0[iat] += sum*tau;
	}//loop-iat

      }//loop-ig3
    }//loop-ig2
  }//loop-ig1



  Scalar_t dir_fac = 2.0*ConvF*Eps;//!< Eps factor for direct force x2
  gamdir = 0.0;
  //summation over the direct lattice
  for(int ix1=-maxx1; ix1<=maxx1; ix1++) {
    for(int ix2=-maxx2; ix2<=maxx2; ix2++) {
      for(int ix3=-maxx3; ix3<=maxx3; ix3++) {
	SingleParticlePos_t xlp =
	  static_cast<Scalar_t>(ix1)*lat.a(0) +
	  static_cast<Scalar_t>(ix2)*lat.a(1) + 
	  static_cast<Scalar_t>(ix3)*lat.a(2);
	for(int iat=0; iat<nat; iat++) {
	  for(int jat=0; jat<nat; jat++) {
	    SingleParticlePos_t xxx = R[iat]-R[jat]-xlp;
	    Scalar_t xxx2 = dot(xxx,xxx);

	    // (l-prim,kappa) = (0,kappa0) is excluded
	    if (xxx2 < 1.0e-8) continue;
	    Scalar_t arg = sqrt(xxx2) * sqeps;
                                                 
	    // neglegible contribution
	    if ( -arg < acclog) continue;
	    // function h(x) = erfc(x)/x
	    Scalar_t hx = erfc(arg)/arg; 

	    //.....gamma-ewald:
	    gamdir += q[iat]*q[jat]*hx;

	    //
	    //stress: 
	    // function h'(x)*x
	    Scalar_t  arg2 = arg*arg, hprim;
	    if(-arg2>acclog)
	      hprim = - hx - 2.0/spi*exp( - arg2);
	    else 
	      hprim = - hx;
                                                            
	    Scalar_t factor = q[iat]*q[jat]*hprim/xxx2;
	    stress_dir += factor*outerProduct(xxx,xxx); 

	    //forces
	    hprim = hprim/arg;//h'(x)
	    factor = hprim*q[jat]/sqrt(xxx2)*dir_fac;

	    // Important: f0 contains conversion factor + sign
	    f0[iat] -= factor*xxx; 
	  }//jat
	}//iat

      }//ix3
    }//ix2
  }//ix1

  //the third and fourth sums
  s3 = 0.0;
  s4 = 0.0;
  for(int iat=0; iat<nat; iat++) {
    s3 += q[iat]*q[iat];
    s4 += q[iat];
  }


  rec_fac = ConvE*pi/(vol*Eps);//!< recip term including 2
  dir_fac = ConvE*sqeps;       //!< direct term including 2
  for(int idir=0; idir<OHMMS_DIM; idir++) {
    for(int jdir=0; jdir<OHMMS_DIM; jdir++) {
      stress(idir,jdir) += 
	rec_fac*stress_rec(idir,jdir) + dir_fac*stress_dir(idir,jdir);
      if(idir==jdir)
	stress(idir,jdir) += rec_fac*s4*s4;
    }
  }

  /*!\latex {
    reciprocal term:
    ${2\pi\over V} {1\over 4\eps} V_G= {\pi\over 2 V \eps}V_G$
    direct term:${\sqrt{\eps}\over 2} V_R$ 
    }
  unit = (au)^2/A
  */
  return 
    rec_fac*gamrec+dir_fac*gamdir-2.0*ConvE*sqeps/spi*s3 - rec_fac*s4*s4;
}

#endif


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
