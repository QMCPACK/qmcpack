//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Ken Esler                            //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: esler@uiuc.edu                                     //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#include "MuffinTin.h"
#include <einspline/bspline_base.h>
#include <einspline/bspline.h>
#include <einspline/multi_bspline.h>
#include <cmath>

namespace qmcplusplus {
  // Fast implementation
  // See Geophys. J. Int. (1998) 135,pp.307-309
  void
  MuffinTinClass::evalYlm (TinyVector<double,3> rhat)
  {
    const double fourPiInv = 0.0795774715459477;
    
    double costheta = rhat[2];
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double cottheta = costheta/sintheta;
    
    double cosphi, sinphi;
    cosphi=rhat[0]/sintheta;
    sinphi=rhat[1]/sintheta;
    
    complex<double> e2iphi(cosphi, sinphi);
    
    
    double lsign = 1.0;
    double dl = 0.0;
    for (int l=0; l<=lMax; l++) {
      double XlmVec[2*l+1], dXlmVec[2*l+1];
      XlmVec[2*l]  = lsign;  
      dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
      XlmVec[0]    = lsign*XlmVec[2*l];
      dXlmVec[0]   = lsign*dXlmVec[2*l];
      double dm = dl;
      double msign = lsign;
      for (int m=l; m>0; m--) {
	double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
	XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
	dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
	// Copy to negative m
	XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
	dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
	msign *= -1.0;
	dm -= 1.0;
      }
      double sum = 0.0;
      for (int m=-l; m<=l; m++) 
	sum += XlmVec[l+m]*XlmVec[l+m];
      // Now, renormalize the Ylms for this l
      double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
      for (int m=-l; m<=l; m++) {
	XlmVec[l+m]  *= norm;
	dXlmVec[l+m] *= norm;
      }
      
      // Multiply by azimuthal phase and store in YlmVec
      complex<double> e2imphi (1.0, 0.0);
      for (int m=0; m<=l; m++) {
	YlmVec[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
	YlmVec[l*(l+1)-m]  =  XlmVec[l-m]*conj(e2imphi);
	dYlmVec[l*(l+1)+m] = dXlmVec[l+m]*e2imphi;
	dYlmVec[l*(l+1)-m] = dXlmVec[l-m]*conj(e2imphi);
	e2imphi *= e2iphi;
      } 
      
      dl += 1.0;
      lsign *= -1.0;
    }
  }
    
  bool
  MuffinTinClass::inside(TinyVector<double,3> r)
  {
    TinyVector<double,3> ru(PrimLattice.toUnit(r-Center));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= round (ru[i]);
    TinyVector<double,3> dr(PrimLattice.toCart(ru));
    return dot(dr,dr) < APWRadius*APWRadius;
  }


  void
  MuffinTinClass::init_APW (double radius, int num_rad_points,
			    int lmax, int numOrbitals)
  {
    lMax = lmax;
    APWRadius = radius;
    NumOrbitals = numOrbitals;
    
    // Create the grid
    Ugrid rgrid;
    rgrid.start = 0.0;
    rgrid.end = radius;
    rgrid.num = num_rad_points;
    
    // Boundary conditions
    BCtype_z rBC;
    rBC.lCode = NATURAL;
    rBC.rCode = NATURAL;
    
    // Create the multi-spline
    int numYlm = (lmax+1)*(lmax+1);
    int numSplines = numYlm*numOrbitals;
    RadialSplines = create_multi_UBspline_1d_z (rgrid, rBC, numSplines);
    
    // Resize internal storage
    YlmVec.resize(numYlm);
    dYlmVec.resize(numYlm);
    RadialVec.resize(numSplines);
    dRadialVec.resize(numSplines);
    d2RadialVec.resize(numSplines);
    kPoints.resize(numOrbitals);
  }
  
  void
  MuffinTinClass::set_APW (int orbNum, TinyVector<double,3> k,
			   Array<complex<double>,2> &u_lm,
			   double Z)
  {
    kPoints[orbNum] = k;
    int numYlm = (lMax+1)*(lMax+1);
    int num_r = u_lm.size(1);
    
    if (numYlm != u_lm.size(0)) 
      app_error() << "Wrong dimension in MuffinTinClass::setAPW.\n";

    /////////////////////////////////////////////////////////// 
    // To get the correct behavior near the nucleus, we will //
    // actually spline u_lm(r)/r^l, and then multiply this   //
    // back on when we evaluate.                             //
    /////////////////////////////////////////////////////////// 
    double dr = APWRadius / (double)(num_r -1);
    Array<complex<double>,1> uvec (num_r);
    for (int ir=1; ir<num_r; ir++) {
      double r = (double)ir * dr;
      double r2l = 1.0;
      for (int l=0; l<=lMax; l++) {
	for (int m=-l; m<=l; m++) {
	  int lm = l*(l+1) + m;
	  u_lm(lm, ir) /= r2l;
	}
	r2l *= r;
      }
    }

    for (int l=0; l<=lMax; l++) 
      for (int m=-l; m<=l; m++) {
	int lm = l*(l+1) + m;
	u_lm(lm, 0)  = 2.0*u_lm(lm,1) - u_lm(lm,2);
      }

	
    for (int lm=0; lm<numYlm; lm++) {
      int spline_num = orbNum*numYlm + lm;
      for (int ir=0; ir<num_r; ir++)
	uvec(ir) = u_lm(lm,ir);
   
      set_multi_UBspline_1d_z (RadialSplines, spline_num,
			       uvec.data());
      BCtype_z rBC;
      rBC.rCode = NATURAL;
      rBC.lCode = DERIV1;
      if (lm == 0) {
	rBC.lVal_r = -Z*uvec(0).real();
	rBC.lVal_i = -Z*uvec(0).imag();
      }
      else {
	rBC.lVal_r = -0.5*Z*uvec(0).real();
	rBC.lVal_i = -0.5*Z*uvec(0).imag();
      }

      //if (lm != 0)
      //rBC.lCode = NATURAL;

      set_multi_UBspline_1d_z_BC (RadialSplines, spline_num, uvec.data(), rBC);
    }
  }
  
  
  void 
  MuffinTinClass::set_lattice (Tensor<RealType,3> lattice)
  {
    PrimLattice.set(lattice);
  }
  
  void
  MuffinTinClass::set_center (TinyVector<double,3> r)
  {
    Center = r;
  }
  


  void
  MuffinTinClass::evaluate (TinyVector<double,3> r, 
			    Vector<complex<double> > &phi)
  {
    TinyVector<double,3> disp, u, dr, L;
    disp = r - Center;
    TinyVector<double,3> ru(PrimLattice.toUnit(disp));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= round (ru[i]);
    dr = PrimLattice.toCart(ru);

    L = disp - dr;

    if (dot(dr,dr) > APWRadius*APWRadius) {
      for (int i=0; i<phi.size(); i++)
	phi[i] = complex<double>();
      return;
    }
    
    double drmag = std::sqrt (dot(dr,dr));
    TinyVector<double,3> drhat = (1.0/drmag)*dr;
    
    // Evaluate the Ylms
    //evalYlm (drhat);
    evalYlm(drhat);
    
    // Evaluate the splines
    eval_multi_UBspline_1d_z (RadialSplines, drmag, RadialVec.data());
    
    // Multiply by r^l term
    int j=0; 
    for (int iorb=0; iorb<NumOrbitals; iorb++) {
      double r2l = 1.0;
      for (int l=0; l<=lMax; l++) {
	for (int m=-l; m<=l; m++) {
	  int lm = l*(l+1) + m;
	  RadialVec[j] *= r2l;
	  j++;
	}
	r2l *= drmag;
      }
    }

    
    int numYlm = (lMax+1)*(lMax+1);
    // Compute phi
    int i=0; 
    for (int iorb=0; iorb<NumOrbitals; iorb++) {
      phi[iorb] = complex<double>();
      for (int lm=0; lm<numYlm; lm++, i++)
	phi[iorb] += RadialVec[i] * YlmVec[lm];
      // Multiply by phase factor for k-point translation
      double phase = -dot(L,kPoints[iorb]);
      double s,c;
      sincos (phase, &s, &c);
      phi[iorb] *= complex<double>(c,s);
    }
  }
  
  
  void 
  MuffinTinClass::evaluateFD (TinyVector<double,3> r,
			      Vector<complex<double> > &phi,
			      Vector<TinyVector<complex<double>,3> > &grad,
			      Vector<complex<double> > &lapl)
  {
    double eps = 1.0e-6;
    TinyVector<double,3> dx(eps, 0.0, 0.0);
    TinyVector<double,3> dy(0.0, eps, 0.0);
    TinyVector<double,3> dz(0.0, 0.0, eps);
    
    int n = phi.size();
    Vector<complex<double> > xplus(n), xminus(n), 
      yplus(n), yminus(n), zplus(n), zminus(n);
    
    evaluate (r, phi);
    evaluate (r+dx, xplus);    evaluate (r-dx, xminus);
    evaluate (r+dy, yplus);    evaluate (r-dy, yminus);
    evaluate (r+dz, zplus);    evaluate (r-dz, zminus);
    for (int i=0; i<n; i++) {
      grad[i][0] = (xplus[i]-xminus[i])/(2.0*eps);
      grad[i][1] = (yplus[i]-yminus[i])/(2.0*eps);
      grad[i][2] = (zplus[i]-zminus[i])/(2.0*eps);
      lapl[i]    = (xplus[i]+xminus[i]+yplus[i]+yminus[i]+zplus[i]+zminus[i]
		    - 6.0*phi[i])/(eps*eps);
    }
  }
  
  
  
  
  void 
  MuffinTinClass::evaluate (TinyVector<double,3> r,
			    Vector<complex<double> > &phi,
			    Vector<TinyVector<complex<double>,3> > &grad,
			    Vector<complex<double> > &lapl)
  {
    TinyVector<double,3> disp, dr, L;
    disp = r - Center;
    TinyVector<double,3> ru(PrimLattice.toUnit(disp));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= round (ru[i]);
    dr = PrimLattice.toCart(ru);
    L = disp - dr;
    
    if (dot(dr,dr) > APWRadius*APWRadius) {
      for (int i=0; i<phi.size(); i++) {
	phi[i] = lapl[i] = complex<double>();
	for (int j=0; j<3; j++)
	  grad[i][j] = complex<double>();
      }
      return;
    }

    TinyVector<double,3> rhat, thetahat, phihat;
    double drmag = std::sqrt (dot(dr,dr));
    rhat = (1.0/drmag)*dr;
    double costheta = rhat[2];
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double cosphi = rhat[0]/sintheta;
    double sinphi = rhat[1]/sintheta;
    thetahat = TinyVector<double,3>(costheta*cosphi,
				    costheta*sinphi,
				    -sintheta);
    phihat   = TinyVector<double,3>(-sinphi,
				    cosphi,
				    0.0 );
    
    // Evaluate the Ylms
    evalYlm(rhat);
    
    // Evaluate the splines
    eval_multi_UBspline_1d_z_vgh (RadialSplines, drmag, 
				  RadialVec.data(), 
				  dRadialVec.data(), 
				  d2RadialVec.data());

    // Multiply by r^l term
    int j=0; 
    for (int iorb=0; iorb<NumOrbitals; iorb++) {
      double r2l = 1.0;
      double r2lm1 = 1.0/drmag;
      double r2lm2 = 1.0/(drmag*drmag);
      for (int l=0; l<=lMax; l++) {
	for (int m=-l; m<=l; m++) {
	  int lm = l*(l+1) + m;
	  complex<double> u = RadialVec[j];
	  complex<double> du = dRadialVec[j];
	  complex<double> d2u = d2RadialVec[j];

	  RadialVec[j] = r2l * u;
	  dRadialVec[j] = (double)l * r2lm1 * u +
	    r2l * du;
	  d2RadialVec[j] = (double)(l*(l-1)) * r2lm2 * u
	    + 2.0 * (double)l*r2lm1 * du + r2l * d2u;
	  
	  j++;
	}
	r2l   *= drmag;
	r2lm1 *= drmag;
	r2lm2 *= drmag;
      }
    }
    
    int numYlm = (lMax+1)*(lMax+1);
    // Compute phi
    int i=0; 
    for (int iorb=0; iorb<NumOrbitals; iorb++) {
      phi[iorb] = complex<double>();
      grad[iorb][0] = grad[iorb][1] = grad[iorb][2] = complex<double>();
      lapl[iorb] = complex<double>();
      int lm=0;
      for (int l=0; l<=lMax; l++)
	for (int m=-l; m<=l; m++, lm++,i++) {
	  complex<double> im(0.0,(double)m);
	  phi[iorb]  += RadialVec[i] * YlmVec[lm];
	  grad[iorb] += 
	    (dRadialVec[i]                 *     YlmVec[lm] * rhat     +
	     RadialVec[i]/drmag            *    dYlmVec[lm] * thetahat +
	     RadialVec[i]/(drmag*sintheta) * im *YlmVec[lm] * phihat);
	  lapl[iorb] += YlmVec[lm] * 
	    (-(double)(l*(l+1))/(drmag*drmag) * RadialVec[i]
	     + d2RadialVec[i] + 2.0/drmag *dRadialVec[i]);
	}
      
      // Multiply by phase factor for k-point translation
      double phase = -dot(L,kPoints[iorb]);
      double s,c;
      sincos (phase, &s, &c);
      phi[iorb]  *= complex<double>(c,s);
      grad[iorb] *= complex<double>(c,s);
      lapl[iorb] *= complex<double>(c,s);
    }
  }



  void
  MuffinTinClass::addCore (int l, int m, double rmax, Vector<double> &g0,
			   TinyVector<double,3> kVec, double Z)
  {
    Ugrid rgrid;
    rgrid.start = 0.0;
    rgrid.end   = rmax;
    rgrid.num = g0.size();
    int N = g0.size();


    // Replace the first few data points with the exact hydrogenic
    // wave function
    const int matchPoint = 4;
    double dr = rmax / (N-1);
    double rMatch = dr * matchPoint;
    double cMatch = g0[matchPoint]/std::exp(-Z*rMatch);
    for (int i=0; i<matchPoint; i++) {
      double r = dr*i;
      g0[i] = cMatch * std::exp(-Z*r);
    }
    
    BCtype_d rBC;
    rBC.lCode = DERIV1;
    rBC.lVal  = -Z*g0[0];
    rBC.rCode = NATURAL;

    // Compute radius at which to truncate the core state
    double norm = 0.0;
    int i=N-1;
    while ( i>0 && norm<1.0e-5) {
      double r = dr*(double) i;
      double u = g0[i];
      norm += u*u*r*r * dr;
      i--;
    }
    double rcut = 1000.0*(i+1)*dr;

    CoreRadii.push_back(rcut);
    UBspline_1d_d *spline = create_UBspline_1d_d (rgrid, rBC, g0.data());
    double u, du, d2u;
    eval_UBspline_1d_d_vgl (spline, 0.0, &u, &du, &d2u);
    fprintf (stderr, "Set boundary value = %1.16e\n", -Z*g0[0]);
    fprintf (stderr, "Evaluated value    = %1.16e\n", du);

    CoreSplines.push_back(spline);
    Core_lm.push_back(TinyVector<int,2>(l,m));
    Core_kVecs.push_back (kVec);
    
    NumCore++;
  }

  void
  MuffinTinClass::evaluateCore (TinyVector<double,3> r, 
				Vector<complex<double> > &phi,
				int first)
  {
    TinyVector<double,3> disp, dr, drhat;
    disp = r - Center;
    TinyVector<double,3> ru(PrimLattice.toUnit(disp));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= round (ru[i]);
    dr = PrimLattice.toCart(ru);

    double drmag = std::sqrt(dot(dr,dr));
    drhat = (1.0/drmag) * dr;

    // This is a slow hack
    evalYlm (drhat);  
  
    for (int i=0; i<CoreSplines.size(); i++) { 
      if (drmag > CoreRadii[i]) 
	phi[first+i] = complex<double>();
      else {
	int l = Core_lm[i][0];
	int m = Core_lm[i][1];
	int lm = l*(l+1)+m;
	complex<double> ylm = YlmVec[lm];
	double u;
	eval_UBspline_1d_d (CoreSplines[i], drmag, &u);
	phi[first+i] = ylm*(u);
	// double phase = dot (r, Core_kVecs[i]);
	// double s, c;
	// sincos(phase, &s, &c);
	// phi[first+i] *= complex<double>(c,s);
      }
    }
  }

 void
  MuffinTinClass::evaluateCore (TinyVector<double,3> r, 
				Vector<complex<double> > &phi,
				Vector<TinyVector<complex<double>,3> > &grad,
				Vector<complex<double> > &lapl, int first)
  {
    TinyVector<double,3> disp, dr;
    disp = r - Center;
    TinyVector<double,3> ru(PrimLattice.toUnit(disp));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= round (ru[i]);
    dr = PrimLattice.toCart(ru);

    TinyVector<double,3> rhat, thetahat, phihat;
    double drmag = std::sqrt (dot(dr,dr));
    rhat = (1.0/drmag)*dr;
    double costheta = rhat[2];
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double cosphi = rhat[0]/sintheta;
    double sinphi = rhat[1]/sintheta;
    thetahat = TinyVector<double,3>(costheta*cosphi,
				    costheta*sinphi,
				    -sintheta);
    phihat   = TinyVector<double,3>(-sinphi,
				    cosphi,
				    0.0 );

    // This is a slow hack
    evalYlm (rhat);  
    for (int i=0; i<CoreSplines.size(); i++) { 
      if (drmag > CoreRadii[i]) {
	phi[first+i]   = complex<double>();
	grad[first+i]  = complex<double>() * rhat;
	lapl[first+i] =  complex<double>();
      }
      else {
	int l = Core_lm[i][0];
	int m = Core_lm[i][1];
	int lm = l*(l+1)+m;
	complex<double> ylm = YlmVec[lm];
	complex<double> im(0.0,(double)m);
	
	double u, du, d2u;
	eval_UBspline_1d_d_vgl (CoreSplines[i], drmag, &u, &du, &d2u);
	
	phi[first+i] = ylm*u;      
	grad[first+i] = (du                 *     YlmVec[lm] * rhat     +
			 u/drmag            *    dYlmVec[lm] * thetahat +
			 u/(drmag*sintheta) * im *YlmVec[lm] * phihat);
	lapl[first+i] = YlmVec[lm] * (-(double)(l*(l+1))/(drmag*drmag) * u
				      + d2u + 2.0/drmag *du );
	
      }
    }
  }
}
