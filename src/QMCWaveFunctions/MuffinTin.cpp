#include "MuffinTin.h"
#include "Ylm.h"
#include <einspline/multi_bspline.h>

// Cheesy, slow implementation for now
void
MuffinTinClass::evalYlm (TinyVector<double,3> rhat)
{
  int i=0;
  for (int l=0; l<=lMax; l++) 
    for (int m=-l; m<=l; m++, i++) 
      YlmVec[i] = Ylm(l,m,rhat);
}

// Fast implementation
// See Geophys. J. Int. (1998) 135,pp.307-309
void
MuffinTinClass::evalYlmFast (TinyVector<double,3> rhat)
{
  const double fourPiInv = 0.0795774715459477;

  double costheta = rhat[2];
  double sintheta = sqrt(1.0-costheta*costheta);
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
      double tmp = sqrt((dl+dm)*(dl-dm+1.0));
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
    double norm = sqrt((2.0*dl+1.0)*fourPiInv / sum);
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
    
  
void
MuffinTinClass::initAPW (double radius, int num_rad_points,
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
MuffinTinClass::setAPW (int orbNum, TinyVector<double,3> k,
			Array<complex<double>,2> &u_lm)
{
  kPoints[orbNum] = k;
  int numYlm = (lMax+1)*(lMax+1);
  int num_r = u_lm.extent(1);

  assert (numYlm == u_lm.extent(0));

  Array<complex<double>,1> uvec (num_r);
  for (int lm=0; lm<numYlm; lm++) {
    int spline_num = orbNum*numYlm + lm;
    for (int ir=0; ir<num_r; ir++)
      uvec(ir) = u_lm(lm,ir);
    
    set_multi_UBspline_1d_z (RadialSplines, spline_num, uvec.data());
  }
}


void 
MuffinTinClass::setLattice (Mat3 lattice)
{
  Lattice.SetDirect(lattice);
}

void
MuffinTinClass::setCenter (TinyVector<double,3> r)
{
  Center = r;
}



void
MuffinTinClass::evaluate (TinyVector<double,3> r, 
			  vector<complex<double> > &phi)
{
  TinyVector<double,3> disp, dr, L;
  disp = r - Center;
  dr = Lattice.MinImage(disp);
  L = disp - dr;
  if (dot(dr,dr) > APWRadius*APWRadius) {
    for (int i=0; i<phi.size(); i++)
      phi[i] = complex<double>();
    return;
  }
  
  double drmag = sqrt (dot(dr,dr));
  TinyVector<double,3> drhat = (1.0/drmag)*dr;

  // Evaluate the Ylms
  //evalYlm (drhat);
  evalYlmFast(drhat);

  // Evaluate the splines
  eval_multi_UBspline_1d_z (RadialSplines, drmag, RadialVec.data());

  int numYlm = (lMax+1)*(lMax+1);
  // Compute phi
  int i=0; 
  for (int iorb=0; iorb<NumOrbitals; iorb++) {
    phi[iorb] = complex<double>();
    for (int lm=0; lm<numYlm; lm++, i++)
      phi[iorb] += RadialVec[i] * YlmVec[lm];
    // Multiply by phase factor for k-point translation
    double phase = dot(L,kPoints[iorb]);
    double s,c;
    sincos (phase, &s, &c);
    phi[iorb] *= complex<double>(c,s);
  }
}


void 
MuffinTinClass::evaluateFD (TinyVector<double,3> r,
			    vector<complex<double> > &phi,
			    vector<TinyVector<complex<double>,3> > &grad,
			    vector<complex<double> > &lapl)
{
  double eps = 1.0e-6;
  TinyVector<double,3> dx(eps, 0.0, 0.0);
  TinyVector<double,3> dy(0.0, eps, 0.0);
  TinyVector<double,3> dz(0.0, 0.0, eps);

  int n = phi.size();
  vector<complex<double> > xplus(n), xminus(n), 
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
			  vector<complex<double> > &phi,
			  vector<TinyVector<complex<double>,3> > &grad,
			  vector<complex<double> > &lapl)
{
  TinyVector<double,3> disp, dr, L;
  disp = r - Center;
  dr = Lattice.MinImage(disp);
  L = disp - dr;
 
  TinyVector<double,3> rhat, thetahat, phihat;
  double drmag = sqrt (dot(dr,dr));
  rhat = (1.0/drmag)*dr;
  double costheta = rhat[2];
  double sintheta = sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  thetahat = TinyVector<double,3>(costheta*cosphi,
				  costheta*sinphi,
				  -sintheta);
  phihat   = TinyVector<double,3>(-sinphi,
				   cosphi,
				      0.0 );

  // Evaluate the Ylms
  //evalYlm (drhat);
  evalYlmFast(rhat);

  // Evaluate the splines
  eval_multi_UBspline_1d_z_vgh (RadialSplines, drmag, 
				RadialVec.data(), dRadialVec.data(), d2RadialVec.data());

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
	grad[iorb] += (dRadialVec[i]                 *     YlmVec[lm] * rhat     +
		       RadialVec[i]/drmag            *    dYlmVec[lm] * thetahat +
 		       RadialVec[i]/(drmag*sintheta) * im *YlmVec[lm] * phihat);
	lapl[iorb] += YlmVec[lm] * 
	  (-(double)(l*(l+1))/(drmag*drmag) * RadialVec[i]
	   + d2RadialVec[i] + 2.0/drmag *dRadialVec[i]);
      }
		     
    // Multiply by phase factor for k-point translation
    double phase = dot(L,kPoints[iorb]);
    double s,c;
    sincos (phase, &s, &c);
    phi[iorb] *= complex<double>(c,s);
  }



}
