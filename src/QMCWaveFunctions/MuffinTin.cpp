//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include <einspline/bspline_base.h>
#include <einspline/nubspline.h>
#include <einspline/multi_nubspline.h>
#include <Numerics/DeterminantOperators.h>
#include <Particle/DistanceTableData.h>
#include <QMCWaveFunctions/MuffinTin.h>
#include <config/stdlib/math.h>


namespace qmcplusplus
{

// M is the number of basis functions.  For each value of x, y
// should contain the values to be fitted.  F should contain
// all the basis functions evaluated at each x.
void
MuffinTinClass::LinFit (std::vector<double> &y,                  // input
                        std::vector<TinyVector<double,2> > &F,   // input
                        TinyVector<double,2> &a )           // output
{
  int M=2;
  int N = F.size();
  if (y.size() != F.size())
    app_error() << "Different number of rows of basis functions than"
                << " of data points in LinFit.  Exitting.\n";
  // Next, construct alpha matrix
  Matrix<double> alpha(M,M), alphaInv(M,M), ident(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
    {
      alpha(k,j) = 0.0;
      for (int i=0; i<N; i++)
        alpha(k,j) += F[i][j] * F[i][k];
    }
  // Next, construct beta vector
  Vector<double> beta(M);
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta[k] += y[i]*F[i][k];
  // Now, invert alpha
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      alphaInv(i,j) = alpha(i,j);
  double det = invert_matrix(alphaInv);
  for (int i=0; i<M; i++)
  {
    a[i] = 0.0;
    for (int j=0; j<M; j++)
      a[i] += alphaInv(i,j) * beta[j];
  }
}


// M is the number of basis functions.  For each value of x, y
// should contain the values to be fitted.  F should contain
// all the basis functions evaluated at each x.
void
MuffinTinClass::LinFit (std::vector<double> &y,                  // input
                        std::vector<TinyVector<double,3> > &F,   // input
                        TinyVector<double,3> &a )           // output
{
  int M=3;
  int N = F.size();
  // Next, construct alpha matrix
  Matrix<double> alpha(M,M), alphaInv(M,M), ident(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
    {
      alpha(k,j) = 0.0;
      for (int i=0; i<N; i++)
        alpha(k,j) += F[i][j] * F[i][k];
    }
  // Next, construct beta vector
  Vector<double> beta(M);
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta[k] += y[i]*F[i][k];
  // Now, invert alpha
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      alphaInv(i,j) = alpha(i,j);
  double det = invert_matrix(alphaInv);
  for (int i=0; i<M; i++)
  {
    a[i] = 0.0;
    for (int j=0; j<M; j++)
      a[i] += alphaInv(i,j) * beta[j];
  }
}


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
  std::complex<double> e2iphi(cosphi, sinphi);
  double lsign = 1.0;
  double dl = 0.0;
  for (int l=0; l<=lMax; l++)
  {
    double XlmVec[2*l+1], dXlmVec[2*l+1];
    XlmVec[2*l]  = lsign;
    dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
    XlmVec[0]    = lsign*XlmVec[2*l];
    dXlmVec[0]   = lsign*dXlmVec[2*l];
    double dm = dl;
    double msign = lsign;
    for (int m=l; m>0; m--)
    {
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
    for (int m=-l; m<=l; m++)
    {
      XlmVec[l+m]  *= norm;
      dXlmVec[l+m] *= norm;
    }
    // Multiply by azimuthal phase and store in YlmVec
    std::complex<double> e2imphi (1.0, 0.0);
    for (int m=0; m<=l; m++)
    {
      YlmVec[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
      YlmVec[l*(l+1)-m]  =  XlmVec[l-m]*qmcplusplus::conj(e2imphi);
      dYlmVec[l*(l+1)+m] = dXlmVec[l+m]*e2imphi;
      dYlmVec[l*(l+1)-m] = dXlmVec[l-m]*qmcplusplus::conj(e2imphi);
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
MuffinTinClass::inside(TinyVector<double,3> r,
                       bool &in, bool &needBlend)
{
  TinyVector<double,3> ru(PrimLattice.toUnit(r-Center));
  for (int i=0; i<OHMMS_DIM; i++)
    ru[i] -= round (ru[i]);
  TinyVector<double,3> dr(PrimLattice.toCart(ru));
  in = dot(dr,dr) < APWRadius*APWRadius;
  if (in)
    needBlend = dot(dr,dr) > BlendRadius*BlendRadius;
}


void
MuffinTinClass::blend_func (double r, double &b)
{
  if (r < BlendRadius)
    b = 0.0;
  else
  {
    double x = (r - BlendRadius)/(APWRadius - BlendRadius);
    b = 1.0 - 10.0 * x*x*x +15.0*x*x*x*x - 6.0*x*x*x*x*x;
  }
}

void
MuffinTinClass::blend_func (double r, double &b, double &db,
                            double &d2b)
{
  if (r < BlendRadius)
    b = db = d2b = 0.0;
  else
  {
    double dr = APWRadius - BlendRadius;
    double drInv = 1.0/dr;
    double x = (r - BlendRadius)*drInv;
    b =  1.0 - 10.0 * x*x*x +15.0*x*x*x*x - 6.0*x*x*x*x*x;
    db = drInv * (-30.0*x*x + 60.0*x*x*x -30.0*x*x*x*x);
    d2b = drInv * drInv * (-60.0*x + 180.0*x*x -120.0*x*x*x);
  }
}


// void
// MuffinTinClass::blend_func(double r, double &b)
// {
//   if (r < BlendRadius)
//     b = 0.0;
//   else {
//     double x = (r - BlendRadius)/(APWRadius - BlendRadius);
//     b =  0.5*(std::cos(M_PI*x)+1.0);
//   }
// }

// void
// MuffinTinClass::blend_func (double r, double &b, double &db,
// 			      double &d2b)
// {
//   if (r < BlendRadius)
//     b = db = d2b = 0.0;
//   else {
//     double dr = APWRadius - BlendRadius;
//     double drInv = 1.0/dr;
//     double x = (r - BlendRadius)*drInv;
//     b =  0.5*(std::cos(M_PI*x)+1.0);
//     db = -0.5*M_PI*std::sin(M_PI*x)*drInv;
//     d2b = -0.5*M_PI*M_PI*std::cos(M_PI*x)*drInv*drInv;
//   }
// }




TinyVector<double,3>
MuffinTinClass::disp(TinyVector<double,3> r)
{
  TinyVector<double,3> ru(PrimLattice.toUnit(r-Center));
  for (int i=0; i<OHMMS_DIM; i++)
    ru[i] -= round (ru[i]);
  return PrimLattice.toCart(ru);
}


void
MuffinTinClass::init_APW (Vector<double> rgrid,
                          int lmax, int numOrbitals)
{
  lMax = lmax;
  APWRadius = rgrid[rgrid.size()-1];
  // HACK HACK HACK
  BlendRadius = APWRadius - 0.0;
  NumOrbitals = numOrbitals;
  // Set rSmall.
  // Find first place where (r[i+1]-r[i]) > 1e-5
  int ir=0;
  while ((rgrid[ir+1]-rgrid[ir]) < drMin)
    ir++;
  iSmall = ir;
  rSmall = rgrid[ir];
  // Create the grid
  RadialGrid =  create_log_grid (rgrid[0], APWRadius, rgrid.size());
  //RadialGrid = create_general_grid (rgrid.data(), rgrid.size());
  for (int i=0; i<rgrid.size(); i++)
    if (std::abs(rgrid[i]-RadialGrid->points[i]) > 1.0e-12)
      app_error() << "Error in creating log grid.\n"
                  << "rgrid[i] = " << rgrid[i] << "   "
                  << "RadialGrid->points[i] = "
                  << RadialGrid->points[i] << std::endl;
  // Boundary conditions
  BCtype_z rBC;
  rBC.lCode = NATURAL;
  rBC.rCode = NATURAL;
  // Create the multi-spline
  int numYlm = (lmax+1)*(lmax+1);
  int numSplines = numYlm*numOrbitals;
  RadialSplines =
    create_multi_NUBspline_1d_z (RadialGrid, rBC, numSplines);
  // Resize internal storage
  YlmVec.resize(numYlm);
  dYlmVec.resize(numYlm);
  RadialVec.resize(numSplines);
  dRadialVec.resize(numSplines);
  d2RadialVec.resize(numSplines);
  Small_r_APW_Fits.resize(numSplines);
  kPoints.resize(numOrbitals);
}

void
MuffinTinClass::set_APW (int orbNum, TinyVector<double,3> k,
                         Array<std::complex<double>,2> &u_lm,
                         Array<std::complex<double>,1> &du_lm_final,
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
  Array<std::complex<double>,1> uvec (num_r);
  double rlast2l = 1.0;
  int lastr = u_lm.size(1)-1;
  double rlast = RadialGrid->points[lastr];
  for (int l=0; l<=lMax; l++)
  {
    for (int m=-l; m<=l; m++)
    {
      int lm = l*(l+1) + m;
      std::complex<double> u = u_lm(lm,lastr);
      std::complex<double> du = du_lm_final(lm);
      du_lm_final(lm) = (1.0/rlast2l) * (du - (double)l/rlast * u);
    }
    rlast2l *= rlast;
  }
  for (int ir=0; ir<num_r; ir++)
  {
    double r = RadialGrid->points[ir];
    double r2l = 1.0;
    for (int l=0; l<=lMax; l++)
    {
      for (int m=-l; m<=l; m++)
      {
        int lm = l*(l+1) + m;
        u_lm(lm, ir) = u_lm(lm,ir)/r2l; //u_lm(lm, ir) /= r2l;
      }
      r2l *= r;
    }
  }
  // Temp vectors for small r fit
  std::vector<std::complex<double> > uSmall(iSmall+1);
  std::vector<double> rSmall(iSmall+1);
  for (int l=0; l<=lMax; l++)
  {
    for (int m=-l; m<=l; m++)
    {
      int lm = l*(l+1) + m;
      int spline_num = orbNum*numYlm + lm;
      for (int ir=0; ir<num_r; ir++)
        uvec(ir) = u_lm(lm,ir);
      // Set small r coefficients
      for (int ir=0; ir<=iSmall; ir++)
      {
        uSmall[ir] = uvec(ir);
        rSmall[ir] = RadialGrid->points[ir];
      }
      Small_r_APW_Fits[spline_num].FitCusp(rSmall, uSmall, -Z/(double)(l+1));
      set_multi_NUBspline_1d_z (RadialSplines, spline_num,
                                uvec.data());
      BCtype_z rBC;
      rBC.rCode = DERIV1;
      rBC.lCode = DERIV1;
      std::complex<double> u0 = uvec(0);
      rBC.lVal_r = -Z*u0.real()/(double)(l+1);
      rBC.lVal_i = -Z*u0.imag()/(double)(l+1);
      rBC.rVal_r = du_lm_final(lm).real();
      rBC.rVal_i = du_lm_final(lm).imag();
      set_multi_NUBspline_1d_z_BC (RadialSplines, spline_num, uvec.data(), rBC);
    }
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
                          Vector<std::complex<double> > &phi)
{
  TinyVector<double,3> disp, u, dr, L;
  disp = r - Center;
  TinyVector<double,3> ru(PrimLattice.toUnit(disp));
  for (int i=0; i<OHMMS_DIM; i++)
    ru[i] -= round (ru[i]);
  dr = PrimLattice.toCart(ru);
  L = disp - dr;
  if (dot(dr,dr) > APWRadius*APWRadius)
  {
    for (int i=0; i<phi.size(); i++)
      phi[i] = std::complex<double>();
    return;
  }
  double drmag = std::sqrt (dot(dr,dr));
  TinyVector<double,3> drhat = (1.0/drmag)*dr;
  // Evaluate the Ylms
  //evalYlm (drhat);
  evalYlm(drhat);
  // Evaluate the splines
  if (drmag > rSmall)
    eval_multi_NUBspline_1d_z (RadialSplines, drmag, RadialVec.data());
  else
    for (int i=0; i<RadialVec.size(); i++)
      Small_r_APW_Fits[i].eval (drmag, RadialVec[i]);
  // Multiply by r^l term
  int j=0;
  for (int iorb=0; iorb<NumOrbitals; iorb++)
  {
    double r2l = 1.0;
    for (int l=0; l<=lMax; l++)
    {
      for (int m=-l; m<=l; m++)
      {
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
  for (int iorb=0; iorb<NumOrbitals; iorb++)
  {
    phi[iorb] = std::complex<double>();
    for (int lm=0; lm<numYlm; lm++, i++)
      phi[iorb] += RadialVec[i] * YlmVec[lm];
    // Multiply by phase factor for k-point translation
    double phase = -dot(L,kPoints[iorb]);
    double s,c;
    sincos (phase, &s, &c);
    phi[iorb] *= std::complex<double>(c,s);
  }
}


void
MuffinTinClass::evaluateFD (TinyVector<double,3> r,
                            Vector<std::complex<double> > &phi,
                            Vector<TinyVector<std::complex<double>,3> > &grad,
                            Vector<std::complex<double> > &lapl)
{
  double eps = 1.0e-6;
  TinyVector<double,3> dx(eps, 0.0, 0.0);
  TinyVector<double,3> dy(0.0, eps, 0.0);
  TinyVector<double,3> dz(0.0, 0.0, eps);
  int n = phi.size();
  Vector<std::complex<double> > xplus(n), xminus(n),
         yplus(n), yminus(n), zplus(n), zminus(n);
  evaluate (r, phi);
  evaluate (r+dx, xplus);
  evaluate (r-dx, xminus);
  evaluate (r+dy, yplus);
  evaluate (r-dy, yminus);
  evaluate (r+dz, zplus);
  evaluate (r-dz, zminus);
  for (int i=0; i<n; i++)
  {
    grad[i][0] = (xplus[i]-xminus[i])/(2.0*eps);
    grad[i][1] = (yplus[i]-yminus[i])/(2.0*eps);
    grad[i][2] = (zplus[i]-zminus[i])/(2.0*eps);
    lapl[i]    = (xplus[i]+xminus[i]+yplus[i]+yminus[i]+zplus[i]+zminus[i]
                  - 6.0*phi[i])/(eps*eps);
  }
}

void
MuffinTinClass::evaluate (TinyVector<double,3> r,
                          Vector<std::complex<double> > &phi,
                          Vector<TinyVector<std::complex<double>,3> > &grad,
                          Vector<Tensor<std::complex<double>,3> > &hess)
{
  APP_ABORT("Hessian not inplemented in MuffinTinClass::evaluate. \n");
}


void
MuffinTinClass::evaluate (TinyVector<double,3> r,
                          Vector<std::complex<double> > &phi,
                          Vector<TinyVector<std::complex<double>,3> > &grad,
                          Vector<std::complex<double> > &lapl)
{
  TinyVector<double,3> disp, dr, L;
  disp = r - Center;
  TinyVector<double,3> ru(PrimLattice.toUnit(disp));
  for (int i=0; i<OHMMS_DIM; i++)
    ru[i] -= round (ru[i]);
  dr = PrimLattice.toCart(ru);
  L = disp - dr;
  if (dot(dr,dr) > APWRadius*APWRadius)
  {
    for (int i=0; i<phi.size(); i++)
    {
      phi[i] = lapl[i] = std::complex<double>();
      for (int j=0; j<3; j++)
        grad[i][j] = std::complex<double>();
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
  if (drmag > rSmall)
    // Evaluate the splines
    eval_multi_NUBspline_1d_z_vgh (RadialSplines, drmag,
                                   RadialVec.data(),
                                   dRadialVec.data(),
                                   d2RadialVec.data());
  else
    for (int i=0; i< RadialVec.size(); i++)
      Small_r_APW_Fits[i].eval (drmag, RadialVec[i], dRadialVec[i], d2RadialVec[i]);
  // Multiply by r^l term
  int num_lm = (lMax+1)*(lMax+1);
  int j = 0;
  for (int iorb=0; iorb<NumOrbitals; iorb++)
  {
    double r2l = 1.0;
    double r2lm1 = 1.0/drmag;
    double r2lm2 = 1.0/(drmag*drmag);
    for (int l=0; l<=lMax; l++)
    {
      for (int m=-l; m<=l; m++)
      {
        int lm = l*(l+1) + m;
        std::complex<double> u = RadialVec[j];
        std::complex<double> du = dRadialVec[j];
        std::complex<double> d2u = d2RadialVec[j];
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
  int lStop = (drmag < rSmall) ? 2 : lMax;
  lStop = lMax;
  // Compute phi
  for (int iorb=0; iorb<NumOrbitals; iorb++)
  {
    int i = numYlm * iorb;
    phi[iorb] = std::complex<double>();
    grad[iorb][0] = grad[iorb][1] = grad[iorb][2] = std::complex<double>();
    lapl[iorb] = std::complex<double>();
    int lm=0;
    for (int l=0; l<=lStop; l++)
      for (int m=-l; m<=l; m++, lm++,i++)
      {
        std::complex<double> im(0.0,(double)m);
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
    phi[iorb]  *= std::complex<double>(c,s);
    grad[iorb] *= std::complex<double>(c,s);
    lapl[iorb] *= std::complex<double>(c,s);
  }
}



void
MuffinTinClass::addCore (int l, int m, Vector<double> &r, Vector<double> &g0,
                         TinyVector<double,3> kVec, double Z)
{
  int N = r.size();
  NUgrid *rgrid =
    create_log_grid (r[0], r[N-1], N);
  // NUgrid *rgrid = create_general_grid (r.data(), N);
  // Compute small-r coefficients
  int irSmall=0;
  while ((r[irSmall+1] - r[irSmall]) < drMin && irSmall<(r.size()-1))
    irSmall++;
  rSmallCore = r[irSmall+1];
  //fprintf (stderr, "rSmallCore = %1.8f  irSmall = %d\n",
  //         rSmallCore, irSmall);
  std::vector<double> vals(irSmall+50), rvals(irSmall+50);
  for (int ir=0; ir<irSmall+50; ir++)
  {
    vals[ir] = g0[ir];
    rvals[ir] = r[ir];
  }
  ExpFitClass<4> smallFit;
  smallFit.FitCusp(rvals, vals, -Z/(double)(l+1));
  Small_r_Core_Fits.push_back(smallFit);
  BCtype_d rBC;
  rBC.lCode = NATURAL;
  rBC.lVal  = -Z*g0[0];
  rBC.rCode = FLAT;
  // Compute radius at which to truncate the core state
  double norm = 0.0;
  int i=N-1;
  while ( i>1 && norm<1.0e-5)
  {
    double u = g0[i];
    double dr = r[i] - r[i-1];
    norm += u*u*r[i]*r[i]*dr;
    i--;
  }
  double rcut = r[i+1];
  CoreRadii.push_back(rcut);
  int jstart = 0;
  while (r[jstart] < 1.0)
    jstart++;
  jstart = std::min(i-30, jstart);
  // Compute large-r coefficients
  std::vector<TinyVector<double,2> > bfuncs(i+1-jstart);
  TinyVector<double,2> largeCoefs;
  vals.resize(i+1-jstart);
  for (int j=0; j<bfuncs.size(); j++)
  {
    bfuncs[j][0] = 1.0;
    bfuncs[j][1] = r[j+jstart];
    vals[j] = std::log(g0[j+jstart]);
  }
  LinFit(vals, bfuncs, largeCoefs);
  LargerCoreCoefs.push_back(largeCoefs);
  // Create nonuniform B-spline.
  NUBspline_1d_d *spline = create_NUBspline_1d_d (rgrid, rBC, g0.data());
  double u, du, d2u;
  eval_NUBspline_1d_d_vgl (spline, r[0], &u, &du, &d2u);
  CoreSplines.push_back(spline);
  Core_lm.push_back(TinyVector<int,2>(l,m));
  Core_kVecs.push_back (kVec);
  NumCore++;
}

void
MuffinTinClass::evaluateCore (TinyVector<double,3> r,
                              Vector<std::complex<double> > &phi,
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
  for (int i=0; i<CoreSplines.size(); i++)
  {
    int l = Core_lm[i][0];
    int m = Core_lm[i][1];
    int lm = l*(l+1)+m;
    std::complex<double> ylm = YlmVec[lm];
    double u;
    if (drmag < rSmallCore)
      Small_r_Core_Fits[i].eval (drmag, u);
    else
      if (drmag < CoreRadii[i])
        eval_NUBspline_1d_d (CoreSplines[i], drmag, &u);
      else
        if (drmag < 2.75)
        {
          double c0 = LargerCoreCoefs[i][0];
          double c1 = LargerCoreCoefs[i][1];
          u   = std::exp(c0 + c1*drmag);
        }
        else
          u = 0.0;
    phi[first+i] = ylm*(u);
    // double phase = dot (r, Core_kVecs[i]);
    // double s, c;
    // sincos(phase, &s, &c);
    // phi[first+i] *= std::complex<double>(c,s);
  }
}

void
MuffinTinClass::evaluateCore (TinyVector<double,3> r,
                              Vector<std::complex<double> > &phi,
                              Vector<TinyVector<std::complex<double>,3> > &grad,
                              Vector<Tensor<std::complex<double>,3> > &hess, int first)
{
}

void
MuffinTinClass::evaluateCore (TinyVector<double,3> r,
                              Vector<std::complex<double> > &phi,
                              Vector<TinyVector<std::complex<double>,3> > &grad,
                              Vector<std::complex<double> > &lapl, int first)
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
  for (int i=0; i<CoreSplines.size(); i++)
  {
    int l = Core_lm[i][0];
    int m = Core_lm[i][1];
    int lm = l*(l+1)+m;
    std::complex<double> ylm = YlmVec[lm];
    std::complex<double> im(0.0,(double)m);
    double u, du, d2u;
    if (drmag < rSmallCore)
      Small_r_Core_Fits[i].eval (drmag, u, du, d2u);
    else
      if (drmag < CoreRadii[i])
        eval_NUBspline_1d_d_vgl (CoreSplines[i], drmag, &u, &du, &d2u);
      else
        if (drmag < 2.75)
        {
          double c0 = LargerCoreCoefs[i][0];
          double c1 = LargerCoreCoefs[i][1];
          u   = std::exp(c0 + c1*drmag);
          du  = c1 * u;
          d2u = c1 * du;
        }
        else
          u = du = d2u = 0.0;
    phi[first+i] = ylm*u;
    grad[first+i] = (du                 *     YlmVec[lm] * rhat     +
                     u/drmag            *    dYlmVec[lm] * thetahat +
                     u/(drmag*sintheta) * im *YlmVec[lm] * phihat);
    lapl[first+i] = YlmVec[lm] * (-(double)(l*(l+1))/(drmag*drmag) * u
                                  + d2u + 2.0/drmag *du );
  }
}


}
