//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "Functionals.h"
#include "../config.h"
#include <cmath>
#include <iostream>
#include "Blitz.h"



/////////////////////////////////////////////////////////////////
//                    LSDA Exchange Functions                  //
/////////////////////////////////////////////////////////////////

inline double f(double zeta)
{
  const double FourThirds = 4.0/3.0;
  const double TwoToOneThird = 1.25992104989487;//pow(2.0,1.0/3.0);
  return ((pow(1.0+zeta, FourThirds) + pow(1.0-zeta, FourThirds)
	   - 2.0) / (2.0 * (TwoToOneThird-1.0)));
}

inline double df_dzeta(double zeta)
{
  double numer = (4.0/3.0)*(pow(1.0+zeta,1.0/3.0)-pow(1.0-zeta,1.0/3.0));
  double denom = 2.0*(pow(2.0,1.0/3.0)-1.0);
  return (numer/denom);
}


void ExchangePotential (double nup, double ndown,
			double &Vup, double &Vdown)
{
  double n = nup+ndown;
  double zeta = (nup-ndown)/n;
  
  const double Third = 1.0/3.0;
  const double TwoToOneThird = 1.25992104989487;//pow(2.0,Third);
  
    double ExP = -3.0*pow(3.0*n/(8.0*M_PI), Third);
  //double ExP = -3.0/2.0*pow(3.0*n/(8.0*M_PI), Third);
  double ExF = TwoToOneThird * ExP;

  const double FourThirds = 4.0/3.0;
  double f = (pow(1.0+zeta, FourThirds) + pow(1.0-zeta, FourThirds)
	      - 2.0) / (2.0 * (TwoToOneThird-1.0));

  double dExP_dn = -pow(3.0/(8.0*M_PI*n*n), Third);
  double dExF_dn = TwoToOneThird * dExP_dn;

  double dEx_dn = dExP_dn + (dExF_dn - dExP_dn) * f;
  double df_dzeta =
    4.0*Third*(pow(1.0+zeta,Third)-pow(1.0-zeta,Third))/
    (2.0*(TwoToOneThird-1.0));

  double Ex = ExP + (ExF - ExP)*f;
  
  //fprintf (stderr, "C++ f = %1.12f\n", f);
  //fprintf (stderr, "C++ zeta = %1.12f\n", zeta);
  //fprintf (stderr, "C++ Ex = %1.12f\n", Ex);

  
  Vup   = Ex + n*dEx_dn - (zeta-1.0)*(ExF - ExP)*df_dzeta;
  Vdown = Ex + n*dEx_dn - (zeta+1.0)*(ExF - ExP)*df_dzeta;
  // Now the original VWN papers were in Rydbergs.  Multiply
  // everything by 0.5 to get in Hartrees
  Vup *= 0.5;
  Vdown *= 0.5;
  Ex *= 0.5;
  if (isnan(Vup))
    Vup = 0.0;
  if (isnan(Vdown))
    Vdown = 0.0;
  if (isnan(Ex))
    Ex = 0.0;
}
  

inline double F(double n, double A, double x0, double b, double c)
{
  const double sixth = 1.0/6.0;
  double x = pow(3.0/(4.0*M_PI*n),sixth);
  double X = x*(x+b) + c;
  double X0 = x0*(x0+b) + c;
  double Q = sqrt(4.0*c-b*b);
  double atanQ = atan(Q/(2.0*x+b));

  double term1 = log(x*x/X);
  double term2 = (2.0*b/Q)*atanQ;
  double term3 = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  double term4 = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  return (A*(term1+term2+term3+term4));
}


double dFterms(double n, double A, double x0, double b, double c)
{
  double eps = 1.0e-6;
  double np = n+eps;
  double nm = n-eps;
  const double sixth = 1.0/6.0;
  double x = pow(3.0/(4.0*M_PI*np),sixth);
  double X = x*(x+b) + c;
  double X0 = x0*(x0+b) + c;
  double Q = sqrt(4.0*c-b*b);
  double atanQ = atan(Q/(2.0*x+b));

  double term1p = log(x*x/X);
  double term2p = (2.0*b/Q)*atanQ;
  double term3p = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  double term4p = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  x = pow(3.0/(4.0*M_PI*nm),sixth);
  X = x*(x+b) + c;
  X0 = x0*(x0+b) + c;
  Q = sqrt(4.0*c-b*b);
  atanQ = atan(Q/(2.0*x+b));

  double term1m = log(x*x/X);
  double term2m = (2.0*b/Q)*atanQ;
  double term3m = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  double term4m = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  //fprintf (stderr, "dterm1 = %1.12f\n", A*(term1p-term1m)/(2.0*eps));
  //fprintf (stderr, "dterm2 = %1.12f\n", A*(term2p-term2m)/(2.0*eps));
  //fprintf (stderr, "dterm3 = %1.12f\n", A*(term3p-term3m)/(2.0*eps));
  //fprintf (stderr, "dterm4 = %1.12f\n", A*(term4p-term4m)/(2.0*eps));

  return (A*(term1p+term2p+term3p+term4p));
}



double dF_dn(double n, double A, double x0, double b, double c)
{
  const double sixth = 1.0/6.0;
  double x = pow(3.0/(4.0*M_PI*n),sixth);
  double X = x*(x+b) + c;
  double X0 = x0*(x0+b) + c;
  double Q = sqrt(4.0*c-b*b);
 
  //fprintf (stderr, "C++ x = %1.12f\n", x);
  //fprintf (stderr, "Q = %1.12f\n", Q);
 
  double n3 = n*n*n;
  double n7 = n3*n3*n;
  double prefactor = -(A/6.0) * pow(3.0/(4.0*M_PI*n7),sixth);
  double bp2x = 2.0*x + b;
  double term1 = 2.0/x - (bp2x)/X;
  double Q2m_bp2x_2_inv = 1.0/(Q*Q + (bp2x*bp2x));
  double term2 = -4.0 * b * Q2m_bp2x_2_inv;
  double term34pre = -b*x0/X0;
  double term3 = 2.0/(x-x0) - bp2x/X;
  double term4 = -4.0*(b+2.0*x0)*Q2m_bp2x_2_inv;

  //dFterms (n, A, x0, b, c);
  //fprintf (stderr, "term1 = %1.12f\n", prefactor*term1);
  //fprintf (stderr, "term2 = %1.12f\n", prefactor*term2);
  //fprintf (stderr, "term3 = %1.12f\n", term34pre*prefactor*term3);
  //fprintf (stderr, "term4 = %1.12f\n", term34pre*prefactor*term4);

  return (prefactor*(term1+term2+term34pre*(term3+term4)));
}


double dF_dn_FD (double n, double A, double x0, double b, double c)
{
  const double eps = 1.0e-6;
  double Fp = F(n+eps, A, x0, b, c);
  double Fm = F(n-eps, A, x0, b, c);
  return ((Fp-Fm)/(2.0*eps));
}



double Ec(double nup, double ndown)
{
  double n = nup + ndown;
  double zeta = (nup-ndown)/n;
  
  double EcP =    F(n, 0.0310907, -0.10498, 3.72744, 12.9352);
  //fprintf (stderr, "EcP = %1.12f\n", EcP);
  double EcF =    F(n, 0.01554535, -0.325, 7.06042, 18.0578);
  //fprintf (stderr, "EcF = %1.12f\n", EcF);
  double alphac = F(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		    1.13107, 13.0045);
  //fprintf (stderr, "alphac = %1.12f\n", alphac);

  double f_zeta = f(zeta);
  double f_doubleprime = 4.0/(9*(pow(2,1.0/3.0)-1.0));
  double beta = f_doubleprime*(EcF-EcP)/alphac -1.0;
  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double deltaEc = (alphac*f_zeta/f_doubleprime) * (1.0 + beta * zeta4);
  
  return (EcP + deltaEc);
}


void CheckCorrelationPotential (double  nup, double ndown)
{
  double eps = 1.0e-4;
  
  double n = nup + ndown;
  double zeta = (nup - ndown) / n;

  double np = n+eps;
  double nm = n-eps;
  double zetap = zeta + eps;
  double zetam = zeta - eps;

  double nupp   = 0.5*(np + np * zeta);
  double ndownp = 0.5*(np - np*zeta);
  double nupm   = 0.5*(nm + nm*zeta);
  double ndownm = 0.5*(nm - nm*zeta); 

  double Ecplus =  Ec(nupp, ndownp);
  double Ecminus = Ec(nupm, ndownm);

  double dEc_dn = (Ecplus - Ecminus)/(2.0*eps);

  double dEcP_dn = dF_dn(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  double ddeltaEc_dn = dEc_dn - dEcP_dn;
  
  //fprintf (stderr, "FD: ddeltaEc_dn = %1.12f\n", ddeltaEc_dn);




  double EcPp =    F(np, 0.0310907, -0.10498, 3.72744, 12.9352);
  double EcFp =    F(np, 0.01554535, -0.325, 7.06042, 18.0578);
  double alphacp = F(np, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		     1.13107, 13.0045);

  double EcPm =    F(nm, 0.0310907, -0.10498, 3.72744, 12.9352);
  double EcFm =    F(nm, 0.01554535, -0.325, 7.06042, 18.0578);
  double alphacm = F(nm, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		     1.13107, 13.0045);

  double f_doubleprime = 4.0/(9*(pow(2.0,1.0/3.0)-1.0));

  double betap = f_doubleprime*(EcFp-EcPp)/alphacp -1.0;
  double betam = f_doubleprime*(EcFm-EcPm)/alphacm -1.0;

  double zeta2 = zeta*zeta;
  double zeta4 = zeta2*zeta2;
  double deltaEcp = alphacp*f(zeta)/f_doubleprime * (1.0+betap*zeta4);
  double deltaEcm = alphacm*f(zeta)/f_doubleprime * (1.0+betam*zeta4);

  fprintf (stderr, "FD2: ddeltaEc_dn = %1.12f\n",
	   (deltaEcp-deltaEcm)/(2.0*eps)); 

  fprintf (stderr, "FD: dbeta_dn = %1.12f\n", (betap-betam)/(2.0*eps));

}



void CorrelationPotential(double  nup, double ndown,
			  double &Vup, double &Vdown)
{
  double EC = Ec(nup, ndown);
  //fprintf (stderr, "C++ EC = %1.12f\n", EC);
  double n = nup + ndown;
  double zeta = (nup - ndown)/n;
  double zeta2 = zeta*zeta;
  double zeta3 = zeta * zeta2;
  double zeta4 = zeta2*zeta2;

  double EcP =    F(n, 0.0310907, -0.10498, 3.72744, 12.9352);
  double EcF =    F(n, 0.01554535, -0.325, 7.06042, 18.0578);
  double alphac = F(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		    1.13107, 13.0045);
  double dEcP_dn = dF_dn(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  //double dEcP_dn = dF_dn_FD(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  double dEcF_dn = dF_dn(n, 0.01554535, -0.325, 7.06042, 18.0578);
  double dalphac_dn = dF_dn(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
			    1.13107, 13.0045);  

  double f_zeta = f(zeta);
  double f_prime = df_dzeta(zeta);
  //fprintf (stderr, "f_prime = %1.12f\n", f_prime);
  double f_doubleprime = 4.0/(9.0*(pow(2.0,1.0/3.0)-1.0));
  double beta = f_doubleprime*(EcF - EcP)/alphac -1.0;
  double dbeta_dn = f_doubleprime * ((dEcF_dn -dEcP_dn)/alphac -
				     (EcF-EcP)*dalphac_dn/(alphac*alphac));

  //fprintf (stderr, "dbeta_dn     = %1.12f\n", dbeta_dn);
  double ddeltaEc_dn = 
    dalphac_dn * f_zeta/f_doubleprime *  (1.0+beta*zeta4) +
    alphac * (f_zeta/f_doubleprime) * dbeta_dn * zeta4;

  //fprintf (stderr, "ddeltaEc_dn =     %1.12f\n", ddeltaEc_dn);
  //CheckCorrelationPotential(nup, ndown);

  double ddeltaEc_dzeta =
    alphac/f_doubleprime*((1.0+beta*zeta4)*f_prime + 4.0*beta*f_zeta*zeta3);
  
  //fprintf (stderr, "ddeltaEc_dzeta = %1.12f\n", ddeltaEc_dzeta);

  Vup   = EC + n * (dEcP_dn + ddeltaEc_dn) - (zeta-1.0) * ddeltaEc_dzeta;
  Vdown = EC + n * (dEcP_dn + ddeltaEc_dn) - (zeta+1.0) * ddeltaEc_dzeta;

  if (isnan(Vup))
    Vup = 0.0;
  if (isnan(Vdown))
    Vdown = 0.0;
}


void CPPExCorr(double nup, double ndown,
	       double &Vup, double &Vdown) {
  
  double corrpu, corrpd;
  double expu, expd;
  CorrelationPotential(nup, ndown, corrpu, corrpd);
  ExchangePotential(nup, ndown, expu, expd);
  Vup = corrpu + expu;
  Vdown = corrpd + expd;
}


/*
#define F77_EXCCOR F77_FUNC(exccor,EXCCOR)


extern "C" void 
F77_EXCCOR(double &n, double &zeta, double &exc, double &vxc, 
	   double &vpol, int &type, int &Macdonald_Vosko);


void FortranExCorr(double  nup, double  ndown,
		   double &Vup, double &Vdown)
{
  double n = nup + ndown;
  double zeta = (nup-ndown)/n;
  
  int type = 4;
  int Macdonald_Vosko=1;//0
  
  double exc, vxc, vpol;
  F77_EXCCOR(n, zeta, exc, vxc, vpol, type, Macdonald_Vosko);

  //  fprintf (stderr, "Fortran Exc = %12.8f\n", exc);

  Vup = vxc + vpol;
  Vdown = vxc - vpol;
}
  
void
FortranExCorr (double n, double &Exc, double &Vxc)
{
  double zeta = 0.0;
  //  int type = 4;
  int type = 3;  // Perdue Zunger
  int Macdonald_Vosko = 0;
  double vpol;
  F77_EXCCOR(n, zeta, Exc, Vxc, vpol, type, Macdonald_Vosko);
}

double FortranXCE (double nup, double ndown)
{
  double n = nup + ndown;
  double zeta = (nup-ndown)/n;
  
  int type = 4;
  int Macdonald_Vosko=1;//0
  
  double exc, vxc, vpol;
  F77_EXCCOR(n, zeta, exc, vxc, vpol, type, Macdonald_Vosko);
  return (exc);
}
  */
  



//  double
//  LDA::ExchangePot(double r, const Array<RadialWF,1> &WFs)
//  {

//    return (0.0);
//  }


//  double 
//  LDA::CorrelationPot(double r, const Array<RadialWF,1> &WFs)
//  {

//    return (0.0);

//  }
