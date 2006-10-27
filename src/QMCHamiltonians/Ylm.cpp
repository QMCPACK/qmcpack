#include "Ylm.h"
#include <iostream>

namespace qmcplusplus {
  double LegendrePll (int l, double x) {
    if (l==0)
      return 1.0;
    else {
      double sqt = std::sqrt(1.0-x)*std::sqrt(1.0+x);
      double val = 1.0;
      double dblfact = 1.0;
      for (int i=1; i<=l; i++) {
	val *= -sqt;
	val *= dblfact;
	dblfact += 2.0;
      }
      return val;
    }
  }

  double LegendrePlm (int l, int m, double x) {
    if (m < 0) {
      m = abs (m);
      double posval = LegendrePlm (l, m, x);
      double sign = (m%2==0) ? 1.0 : -1.0;
      double mfact = 1.0;
      for (int i=2; i<=(l-m); i++)
	mfact *= (double)i;
      double pfact = 1.0;
      for (int i=2; i<=(l+m); i++)
	pfact *= (double)i;
      return posval * sign*mfact/pfact;
    }
    // Now we can assume that m is not negative
    double pmm = LegendrePll (m, x);
    double pmp1m = x*(2*m+1)*pmm;

    if (m == l) 
      return pmm;
    else if (l==(m+1))
      return pmp1m;
    else { // Use recursive formula
      double Plm2m = pmm;
      double Plm1m = pmp1m;
      double Pl;
      for (int i=m+2; i<=l;  i++) {
	Pl = (1.0/(double)(i-m)) *
	  (x*(2*i-1)*Plm1m - (i+m-1)*Plm2m);
	Plm2m = Plm1m;
	Plm1m = Pl;
      }
      return Pl;
    }
  }


  complex<double> Ylm(int l, int m, TinyVector<double,3> r)
  {
    double costheta, phi;
    costheta = r[0];
    phi   = std::atan2(r[2],r[1]);
    int lmm = l - m;
    int lpm = l + m;
    double mfact = 1.0;
    double pfact = 1.0;
    for (int i=lmm; i>0; i--)
      mfact *= (double)i;
    for (int i=lpm; i>0; i--)
      pfact *= (double)i;
    double prefactor = std::sqrt ((double)(2*l+1)*mfact/(4.0*M_PI*pfact));
    double Plm = LegendrePlm (l, m, costheta);
    complex<double> e2imphi (std::cos(m*phi), std::sin(m*phi));
    return prefactor * Plm * e2imphi;
  }
}
