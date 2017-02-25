#ifndef QMCPLUSPLUS_BESSEL_H
#define QMCPLUSPLUS_BESSEL_H

namespace qmcplusplus
{

template<typename T>
void bessel_steed_array_cpu(const int lmax, const T x, T* jl_x) {
// use steed/barnett method to calculate regular array of spherical bessel
// function from 0 to l.  This follows the gsl implementation
   if (lmax < 0 || x < 0.0) {
      int j;
      for (j = 0; j<=lmax; j++) jl_x[j] = 0.0;
   } else if (x == 0.0) {
      int j;
      for (j = 1; j<=lmax; j++) jl_x[j] = 0.0;
      jl_x[0] = 1.0;
   } else if (x < 0.00024) {
      double inv_fact = 1.0;
      double x_l = 1.0;
      int l;
      for(l=0; l<=lmax; l++) {
	 jl_x[l]  = x_l * inv_fact;
	 jl_x[l] *= 1.0 - 0.5*x*x/(2.0*l+3.0);
	 inv_fact /= 2.0*l+3.0;
	 x_l *= x;
      }
   } else {
      // Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)]
      double x_inv = 1.0/x;
      double W = 2.0*x_inv;
      double F = 1.0;
      double FP = (lmax+1.0) * x_inv;
      double B = 2.0*FP + x_inv;
      double D = 1.0/B;
      double del = -D;
      
      FP += del;
      
      /* continued fraction */
      do {
	 B += W;
	 D = 1.0/(B-D);
	 del *= (B*D - 1.);
	 FP += del;
	 if(D < 0.0) F = -F;
      } while(fabs(del) >= fabs(FP) * 1.19209e-07);
      
      FP *= F;
      
      if(lmax > 0) {
	 /* downward recursion */
	 double XP2 = FP;
	 double PL = lmax * x_inv;
	 int L  = lmax;
	 int LP;
	 jl_x[lmax] = F;
	 for(LP = 1; LP<=lmax; LP++) {
	    jl_x[L-1] = PL * jl_x[L] + XP2;
	    FP = PL*jl_x[L-1] - jl_x[L];
	    XP2 = FP;
	    PL -= x_inv;
	    --L;
	 }
	 F = jl_x[0];
      }  
      
      /* normalization */
      W = x_inv / hypot(FP, F);
      jl_x[0] = W*F;
      if(lmax > 0) {
	int L;
	for(L=1; L<=lmax; L++) {
          jl_x[L] *= W;
	}
      }
   }
}

}

#endif
