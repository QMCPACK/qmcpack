#ifndef QMCPLUSPLUS_NEW_BESSEL_H
#define QMCPLUSPLUS_NEW_BESSEL_H

#include <cmath>
namespace qmcplusplus
{

/** @ingroup Numerics
 * @brief Compute spherical bessel funciton from 0 to lmax
 *
 * Using Steed/Barnett algorithm from Comp. Phys. Comm. 21, 297 (1981)
 */
template<typename T>
void bessel_steed_array_cpu(const int lmax, const T x, T* jl)
{
  if (lmax<0)
    throw std::runtime_error("Negative lmax not allowed!");
  else if (x<0)
    throw std::runtime_error("Negative x not allowed!");
  else if (x==0.0)
  {
    std::fill(jl,jl+lmax+1,T(0.0));
    jl[0] = T(1.0);
  }
  else if (x<0.00024)
  {
    const double cone(1);
    double inv_accumuated = cone;
    double x_l = cone;
    for(int l=0; l<=lmax; l++)
    {
      jl[l] = x_l*inv_accumuated;
      const double inv=cone/(2*l+3);
      jl[l] *= cone-0.5*x*x*inv;
      inv_accumuated *= inv;
      x_l *= x;
    }
  }
  else
  {
    const double cone(1);
    const double ctwo(2);
    double XI(cone/x);
    double W(XI*ctwo);
    double F(cone);
    double FP((lmax+1)*XI);
    double B(FP*ctwo+XI);
    double D(cone/B);
    double DEL(-D);
    FP += DEL;

    do
    {
      B += W;
      D = cone/(B-D);
      DEL *= (B*D-cone);
      FP += DEL;
      if (D<0.0) F = -F;
    } while (std::fabs(DEL)>=std::fabs(FP)*1.19209e-07);

    FP *= F;
    jl[0] = F;

    if (lmax>0)
    {
      double PL = lmax*XI;
      jl[lmax] = F;
      for(int L=lmax; L>=1; L--)
      {
        jl[L-1] = PL*jl[L]+FP;
        FP = PL*jl[L-1]-jl[L];
        PL -= XI;
      }
      F = jl[0];
    }

    // using hypot instead of sqrt to avoid overflow/underflow
    W = XI/std::hypot(FP,F);
    for(int L=0; L<=lmax; L++)
      jl[L] *= W;
  }
}

}

#endif
