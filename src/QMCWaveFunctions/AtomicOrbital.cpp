//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "AtomicOrbital.h"

namespace qmcplusplus
{

template<> void
AtomicOrbital<std::complex<double> >::allocate()
{
  Numlm = (lMax+1)*(lMax+1);
  YlmVec.resize(Numlm);
  dYlm_dthetaVec.resize(Numlm);
  dYlm_dphiVec.resize(Numlm);
  ulmVec.resize  (Numlm*NumBands);
  dulmVec.resize (Numlm*NumBands);
  d2ulmVec.resize(Numlm*NumBands);
  PolyCoefs.resize(PolyOrder+1, NumBands, Numlm);
  BCtype_z bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  Ugrid grid;
  grid.start = 0.0;
  grid.end = SplineRadius;
  grid.num = SplinePoints;
  // if (RadialSpline) destroy_Bspline (RadialSpline);
  RadialSpline = create_multi_UBspline_1d_z (grid, bc, Numlm*NumBands);
  TwistAngles.resize(NumBands);
}

template<> void
AtomicOrbital<double>::allocate()
{
  Numlm = (lMax+1)*(lMax+1);
  YlmVec.resize(Numlm);
  dYlm_dthetaVec.resize(Numlm);
  dYlm_dphiVec.resize(Numlm);
  ulmVec.resize  (Numlm*NumBands);
  dulmVec.resize (Numlm*NumBands);
  d2ulmVec.resize(Numlm*NumBands);
  PolyCoefs.resize(PolyOrder+1, NumBands, Numlm);
  BCtype_d bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  Ugrid grid;
  grid.start = 0.0;
  grid.end = SplineRadius;
  grid.num = SplinePoints;
  RadialSpline = create_multi_UBspline_1d_d (grid, bc, Numlm*NumBands);
  TwistAngles.resize(NumBands);
}



template<> void
AtomicOrbital<std::complex<double> >::set_band(int band, Array<std::complex<double>,2> &spline_data,
    Array<std::complex<double>,2> &poly_coefs,
    PosType twist)
{
  std::vector<std::complex<double> > one_spline(SplinePoints);
  for (int lm=0; lm<Numlm; lm++)
  {
    int index = band*Numlm + lm;
    for (int i=0; i<SplinePoints; i++)
      one_spline[i] = spline_data(i, lm);
    set_multi_UBspline_1d_z (RadialSpline, index, &one_spline[0]);
    for (int n=0; n<=PolyOrder; n++)
      PolyCoefs(n,band,lm) = poly_coefs (n,lm);
  }
  TwistAngles[band] = twist;
}


// Here, we convert the complex Ylm representation to the real Ylm representation
template<> void
AtomicOrbital<double>::set_band (int band, Array<std::complex<double>,2> &spline_data,
                                 Array<std::complex<double>,2> &poly_coefs,
                                 PosType twist)
{
  std::vector<double> one_spline(SplinePoints);
  for (int l=0; l<=lMax; l++)
  {
    // Set spline for m=0
    for (int i=0; i<SplinePoints; i++)
      one_spline[i] = spline_data(i, l*(l+1)).real();
    int index = band*Numlm + l*(l+1);
    set_multi_UBspline_1d_d (RadialSpline, index, &one_spline[0]);
    // Set poly ofr m=0
    for (int n=0; n<=PolyOrder; n++)
      PolyCoefs(n,band,l*(l+1)) = poly_coefs (n,l*(l+1)).real();
    // Set spline and poly for |m| > 0
    double minus_1_to_m = -1.0;
    for (int m=1; m<=l; m++)
    {
      int lmp = l*(l+1) + m;
      int lmm = l*(l+1) - m;
      index = band*Numlm + lmp;
      for (int i=0; i<SplinePoints; i++)
        one_spline[i] = (spline_data(i, lmp).real() +
                         minus_1_to_m * spline_data(i, lmm).real());
      set_multi_UBspline_1d_d (RadialSpline, index, &one_spline[0]);
      index = band*Numlm + lmm;
      for (int i=0; i<SplinePoints; i++)
        one_spline[i] = (-spline_data(i, lmp).imag() +
                         minus_1_to_m * spline_data(i, lmm).imag());
      set_multi_UBspline_1d_d (RadialSpline, index, &one_spline[0]);
      for (int n=0; n<=PolyOrder; n++)
      {
        PolyCoefs(n,band,lmp) = (poly_coefs (n,lmp).real() +
                                 minus_1_to_m * poly_coefs(n,lmm).real());
        PolyCoefs(n,band,lmm) = (-poly_coefs (n,lmp).imag() +
                                 minus_1_to_m * poly_coefs(n,lmm).imag());
      }
      minus_1_to_m *= -1.0;
    }
  }
  TwistAngles[band] = twist;
  // AtomicOrbital<std::complex<double> > zorb;
  // zorb.set_pos (Pos);
  // zorb.set_lmax(lMax);
  // zorb.set_cutoff(CutoffRadius);
  // zorb.set_spline(SplineRadius, SplinePoints);
  // zorb.set_polynomial (PolyRadius, PolyOrder);
  // zorb.set_num_bands(NumBands);
  // zorb.allocate();
  // zorb.set_band(band, spline_data, poly_coefs, twist);
  // PosType dir(0.324, -0.8, 1.3);
  // dir = (1.0/std::sqrt(dot(dir,dir)))*dir;
  // std::ostringstream fname;
  // fname << "TestAtomic_" << band << ".dat";
  // FILE *fout = fopen (fname.str().c_str(), "w");
  // Vector<double> zval(NumBands), val(NumBands);
  // Vector<double> zlapl(NumBands), lapl(NumBands);
  // Vector<PosType> zgrad(NumBands), grad(NumBands);
  // for (double u=-1.00001; u<=1.0; u+= 0.001) {
  //   PosType r = u*CutoffRadius * dir + Pos;
  //   zorb.evaluate(r, zval, zgrad, zlapl);
  //   evaluate(r, val, grad, lapl);
  //   fprintf (fout, "%12.8f %12.8f %12.8f  %14.8e %14.8e\n",
  // 	       r[0], r[1], r[2], lapl[band], lapl[band]);
  // }
  // fclose (fout);
}
}
