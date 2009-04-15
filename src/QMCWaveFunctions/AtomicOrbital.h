#ifndef ATOMIC_ORBITAL_H
#define ATOMIC_ORBITAL_H

#include <vector>
#include <einspline/bspline.h>

class ComplexAtomicOrbital
{
private:
  Vec3 Pos;
  double CutoffRadius, SplineRadius, PolyRadius;
  int SplinePoints;
  int PolyOrder;
  int lMax, Numlm;
  // Store in order 
  // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
  vector<complex<double> > YlmVec, dYlmVec;
  void CalcYlm(Vec3 rhat);
  multi_UBspline_1d_z *radial_spline;
  // The first index is n in r^n, the second is lm = l*(l+1)+m
  Array<complex<double>,2> PolyCoefs;
public:
  void Setup (Vec3 pos, int lmax, double cutoff_radius, 
	      double spline_radius, int spline_points, 
	      double poly_radius, int poly_order,
	      int numBands)
  {
    Pos = pos;
    lMax = lmax;
    Numlm = (lmax+1)*(lmax+1);
    CutoffRadius = cutoff_radius;
    SplineRadius = spline_radius; SplinePoints = spline_points;
    PolyRadius = poly_radius;     PolyOrder  = poly_order;
    
    YlmVec.resize(numlm);
    dYlmVec.resize(numlm);
    PolyCoefs.resize(PolyOrder+1, Numlm);

    BCtype_z bc;
    bc.lCode = FLAT;  bc.rCode = NATURAL;
    Ugrid grid;
    grid.start = 0.0;  grid.end = SplineRadius;  grid.num = SplinePoints;
    if (!Spline)
      Spline = create_multi_UBspline_1d_z (grid, bc, Numlm*NumBands);
  }

  void SetBand (int Band, Array<complex<double>,2> &spline_coefs,
		Array<complex<double>,2> &poly_coefs);

  bool eval (Vec3 r, vector<complex<double> > &vals);
  bool eval (Vec3 r, vector<complex<double> > &val,
	     vector<TinyVector<complex<double>,3> > &grad,
	     vector<complex<double> > &lapl);

  AtomicOrbital() : Spline(NULL)
  {
    // Nothing else for now
  }
};

#endif
