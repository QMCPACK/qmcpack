#ifndef ATOMIC_ORBITAL_H
#define ATOMIC_ORBITAL_H

#include <vector>
#include <einspline/bspline.h>
#include <einspline/multi_bspline.h>
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Lattice/CrystalLattice.h"
#include "Configuration.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus {

  template<typename StorageType>
  class AtomicOrbital
  {
  private:
    typedef QMCTraits::PosType                    PosType;
    typedef CrystalLattice<RealType,OHMMS_DIM>    UnitCellType;
    typedef Vector<double>                        RealValueVector_t;
    typedef Vector<TinyVector<double,OHMMS_DIM> > RealGradVector_t;
    typedef Vector<complex<double> >              ComplexValueVector_t;
    typedef Vector<TinyVector<complex<double>,OHMMS_DIM> > ComplexGradVector_t;
    // Store in order 
    // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
    vector<complex<double> > YlmVec, dYlmVec, ulmVec, dulmVec, d2ulmVec;
    inline void CalcYlm(PosType rhat);
    multi_UBspline_1d_z *RadialSpline;
    // The first index is n in r^n, the second is lm = l*(l+1)+m
    Array<complex<double>,3> PolyCoefs;
    NewTimer YlmTimer, SplineTimer, SumTimer;
    RealType rmagLast;
  public:
    PosType Pos;
    RealType CutoffRadius, SplineRadius, PolyRadius;
    int SplinePoints;
    int PolyOrder;
    int lMax, Numlm, NumBands;
    UnitCellType Lattice;

    inline void set_pos  (PosType pos){ Pos = pos;   } 
    inline void set_lmax (int lmax)   { lMax = lmax; }
    inline void set_cutoff (RealType cutoff) 
    { CutoffRadius = cutoff; }
    inline void set_spline (RealType radius, int points)
    { SplineRadius = radius;  SplinePoints = points; }
    inline void set_polynomial (RealType radius, int order)
    { PolyRadius = radius; PolyOrder = order;        }
    inline void set_num_bands (int num_bands) 
    { NumBands = num_bands;                          }

    inline void registerTimers()
    {
      YlmTimer.reset();
      SplineTimer.reset();
      TimerManager.addTimer (&YlmTimer);
      TimerManager.addTimer (&SplineTimer);
      TimerManager.addTimer (&SumTimer);
    }

    void allocate()
    {
      Numlm = (lMax+1)*(lMax+1);
      YlmVec.resize(Numlm);  dYlmVec.resize(Numlm);
      ulmVec.resize  (Numlm*NumBands);  
      dulmVec.resize (Numlm*NumBands); 
      d2ulmVec.resize(Numlm*NumBands);
      cerr << "NumBands = " << NumBands << endl;
      PolyCoefs.resize(PolyOrder+1, NumBands, Numlm);
      BCtype_z bc;
      bc.lCode = NATURAL;  bc.rCode = NATURAL;
      Ugrid grid;
      grid.start = 0.0;  grid.end = SplineRadius;  grid.num = SplinePoints;
//       if (RadialSpline)
//       	destroy_Bspline (RadialSpline);
      RadialSpline = create_multi_UBspline_1d_z (grid, bc, Numlm*NumBands);
    }

    void SetBand (int band, Array<complex<double>,2> &spline_data,
		  Array<complex<double>,2> &poly_coefs)
    {
      vector<complex<double> > one_spline(SplinePoints);
      for (int lm=0; lm<Numlm; lm++) {
	int index = band*Numlm + lm;
	for (int i=0; i<SplinePoints; i++)
	  one_spline[i] = spline_data(i, lm);
	set_multi_UBspline_1d_z (RadialSpline, index, &one_spline[0]);
	for (int n=0; n<=PolyOrder; n++)
	  PolyCoefs(n,band,lm) = poly_coefs (n,lm);
      }
    }
    
    inline bool evaluate (PosType r, ComplexValueVector_t &vals);
    inline bool evaluate (PosType r, ComplexValueVector_t &val,
			  ComplexGradVector_t &grad,
			  ComplexValueVector_t &lapl);
    inline bool evaluate (PosType r, RealValueVector_t &vals);
    inline bool evaluate (PosType r, RealValueVector_t &val,
			  RealGradVector_t &grad,
			  RealValueVector_t &lapl);

    
    AtomicOrbital() : RadialSpline(NULL), 
		      YlmTimer("AtomicOrbital::CalcYlm"),
		      SplineTimer("AtomicOrbital::1D spline"),
		      SumTimer("AtomicOrbital::Summation"),
		      rmagLast(1.0e50)
    {
      // Nothing else for now
    }
  };

  
  template<typename StorageType> bool
  AtomicOrbital<StorageType>::evaluate (PosType r, ComplexValueVector_t &vals)
  {
    PosType dr = r - Pos;
    PosType u = Lattice.toUnit(dr);
    for (int i=0; i<OHMMS_DIM; i++)
      u[i] -= round(u[i]);
    dr = Lattice.toCart(u);
    double r2 = dot(dr,dr);
    if (r2 > CutoffRadius*CutoffRadius)
      return false;
    
    double rmag = std::sqrt(r2);
    PosType rhat = (1.0/rmag)*dr;
    
    // Evaluate Ylm's
    CalcYlm (rhat);

    if (std::fabs(rmag - rmagLast) > 1.0e-6) {
      // Evaluate radial functions
      if (rmag > PolyRadius)
	eval_multi_UBspline_1d_z (RadialSpline, rmag, &(ulmVec[0])); 
      else {
	double r2n = 1.0;
	for (int n=0; n <= PolyOrder; n++) {
	  int index = 0;
	  for (int i=0; i<vals.size(); i++)
	    for (int lm=0; lm<Numlm; lm++)
	      ulmVec[index++] += r2n*PolyCoefs(n,i,lm);
	  r2n *= rmag;
	}
      }
      rmagLast = rmag;
    }
    SumTimer.start();
    int index = 0;
    for (int i=0; i<vals.size(); i++) {
      vals[i] = complex<double>();
      for (int lm=0; lm < Numlm; lm++)
	vals[i] += ulmVec[index++] * YlmVec[lm];
    }
    SumTimer.stop();
    return true;
  }


  template<typename StorageType> bool
  AtomicOrbital<StorageType>::evaluate (PosType r, RealValueVector_t &vals)
  {
    // FILE *fout = fopen ("radial_splines.dat", "w");
    // for (double rmag=0.0; rmag<1.3; rmag+=0.001) {
    //   eval_multi_UBspline_1d_z (RadialSpline, rmag, &(ulmVec[0]));
    //   for (int i=0; i<Numlm*NumBands; i++)
    // 	fprintf (fout, "%1.8e ", ulmVec[i].real());
    //   fprintf (fout, "\n");
    // }
    // fclose(fout);
    // abort();

    PosType dr = r - Pos;
    PosType u = Lattice.toUnit(dr);
    for (int i=0; i<OHMMS_DIM; i++)
      u[i] -= round(u[i]);
    dr = Lattice.toCart(u);
    double r2 = dot(dr,dr);
    if (r2 > CutoffRadius*CutoffRadius)
      return false;
    
    double rmag = std::sqrt(r2);
    PosType rhat = (1.0/rmag)*dr;
    
    // Evaluate Ylm's
    CalcYlm (rhat);

    if (std::fabs(rmag - rmagLast) > 1.0e-6) {
      // Evaluate radial functions
      if (rmag > PolyRadius) {
	SplineTimer.start();
	eval_multi_UBspline_1d_z (RadialSpline, rmag, &(ulmVec[0])); 
	SplineTimer.stop();
      }
      else {
	for (int index=0; index<ulmVec.size(); index++)
	  ulmVec[index] = complex<double>();
	double r2n = 1.0;
	for (int n=0; n <= PolyOrder; n++) {
	  int index = 0;
	  for (int i=0; i<vals.size(); i++)
	    for (int lm=0; lm<Numlm; lm++)
	      ulmVec[index++] += r2n*PolyCoefs(n,i,lm);
	  r2n *= rmag;
	}
      }
      rmagLast = rmag;
    }
    SumTimer.start();
    int index = 0;
    for (int i=0; i<vals.size(); i++) {
      vals[i] = 0.0;
      for (int lm=0; lm < Numlm; lm++, index++)
    	//vals[i] += real(ulmVec[index++] * YlmVec[lm]);
    	vals[i] += (ulmVec[index].real() * YlmVec[lm].real() -
		    ulmVec[index].imag() * YlmVec[lm].imag());
    }
    SumTimer.stop();
    return true;
  }


  template<typename StorageType> bool
  AtomicOrbital<StorageType>::evaluate (PosType r, RealValueVector_t &vals,
					RealGradVector_t &grads,
					RealValueVector_t &lapl)
  {
    PosType dr = r - Pos;
    PosType u = Lattice.toUnit(dr);
    for (int i=0; i<OHMMS_DIM; i++)
      u[i] -= round(u[i]);
    dr = Lattice.toCart(u);
    double r2 = dot(dr,dr);
    if (r2 > CutoffRadius*CutoffRadius)
      return false;

    double rmag = std::sqrt(r2);
    double rInv = 1.0/rmag;
    PosType rhat = rInv*dr;
    double costheta = rhat[2];
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double cosphi = rhat[0]/sintheta;
    double sinphi = rhat[1]/sintheta;
    PosType thetahat = PosType(costheta*cosphi,
		       costheta*sinphi,
		       -sintheta);
    PosType phihat   = PosType(-sinphi, cosphi, 0.0 );

    
    // Evaluate Ylm's
    CalcYlm (rhat);

    // Evaluate radial functions
    if (rmag > PolyRadius) {
      SplineTimer.start();
      eval_multi_UBspline_1d_z_vgl 
	(RadialSpline, rmag, &(ulmVec[0]), &(dulmVec[0]), &(d2ulmVec[0])); 
      SplineTimer.stop();
    }
    else {
      for (int index=0; index<ulmVec.size(); index++) {
      	ulmVec[index]  = complex<double>();
	dulmVec[index]  = complex<double>();
	d2ulmVec[index] = complex<double>();
      }

      double r2n = 1.0, r2nm1=0.0, r2nm2=0.0;
      double dn = 0.0;
      double dnm1 = -1.0;
      for (int n=0; n <= PolyOrder; n++) {
      	int index = 0;
      	for (int i=0; i<vals.size(); i++)
      	  for (int lm=0; lm<Numlm; lm++,index++) {
	    complex<double> c = PolyCoefs(n,i,lm);
      	    ulmVec[index]   += r2n*c;
	    dulmVec[index]  += dn * r2nm1 * c;
	    d2ulmVec[index] += dn*dnm1 * r2nm2 * c;
	  }
	dn += 1.0;
	dnm1 += 1.0;
	r2nm2 = r2nm1;
	r2nm1 = r2n;
      	r2n *= rmag;
      }
    }

    SumTimer.start();
    int index = 0;
    for (int i=0; i<vals.size(); i++) {
      vals[i] = 0.0;
      for (int j=0; j<OHMMS_DIM; j++) grads[i][j] = 0.0;
      lapl[i] = 0.0;
      int lm=0;
      double grad_rhat=0.0, grad_thetahat=0.0, grad_phihat=0.0;
      for (int l=0; l<= lMax; l++)
	for (int m=-l; m<=l; m++,lm++,index++) {
	  complex<double> im(0.0,(double)m);

	  vals[i]  += real(ulmVec[index] * YlmVec[lm]);
	  grad_rhat     += real(dulmVec[index] * YlmVec[lm]);
	  grad_thetahat += real(ulmVec[index] * rInv * dYlmVec[lm]);
	  grad_phihat   += real(ulmVec[index] * im *YlmVec[lm])/(rmag*sintheta);

// 	  grads[i] += real
// 	    (dulmVec[index]                *     YlmVec[lm] * rhat     +
// 	     ulmVec[index]*rInv            *    dYlmVec[lm] * thetahat +
// 	     ulmVec[index]/(rmag*sintheta) * im *YlmVec[lm] * phihat);
	  
	  lapl[i] += real(YlmVec[lm] * 
	    (-(double)(l*(l+1))*rInv*rInv * ulmVec[index]
	     + d2ulmVec[index] + 2.0*rInv *dulmVec[index]));
	}
      grads[i] = (grad_rhat * rhat + grad_thetahat * thetahat +
		  grad_phihat * phihat);
    }
    SumTimer.stop();
    rmagLast = rmag;
    return true;
  }



  template<typename StorageType> bool
  AtomicOrbital<StorageType>::evaluate (PosType r, ComplexValueVector_t &vals,
					ComplexGradVector_t &grads,
					ComplexValueVector_t &lapl)
  {
    PosType dr = r - Pos;
    PosType u = Lattice.toUnit(dr);
    for (int i=0; i<OHMMS_DIM; i++)
      u[i] -= round(u[i]);
    dr = Lattice.toCart(u);
    double r2 = dot(dr,dr);
    if (r2 > CutoffRadius*CutoffRadius)
      return false;

    double rmag = std::sqrt(r2);
    double rInv = 1.0/rmag;
    PosType rhat = rInv*dr;
    double costheta = rhat[2];
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double cosphi = rhat[0]/sintheta;
    double sinphi = rhat[1]/sintheta;
    PosType thetahat = PosType(costheta*cosphi,
		       costheta*sinphi,
		       -sintheta);
    PosType phihat   = PosType(-sinphi, cosphi, 0.0 );

    
    // Evaluate Ylm's
    CalcYlm (rhat);

    // Evaluate radial functions
    if (rmag > PolyRadius) {
      SplineTimer.start();
      eval_multi_UBspline_1d_z_vgl 
	(RadialSpline, rmag, &(ulmVec[0]), &(dulmVec[0]), &(d2ulmVec[0])); 
      SplineTimer.stop();
    }
    else {
      for (int index=0; index<ulmVec.size(); index++) {
      	ulmVec[index]  = complex<double>();
	dulmVec[index]  = complex<double>();
	d2ulmVec[index] = complex<double>();
      }

      double r2n = 1.0, r2nm1=0.0, r2nm2=0.0;
      double dn = 0.0;
      double dnm1 = -1.0;
      for (int n=0; n <= PolyOrder; n++) {
      	int index = 0;
      	for (int i=0; i<vals.size(); i++)
      	  for (int lm=0; lm<Numlm; lm++,index++) {
	    complex<double> c = PolyCoefs(n,i,lm);
      	    ulmVec[index]   += r2n*c;
	    dulmVec[index]  += dn * r2nm1 * c;
	    d2ulmVec[index] += dn*dnm1 * r2nm2 * c;
	  }
	dn += 1.0;
	dnm1 += 1.0;
	r2nm2 = r2nm1;
	r2nm1 = r2n;
      	r2n *= rmag;
      }
    }

    SumTimer.start();
    int index = 0;
    for (int i=0; i<vals.size(); i++) {
      vals[i] = 0.0;
      for (int j=0; j<OHMMS_DIM; j++) grads[i][j] = 0.0;
      lapl[i] = 0.0;
      int lm=0;
      complex<double> grad_rhat, grad_thetahat, grad_phihat;
      for (int l=0; l<= lMax; l++)
	for (int m=-l; m<=l; m++,lm++,index++) {
	  complex<double> im(0.0,(double)m);

	  vals[i]  += ulmVec[index] * YlmVec[lm];
	  grad_rhat     += dulmVec[index] * YlmVec[lm];
	  grad_thetahat += ulmVec[index] * rInv * dYlmVec[lm];
	  grad_phihat   += (ulmVec[index] * im*YlmVec[lm])/(rmag*sintheta);
	  lapl[i] += YlmVec[lm] * 
	    (-(double)(l*(l+1))*rInv*rInv * ulmVec[index]
	     + d2ulmVec[index] + 2.0*rInv *dulmVec[index]);
	}
      for (int j=0; j<OHMMS_DIM; j++)
	grads[i][j] = grad_rhat*rhat[j] + grad_thetahat*thetahat[j] 
	  + grad_phihat * phihat[j];
    }
    SumTimer.stop();
    rmagLast = rmag;
    return true;
  }





  // Fast implementation
  // See Geophys. J. Int. (1998) 135,pp.307-309
  template<typename StorageType> void
  AtomicOrbital<StorageType>::CalcYlm (PosType rhat)
  {
    YlmTimer.start();
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
    double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
    for (int l=0; l<=lMax; l++) {
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
    YlmTimer.stop();
  }



}
#endif
