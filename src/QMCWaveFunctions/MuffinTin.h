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

#ifndef MUFFIN_TIN_H
#define MUFFIN_TIN_H

#include <vector>
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Optimize/VarList.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Lattice/CrystalLattice.h"
#include "QMCWaveFunctions/MuffinTin.h"
#include <einspline/bspline_base.h>
#include <einspline/nubspline_structs.h>
#include <einspline/multi_nubspline_structs.h>
#include "Configuration.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus {
  template<int M>
  class ExpFitClass 
  {
  private:
    TinyVector<double,M> Coefs, dCoefs, d2Coefs;
  public:
    inline void Fit (vector<double> &r, vector<double> &u);
    inline void FitCusp(vector<double> &r, vector<double> &u, double cusp);
    inline void eval (double r, double &u);
    inline void eval (double r, double &u, double &du, double &d2u);
  };

  template<int M> void 
  ExpFitClass<M>::FitCusp (vector<double> &r, vector<double> &u, double cusp)
  {
    int N = r.size();

    if (r.size() != u.size()) 
      app_error() << "Different number of rows of basis functions than"
		  << " of data points in LinFit.  Exitting.\n";
    vector<TinyVector<double,M-1> > F(N);
    vector<double> log_u(N);
    for (int i=0; i<N; i++) {
      log_u[i] = std::log (u[i]) - cusp * r[i];
      double r2jp1 = r[i]*r[i];
      F[i][0] = 1.0;
      for (int j=1; j<M-1; j++) {
	F[i][j] = r2jp1;
	r2jp1 *= r[i];
      }
    }
      
    // Next, construct alpha matrix
    Matrix<double> alpha(M-1,M-1), alphaInv(M-1,M-1);
    alpha = 0.0;
    for (int j=0; j<M-1; j++)
      for (int k=0; k<M-1; k++) {
	alpha(k,j) = 0.0;
	for (int i=0; i<N; i++)
	  alpha(k,j) += F[i][j] * F[i][k];
      }
    
    // Next, construct beta vector
    TinyVector<double,M-1> beta;
    beta = 0.0;
    for (int k=0; k<M-1; k++)
      for (int i=0; i<N; i++)
	beta[k] += log_u[i]*F[i][k];
    
    // Now, invert alpha
    for (int i=0; i<M-1; i++)
      for (int j=0; j<M-1; j++)
	alphaInv(i,j) = alpha(i,j);
    
    double det = invert_matrix(alphaInv);
	    
    TinyVector<double,M-1> c;

    for (int i=0; i<M-1; i++) {
      c[i] = 0.0;
      for (int j=0; j<M-1; j++)
	c[i] += alphaInv(i,j) * beta[j];
    }
    Coefs[0] = c[0];
    Coefs[1] = cusp;
    for (int i=2; i<M; i++)
      Coefs[i] = c[i-1];
    dCoefs  = 0.0;
    d2Coefs = 0.0;
    for (int i=0; i<M-1; i++) 
      dCoefs[i] = (double)(i+1) * Coefs[i+1];
    for (int i=0; i<M-2; i++)
      d2Coefs[i] = (double)(i+1) * dCoefs[i+1];
  }

  template<int M> void 
  ExpFitClass<M>::Fit (vector<double> &r, vector<double> &u)
  {
    int N = r.size();

    if (r.size() != u.size()) 
      app_error() << "Different number of rows of basis functions than"
		  << " of data points in LinFit.  Exitting.\n";
    vector<TinyVector<double,M> > F(N);
    vector<double> log_u(N);
    for (int i=0; i<N; i++) {
      log_u[i] = std::log (u[i]);
      double r2j 1.0;
      for (int j=0; j<M; j++) {
	F[i][j] = r2j;
	r2j *= r[i];
      }
    }
      
    // Next, construct alpha matrix
    Matrix<double> alpha(M,M), alphaInv(M,M);
    alpha = 0.0;
    for (int j=0; j<M; j++)
      for (int k=0; k<M; k++) {
	alpha(k,j) = 0.0;
	for (int i=0; i<N; i++)
	  alpha(k,j) += F[i][j] * F[i][k];
      }
    
    // Next, construct beta vector
    TinyVector<double,M> beta;
    beta = 0.0;
    for (int k=0; k<M; k++)
      for (int i=0; i<N; i++)
	beta[k] += log_u[i]*F[i][k];
    
    // Now, invert alpha
    for (int i=0; i<M; i++)
      for (int j=0; j<M; j++)
	alphaInv(i,j) = alpha(i,j);
    
    double det = invert_matrix(alphaInv);
	    
    for (int i=0; i<M; i++) {
      Coefs[i] = 0.0;
      for (int j=0; j<M; j++)
	Coefs[i] += alphaInv(i,j) * beta[j];
    }
    dCoefs  = 0.0;
    d2Coefs = 0.0;
    for (int i=0; i<M-1; i++) 
      dCoefs[i] = (double)(i+1) * Coefs[i+1];
    for (int i=0; i<M-2; i++)
      d2Coefs[i] = (double)(i+1) * dCoefs[i+1];
  }
    


  
  template<int M> void
  ExpFitClass<M>::eval (double r, double &u)
  {
    double r2j = 1.0;
    double poly = 0.0;
    for (int j=0; j<M; j++) {
      poly += Coefs[j] * r2j;
      rj2 *= r;
    }
    return std::exp(poly);
  }

  template<int M> void
  ExpFitClass<M>::eval (double r, double &u, double &du, double &d2u) {
    double r2j = 1.0;
    double P=0.0, dP=0.0, d2P=0.0;
    for (int j=0; j<M; j++) {
      P   +=   Coefs[j] * r2j;
      dP  +=  dCoefs[j] * r2j;      
      d2P += d2Coefs[j] * r2j;
      r2j *= r;
    }
    u   = std::exp (P);
    du  = dP * u;
    d2u = (d2P + dP*dP)*u;
  }




  // This class stores and evaluates LAPW+LO type functions inside the
  // muffin tin for a particular atom
  class MuffinTinClass
  {
  private:
    typedef QMCTraits::RealType RealType;
    typedef CrystalLattice<RealType,OHMMS_DIM> UnitCellType;
    UnitCellType PrimLattice;
    TinyVector<double,3> Center;
    // Index is the orbital number
    vector<TinyVector<double,3> > kPoints;
    double APWRadius;
    // This is the minimum grid delta.  For grid points spaced closer
    // than this value, the second derivative on the spline is
    // numerically unstable
    double drMin;
    
    int NumOrbitals;
    
    // The maximum l-channel in the sum
    int lMax;
    // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
    vector<complex<double> > YlmVec, dYlmVec;
    
    // The nonuniform radial grid for the APW splines
    NUgrid *RadialGrid;

    // There are NumOrbitals * Num_Ylm splines.  One can think of this
    // as a matrix of splines.  These splines include both the APW and
    // local orbital contribtions.
    multi_NUBspline_1d_z *RadialSplines;

    // For r smaller than rSmall, we use the polynomial fit below
    int iSmall;
    double rSmall;
    // These are coefficients of a quadratic polynomial used to
    // replace the radial splines at very small r.
    vector<TinyVector<complex<double>,3> > SmallrCoefs;
    
    // This is a helper function for fitting the small-r values
    void LinFit (vector<double> &y,  vector<TinyVector<double,2> > &F,
		 TinyVector<double,2> &a );
    void LinFit (vector<double> &y,  vector<TinyVector<double,3> > &F,
		 TinyVector<double,3> &a );

    // Temporary store for evaluating the splines
    vector<complex<double> > RadialVec, dRadialVec, d2RadialVec;
    // Evaluates all the Ylm's up to lMax
    void evalYlm(TinyVector<double,3> rhat);
    
    /////////////////
    // Core states //
    /////////////////
    // The number of core-state orbitals
    int NumCore;
    // Nonuniform spline for storing core orbitals
    vector<NUBspline_1d_d*> CoreSplines;
    // This is the radius below which we will use the polynomial fit.
    double rSmallCore;
    // Exponential fits for small and large r
    vector<ExpFitClass<4> > Small_r_Fits;
    vector<ExpFitClass<2> > Large_r_Fits;
    // Stores the polynomial fit for small r
    vector<TinyVector<double,3> > SmallrCoreCoefs;
    // Stores the expontential fit for large r
    vector<TinyVector<double,2> > LargerCoreCoefs;
    // Stores the l and m for each core state
    vector<TinyVector<int,2> > Core_lm;
    // Stores the k-vector for the core states
    vector<TinyVector<double,3> > Core_kVecs;
    // Outside this radials, the orbital is zero
    vector<double> CoreRadii;
    
  public:
    // The ion particle set
    ParticleSet *IonSet;
    // The electron-ion distance table
    DistanceTableData *ElectronIonTable;
    // Which atom this tin corresponds to
    int Atom;

    ///////////////////////////////////
    // Augmented plane-wave routines //
    ///////////////////////////////////
    void set_lattice (Tensor<RealType,3> lattice);
    void set_center  (TinyVector<double,3> center);
    void set_APW_radius (RealType radius);
    void set_APW_num_points (int num_points);
    void init_APW (Vector<double> rgrid, 
		   int lmax, int numOrbitals);
    // The first index of u_lm is l*(l+1)+m.  The second is the radial index.
    void set_APW (int orbNum, TinyVector<double,3> k,
		  Array<complex<double>,2> &u_lm, 
		  double Z);
    
    bool inside (TinyVector<double,3> r);
    void evaluate (TinyVector<double,3> r, Vector<complex<double> > &phi);
    void evaluate (TinyVector<double,3> r,
		   Vector<complex<double> > &phi,
		   Vector<TinyVector<complex<double>,3> > &grad,
		   Vector<complex<double> > &lapl);
    void evaluateFD (TinyVector<double,3> r,
		     Vector<complex<double> > &phi,
		     Vector<TinyVector<complex<double>,3> > &grad,
		     Vector<complex<double> > &lapl);
    inline int get_num_orbitals() { return NumOrbitals; }
    
    /////////////////////////
    // Core state routines //
    /////////////////////////
    inline int get_num_core() { return NumCore; }
    void addCore (int l, int m, Vector<double> &r, Vector<double> &g0,
		  TinyVector<double,3> k, double Z);
    void evaluateCore (TinyVector<double,3> r, 
		       Vector<complex<double> > &phi, int first=0);
    void evaluateCore (TinyVector<double,3> r, 
		       Vector<complex<double> > &phi,
		       Vector<TinyVector<complex<double>,3> > &grad,
		       Vector<complex<double> > &lapl,
		       int first=0);

    friend class LAPWClass;
    MuffinTinClass() : RadialSplines(NULL), CoreSplines(NULL),
		       APWRadius(0.0), NumOrbitals(0), NumCore(0),
		       lMax(0), drMin(1.0e-4)
      
    {
      
    }
    ~MuffinTinClass() 
    {
      if (RadialSplines)
	destroy_Bspline (RadialSplines);
      for (int i=0; i<CoreSplines.size(); i++)
      if (CoreSplines[i])
	destroy_Bspline (CoreSplines[i]);
    }
  };
}


#endif 
