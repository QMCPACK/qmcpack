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
#include <einspline/bspline_structs.h>
#include <einspline/multi_bspline_structs.h>
#include "Configuration.h"

namespace qmcplusplus {
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
    
    int NumOrbitals;
    
    // The maximum l-channel in the sum
    int lMax;
    // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
    vector<complex<double> > YlmVec, dYlmVec;
    
    // There are NumOrbitals * Num_Ylm splines.  One can think of this
    // as a matrix of splines.  These splines include both the APW and
    // local orbital contribtions.
    multi_UBspline_1d_z *RadialSplines;
    
    // Temporary store for evaluating the splines
    vector<complex<double> > RadialVec, dRadialVec, d2RadialVec;
    // Evaluates all the Ylm's up to lMax
    void evalYlm(TinyVector<double,3> rhat);
    
    /////////////////
    // Core states //
    /////////////////
    // The number of core-state orbitals
    int NumCore;
    vector<UBspline_1d_d*> CoreSplines;
    // Stores the l and m for each core state
    vector<TinyVector<int,2> > Core_lm;
    // Stores the k-vector for the core states
    vector<TinyVector<double,3> > Core_kVecs;
    // Outside this radials, the orbital is zero
    double CoreRadius;
    
  public:
    ///////////////////////////////////
    // Augmented plane-wave routines //
    ///////////////////////////////////
    void set_lattice (Tensor<RealType,3> lattice);
    void set_center  (TinyVector<double,3> center);
    void set_APW_radius (RealType radius);
    void set_APW_num_points (int num_points);
    void init_APW (double radius, int num_rad_points, 
		   int lmax, int numOrbitals);
    // The first index of u_lm is l*(l+1)+m.  The second is the radial index.
    void set_APW (int orbNum, TinyVector<double,3> k,
		  Array<complex<double>,2> &u_lm);
    
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
    void addCore (int l, int m, double rmax, Vector<double> &g0,
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
		       lMax(0), CoreRadius(0.0)
      
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
