//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#ifndef OHMMS_QMC_RADIALGRIDFUNCTOR_GAUSSIAN_H
#define OHMMS_QMC_RADIALGRIDFUNCTOR_GAUSSIAN_H

#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBasis.h"
#include "QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.h"

namespace ohmmsqmc {

  template<class T>
  struct BasicGaussian {
    T Sigma;
    T Coeff;
    T CoeffP;
    BasicGaussian(): Sigma(1.0), Coeff(1.0) { } 
    inline BasicGaussian(T sig, T c) { 
      reset(sig,c);
    } 
    void reset(T sig, T c);
    inline void setgrid(T r) { }
    inline T f(T r2) {
      return Coeff*exp(-Sigma*r2);
    }
    inline T df(T r, T r2) {
      return CoeffP*r*exp(-Sigma*r2);
    }
  };

  template<class T>
  struct GaussianCombo {
    typedef T value_type;
    ///Boolean
    bool Normalized;

    T L;
    T NormL;
    T NormPow;
    std::vector<xmlNodePtr> InParam;
    std::vector<BasicGaussian<T>* > gset;
    explicit GaussianCombo(int l, bool normalized);
    void reset();
    inline value_type f(value_type r) {
      value_type res=0;
      value_type r2 = r*r;
      for(int i=0; i<gset.size(); i++) res += gset[i]->f(r2);
      return res;
    }
    inline value_type df(value_type r) {
      value_type res=0;
      value_type r2 = r*r;
      for(int i=0; i<gset.size(); i++) res += gset[i]->df(r,r2);
      return res;
    }
    bool put(xmlNodePtr cur);
  };


  /**Class to convert GaussianTypeOrbital to a radial orbital on a log grid.
   *
   * For a center,
   *   - only one grid is used
   *   - any number of radial orbitals 
   */
  struct GTO2GridBuilder: public RGFBuilderBase {
    ///Boolean
    bool Normalized;
    ///constructor
    GTO2GridBuilder(bool normalized=false):Normalized(normalized){}
    //bool addGrid(xmlNodePtr cur);
    bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

    bool addGrid(xmlNodePtr cur);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
