//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////
     
/** @file CasinoJastrowFunctors.h
 * @brief Functors which implement jastrow factor from Drummond, Towler and Needs, PRB 70, 235119 (2004)
 * \f[ u(r) = (r - L)^c Theta(L - r) (alpha_0 + [Gamma / (-L)^c + alpha_0 * c / L]*r + \sum(l=2,N) alpha_l * r^l)
 * in the code L == cuotff_
 *             c == derivContinuityOrder_
 */

#ifndef QMCPLUSPLUS_CASINOFUNCTORS_H
#define QMCPLUSPLUS_CASINOFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
#include "simd/allocator.hpp"
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{

template<class T>
struct casinoFunctor:public OptimizableFunctorBase
{
 private:
  /// if true can optimize the coefficients
  bool optParams_;
  /// parameter giving number of differentials for continuity at the cutoff
  int derivContinuityOrder_;
  /// parameter giving the cutoff
  real_type cutoff_;
  /// Value for the cusp
  real_type cusp_;
  /// collection of values that are the prefactor to r
  real_type prefactor_;
  
  /// Values for alpha coefficients (note alpha_1 is not present, so after 0 we jump straight to alpha_2)
  std::vector<real_type> alphas_;
  /// Names for alpha coefficients
  std::vector<std::string> parameterNames_;

 public:
  bool periodic_;
  /// Default constructor
 casinofunctor(): { }
  
  // constructor
  explicit casinofunctor(real_type cusp, int cont) : optParams_(true), derivContinuityOrder_(cont), cutoff_(-1), cusp_(cusp), periodic_(true)
  {
    reset();
  }
    
  OptimizableFunctorBase* makeClone() const
  {
    return new casinofunctor(*this);
  }

  void reset()
  {
    if (alphas_.size() > 0) {
      prefactor_ = cusp_ / pow(-cutoff_,derivContinuityOrder_) + alphas_[0] * derivContinuityOrder_ / cutoff_;
    }
  }

  
  real_type series(real_type r) const
  {
    real_type result = 0;
    real_type cumr = r*r;
    for (int i = 1; i < alphas_.size(); i++) 
      {
	result += alphas_[i]*cumr;
	cumr *= r;
      }
    return result;
  }

  real_type series(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    // recall that alphas_[1] corresponds to a_2
    real_type result = alphas_[1]*r*r;
    dudr = 2*alphas_[1]*r;
    d2udr2 = 2*alphas_[1];

    real_type cumr = r;
    // first term in the series is for alpha_3 * r^3
    for (int i = 2; i < alphas_.size(); i++) 
    {
      result += alphas_[i]*cumr*r*r;
      dudr   += (i+1)*alphas_[i]*cumr*r;
      d2udr2 += (i+1)*i*alphas_[i]*cumr;
      cumr *= r;
    }
    return result;
  }

  real_type series(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    // recall that alphas_[1] corresponds to a_2
    real_type result = alphas_[1]*r*r;
    dudr = 2*alphas_[1]*r;
    d2udr2 = 2*alphas_[1];
    du3dr3 = 0;

    real_type cumr = 1;
    const real_type r3 = r*r*r;
    const real_type r2 = r*r;

    // first term in the series is for alpha_3 * r^3
    for (int i = 2; i < alphas_.size(); i++) 
    {
      result += alphas_[i]*cumr*r3;
      dudr   += (i+1)*alphas_[i]*cumr*r2;
      d2udr2 += (i+1)*i*alphas_[i]*cumr*r;
      du3dr3 += (i+1)*i*(i-1)*alphas_[i]*cumr;
      cumr *= r;
    }
    return result;
  }


  inline real_type evaluate(real_type r) const
  {
    if (r >= cutoff_)
      {
      return 0.0;
      }
    return pow(r-cutoff_,derivContinuityOrder_)*(alphas_[0]+prefactor_*r+series(r));
  }
  
  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    if (r >= cutoff_) 
    {
      dudr=d2udr2=0.0;
      return 0.0;
    }
    // Need to fill in derivatives here...
    const real_type rl = r-cutoff_;
    real_type ser, dser, d2ser;
    ser = series(r, dser, d2ser);
    
    dudr = pow(rl,derivContinuityOrder_-1) * ( alphas_[0]*derivContinuityOrder_ + prefactor_*derivContinuityOrder_*r + (dser + prefactor_)*rl );
    du2dr2 = pow(rl,derivContinuityOrder_-2) * ( 2*derivContinuityOrder_*rml*(dser+prefactor_) + derivContinuityOrder_*(derivContinuityOrder_-1)*(alphas_[0]+prefactor_*r+ser) + d2ser*rml*rml);
    return pow(rl,derivContinuityOrder_)*(alphas_[0]+prefactor_*r+ser);
  }    

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    if (r >= cutoff_) 
    {
      dudr=d2udr2=d3udr3=0.0;
      return 0.0;
    }
    // Need to fill in derivatives here...
    const real_type rl = r-cutoff_;
    real_type ser, dser, d2ser, d3ser;
    ser = series(r, dser, d2ser, d3ser);
    const real_type X = alphas_[0]+prefactor_*r+ser;

    dudr = pow(rl,derivContinuityOrder_-1) * ( alphas_[0]*derivContinuityOrder_ + prefactor_*derivContinuityOrder_*r + (dser + prefactor_)*rl );
    du2dr2 = pow(rl,derivContinuityOrder_-2) * ( 2*derivContinuityOrder_*rl*(dser+prefactor_) + derivContinuityOrder_*(derivContinuityOrder_-1)*X + d2ser*rl*rl);
    du3dr3 = pow(rl,derivContinuityOrder_-3) * ( 3*derivContinuityOrder_*rl*rl*d2ser + 3*derivContinuityOrder_*rl*(derivContinuityOrder_-1)*(dser+prefactor_) + d3ser*rl*rl*rl + derivContinuityOrder_*X*(derivContinuityOrder_*derivContinuityOrder_-3*derivContinuityOrder_+2));
    return pow(rl,derivContinuityOrder_)*(alphas_[0]+prefactor_*r+ser);
  }    

  inline real_type evaluateV(const int iat, const int iStart, const int iEnd,
    const T* restrict _distArray, T* restrict distArrayCompressed ) const
  {
    real_type sum(0);
    for(int idx=iStart; idx<iEnd; idx++)
      if (idx!=iat) sum += evaluate(_distArray[idx]);
    return sum;
  }

  inline void evaluateVGL(const int iat, const int iStart, const int iEnd,
    const T* distArray,  T* restrict valArray,
    T* restrict gradArray, T* restrict laplArray,
    T* restrict distArrayCompressed, int* restrict distIndices ) const
  {
    for(int idx=iStart; idx<iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if ( iat>=iStart && iat<iEnd )
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    if(optParams_)
    {

    if (r >= cutoff_) 
    {
      for (int i = 0; i < derivs.size(); i++)
      {
	for (int j = 0; j < 3; j++)
	{
	  derivs[i][j] = 0.0;
	}
      }
      return true;
    }

    // translation: rlTcm1 == (r minus L) To the power c minus 1
    const real_type rl = r-cutoff_;
    const real_type rlTc = pow(rl,derivContinuityOrder_);
    const real_type rlTcm1 = pow(rl,derivContinuityOrder_-1);
    const real_type rlTcm2 = pow(rl,derivContinuityOrder_-2);
    derivs[0][0] = rlTc*(1.0+derivContinuityOrder_*r/cutoff_);  // df / d alpha_0
    derivs[0][1] = rlTcm1*(derivContinuityOrder_+1)*(derivContinuityOrder_*r/cutoff_); // d / dr (df / d alpha_0)
    derivs[0][2] = rlTcm2*(derivContinuityOrder_/cutoff_)*(-cutoff_ - derivContinuityOrder_*r + 2*derivContinuityOrder_*rl + derivContinuityOrder_*(l+derivContinuityOrder_*r)); // d2 / dr2 (df / d alpha_0)

    // now do all of the derivatives with respect to alphas__i where i >=2
    for (int i = 1; i < derivs.size(); i++) {
      derivs[i][0] = rlTc * pow(r,i+1); 
      derivs[i][1] = rlTcm1 * pow(r,i) * ( (i+1)*r + derivContinuityOrder_*r - (i+1)*cutoff_ );
      derivs[i][2] = rlTcm2 * pow(r,i-1) * ( derivContinuityOrder_*(derivContinuityOrder_-1)*r*r + 2*(i+1)*derivContinuityOrder_*r*rl + (i+1)*i*rl*rl );
    }

    }
    return true;
  }
    

  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    if(optParams_)
    {

    if (r >= cutoff_) 
    {
      for (int i = 0; i < derivs.size(); i++)
      {
	derivs[i] = 0.0;
      }
      return true;
    }

    
    const real_type rl = r-cutoff_;
    const real_type rlTc = pow(rl,derivContinuityOrder_);
    derivs[0] = rlTc*(1.0+derivContinuityOrder_*r/cutoff_);  // df / d alpha_0

    // now do all of the derivatives with respect to alphas__i where i >=2
    for (int i = 1; i < derivs.size(); i++) 
    {
      derivs[i][0] = rlTc * pow(r,i+1); 
    }

    }
    return true;
  }
   
  bool put(xmlNodePtr cur)
  {
    ReportEngine PRE("CasinoJastrowFunctor", "put(xmlNodePtr)");
    int numParams = 0;
    OhmmsAttributeSet rAttrib;
    double locCutoff = -1.0;
    derivContinuityOrder_ = 2;
    rAttrib.add(numParams,             "size");
    rAttrib.add(locCutoff,            "rcut");
    rAttrib.add(locCutoff,            "cutoff");
    rAttrib.add(derivContinuityOrder_, "cont");
    rAttrib.put(cur);

    if (locCutoff < 0.0)
      if (periodic)
        app_log() << "  Jastrow cutoff unspecified.  Setting to Wigner-Seitz radius = " << cutoff_ << ".\n";
      else
      {
        APP_ABORT("  Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
      }
    else
      if (periodic && locCutoff > cutoff_)
      {
        if (locCutoff - cutoff_ > 1e-4)
        {
          APP_ABORT( "  The Jastrow cutoff specified should not be larger than Wigner-Seitz radius.");
        }
        else
        {
          app_log() << "  The Jastrow cutoff specified is slightly larger than the Wigner-Seitz radius.";
          app_log() << "  Setting to Wigner-Seitz radius = " << cutoff_ << ".\n";
        }
      }
      else
        cutoff_ = locCutoff;
    if (numParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Casino jastrow function.",true);
    }
    app_log() << " size = " << numParams << " parameters " << std::endl;
    app_log() << " cusp = " << cusp_ << std::endl;
    app_log() << " rcut = " << cutoff_ << std::endl;
    alphas_.resize(numParams);





  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  void resetParameters(const opt_variables_type& active)
  {

    // look at bspline jastrow for this

  }
};  
  

    
