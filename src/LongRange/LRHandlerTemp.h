//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRHANLDERTEMP_H
#define QMCPLUSPLUS_LRHANLDERTEMP_H
#include "coulomb_types.h"
#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{

/* Templated LRHandler class
 *
 * LRHandlerTemp<Func,BreakupBasis> is a modification of LRHandler
 * and a derived class from LRHanlderBase.
 * The first template parameter Func is a generic functor, e.g., CoulombFunctor.
 * The second template parameter is a BreakupBasis and the default is set to LPQHIBasis.
 * LRHandlerBase is introduced to enable run-time options. See RPAContstraints.h
 */
template<class Func, class BreakupBasis=LPQHIBasis>
class LRHandlerTemp: public LRHandlerBase
{

public:
  //Typedef for the lattice-type.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  typedef BreakupBasis BreakupBasisType;

  bool FirstTime;
  mRealType rs;
  BreakupBasisType Basis; //This needs a Lattice for the constructor...
  Func myFunc;

  //Constructor
  LRHandlerTemp(ParticleSet& ref, mRealType kc_in=-1.0):
    LRHandlerBase(kc_in),FirstTime(true), Basis(ref.LRBox)
  {
    myFunc.reset(ref);
  }

  //LRHandlerTemp(ParticleSet& ref, mRealType rs, mRealType kc=-1.0): LRHandlerBase(kc), Basis(ref.LRBox)
  //{
  //  myFunc.reset(ref,rs);
  //}

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  LRHandlerTemp(const LRHandlerTemp& aLR, ParticleSet& ref):
    LRHandlerBase(aLR), FirstTime(true), Basis(aLR.Basis, ref.LRBox)
  {
    myFunc.reset(ref);
    fillFk(ref.SK->KLists);
  }

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new LRHandlerTemp<Func,BreakupBasis>(*this,ref);
  }

  void initBreakup(ParticleSet& ref)
  {
    InitBreakup(ref.LRBox,1);
    fillFk(ref.SK->KLists);
    LR_rc=Basis.get_rc();
  }

  void Breakup(ParticleSet& ref, mRealType rs_ext)
  {
    //ref.LRBox.Volume=ref.getTotalNum()*4.0*M_PI/3.0*rs*rs*rs;
    rs=rs_ext;
    myFunc.reset(ref,rs);
    InitBreakup(ref.LRBox,1);
    fillFk(ref.SK->KLists);
    LR_rc=Basis.get_rc();
  }

  void resetTargetParticleSet(ParticleSet& ref)
  {
    myFunc.reset(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref, mRealType rs)
  {
    myFunc.reset(ref,rs);
  }

  inline mRealType evaluate(mRealType r, mRealType rinv)
  {
    mRealType v=myFunc(r,rinv);
    for(int n=0; n<coefs.size(); n++)
      v -= coefs[n]*Basis.h(n,r);
    return v;
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv)
  {
    mRealType df = myFunc.df(r);
    //RealType df = myFunc.df(r, rinv);
    for(int n=0; n<coefs.size(); n++)
      df -= coefs[n]*Basis.df(n,r);
    return df;
  }


  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r)
  {
    mRealType v=0.0;
    for(int n=0; n<coefs.size(); n++)
      v -= coefs[n]*Basis.h(n,r);
    return v;
  }

  inline mRealType evaluateSR_k0()
  {
    mRealType v0=myFunc.integrate_r2(Basis.get_rc());
    for(int n=0; n<coefs.size(); n++)
      v0 -= coefs[n]*Basis.hintr2(n);
    return v0*2.0*TWOPI/Basis.get_CellVolume();
  }

  inline mRealType evaluateLR_r0()
  {
    mRealType v0=0.0;
    for(int n=0; n<coefs.size(); n++)
      v0 += coefs[n]*Basis.h(n,0.0);
    return v0;
  }

private:

  inline mRealType evalFk(mRealType k)
  {
    //FatK = 4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
    mRealType FatK=myFunc.Fk(k,Basis.get_rc());
    for(int n=0; n<Basis.NumBasisElem(); n++)
      FatK += coefs[n]*Basis.c(n,k);
    return FatK;
  }
  inline mRealType evalXk(mRealType k)
  {
    //RealType FatK;
    //FatK = -4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
    //return (FatK);
    return myFunc.Xk(k,Basis.get_rc());
  }

  /** Initialise the basis and coefficients for the long-range beakup.
   *
   * We loocally create a breakup handler and pass in the basis
   * that has been initialised here. We then discard the handler, leaving
   * basis and coefs in a usable state.
   * This method can be re-called later if lattice changes shape.
   */
  void InitBreakup(ParticleLayout_t& ref,int NumFunctions)
  {
    //First we send the new Lattice to the Basis, in case it has been updated.
    Basis.set_Lattice(ref);
    //Compute RC from box-size - in constructor?
    //No here...need update if box changes
    int NumKnots(15);
    Basis.set_NumKnots(NumKnots);
    Basis.set_rc(ref.LR_rc);
    //Initialise the breakup - pass in basis.
    LRBreakup<BreakupBasis> breakuphandler(Basis);
    //Find size of basis from cutoffs
    mRealType kc = (LR_kc<0)?ref.LR_kc:LR_kc;
    //RealType kc(ref.LR_kc); //User cutoff parameter...
    //kcut is the cutoff for switching to approximate k-point degeneracies for
    //better performance in making the breakup. A good bet is 30*K-spacing so that
    //there are 30 "boxes" in each direction that are treated with exact degeneracies.
    //Assume orthorhombic cell just for deriving this cutoff - should be insensitive.
    //K-Spacing = (kpt_vol)**1/3 = 2*pi/(cellvol**1/3)
    mRealType kcut = 60*M_PI*std::pow(Basis.get_CellVolume(),-1.0/3.0);
    //Use 3000/LMax here...==6000/rc for non-ortho cells
    mRealType kmax(6000.0/ref.LR_rc);
    MaxKshell = static_cast<int>(breakuphandler.SetupKVecs(kc,kcut,kmax));
    if(FirstTime)
    {
      app_log() <<" finding kc:  "<<ref.LR_kc<<" , "<<LR_kc<< std::endl;
      app_log() << "  LRBreakp parameter Kc =" << kc << std::endl;
      app_log() << "    Continuum approximation in k = [" << kcut << "," << kmax << ")" << std::endl;
      FirstTime=false;
    }
    //Set up x_k
    //This is the FT of -V(r) from r_c to infinity.
    //This is the only data that the breakup handler needs to do the breakup.
    //We temporarily store it in Fk, which is replaced with the full FT (0->inf)
    //of V_l(r) after the breakup has been done.
    fillXk(breakuphandler.KList);
    //Allocate the space for the coefficients.
    coefs.resize(Basis.NumBasisElem()); //This must be after SetupKVecs.
    
    mRealType chisqr(0.0);
    chisqr=breakuphandler.DoBreakup(Fk.data(),coefs.data()); //Fill array of coefficients.
    app_log()<<"\n   LR Breakup chi^2 = "<<chisqr<<std::endl;
  }

  void fillXk(std::vector<TinyVector<mRealType,2> >& KList)
  {
    Fk.resize(KList.size());
    for(int ki=0; ki<KList.size(); ki++)
    {
      mRealType k=KList[ki][0];
      Fk[ki] = evalXk(k); //Call derived fn.
    }
  }

  void fillFk(KContainer& KList)
  {
    Fk.resize(KList.kpts_cart.size());
    const std::vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    Fk_symm.resize(MaxKshell);
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      mRealType uk=evalFk(std::sqrt(KList.ksq[ki]));
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
    }
    //for(int ki=0; ki<KList.kpts_cart.size(); ki++){
    //  RealType k=dot(KList.kpts_cart[ki],KList.kpts_cart[ki]);
    //  k=std::sqrt(k);
    //  Fk[ki] = evalFk(k); //Call derived fn.
    //}
  }
};
}
#endif
