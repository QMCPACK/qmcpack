//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRBFEEHANLDERTEMP_H
#define QMCPLUSPLUS_LRBFEEHANLDERTEMP_H

#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"

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
struct LRRPABFeeHandlerTemp: public LRHandlerBase
{

  DECLARE_COULOMB_TYPES

  //Typedef for the lattice-type.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  typedef BreakupBasis BreakupBasisType;

  bool FirstTime;
  mRealType rs;
  mRealType kc;
  BreakupBasisType Basis; //This needs a Lattice for the constructor...
  Func myFunc;

  //Constructor
  LRRPABFeeHandlerTemp(ParticleSet& ref, mRealType kc_in=-1.0):
    LRHandlerBase(kc_in),FirstTime(true), Basis(ref.Lattice)
  {
    myFunc.reset(ref);
  }

  //LRHandlerTemp(ParticleSet& ref, mRealType rs, mRealType kc=-1.0): LRHandlerBase(kc), Basis(ref.Lattice)
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
  LRRPABFeeHandlerTemp(const LRRPABFeeHandlerTemp& aLR, ParticleSet& ref):
    LRHandlerBase(aLR), FirstTime(true), Basis(aLR.Basis, ref.Lattice)
  {
    myFunc.reset(ref);
    fillFk(ref.SK->KLists);
  }

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new LRRPABFeeHandlerTemp<Func,BreakupBasis>(*this,ref);
  }

  void initBreakup(ParticleSet& ref)
  {
    InitBreakup(ref.Lattice,1);
    fillFk(ref.SK->KLists);
    LR_rc=Basis.get_rc();
  }

  void Breakup(ParticleSet& ref, mRealType rs_ext)
  {
    //ref.Lattice.Volume=ref.getTotalNum()*4.0*M_PI/3.0*rs*rs*rs;
    rs=rs_ext;
    myFunc.reset(ref,rs);
    InitBreakup(ref.Lattice,1);
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
    mRealType v=0.0;
    for(int n=0; n<coefs.size(); n++)
      v += coefs[n]*Basis.h(n,r);
    return v;
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv)
  {
    mRealType df = 0.0;
    //mRealType df = myFunc.df(r, rinv);
    for(int n=0; n<coefs.size(); n++)
      df += coefs[n]*Basis.df(n,r);
    return df;
  }


  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r)
  {
    mRealType vk=0.0;
    return vk;
//       for(int n=0; n<coefs.size(); n++) v -= coefs[n]*Basis.h(n,r);
  }

  /** evaluate \f$\sum_k F_{k} \rho^1_{-{\bf k} \rho^2_{\bf k}\f$
   * @param kshell degeneracies of the vectors
   * @param rk1 starting address of \f$\rho^1_{{\bf k}\f$
   * @param rk2 starting address of \f$\rho^2_{{\bf k}\f$
   *
   * Valid for the strictly ordered k and \f$F_{k}\f$.
   */
  inline mRealType evaluate(const std::vector<int>& kshell,
                           const pComplexType* restrict rk1, const pComplexType* restrict rk2)
  {
    mRealType vk=0.0;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType u=0;
      for(; ki<kshell[ks+1]; ki++,rk1++,rk2++)
        u += ((*rk1).real()*(*rk2).real()+(*rk1).imag()*(*rk2).imag());
      vk += Fk_symm[ks]*u;
    }
    //for(int ki=0; ki<Fk.size(); ki++) {
    //  //vk += (rk1[ki]*rk2[minusk[ki]]).real()*Fk[ki];
    //  vk += (rk1[ki].real()*rk2[ki].real()+rk1[ki].imag()*rk2[ki].imag())*Fk[ki];
    //} //ki
    return vk;
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
    //mRealType FatK;
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
//       std::cout<<" finding kc:  "<<ref.LR_kc<<" , "<<LR_kc<< std::endl;
    //Find size of basis from cutoffs
    kc = (LR_kc<0)?ref.LR_kc:LR_kc;
    //mRealType kc(ref.LR_kc); //User cutoff parameter...
    //kcut is the cutoff for switching to approximate k-point degeneracies for
    //better performance in making the breakup. A good bet is 30*K-spacing so that
    //there are 30 "boxes" in each direction that are treated with exact degeneracies.
    //Assume orthorhombic cell just for deriving this cutoff - should be insensitive.
    //K-Spacing = (kpt_vol)**1/3 = 2*pi/(cellvol**1/3)
    mRealType kcut = (30.0)*2*M_PI*std::pow(Basis.get_CellVolume(),-1.0/3.0);
    //Use 3000/LMax here...==6000/rc for non-ortho cells
    mRealType kmax(6000.0/ref.LR_rc);
//       std::cout<<"K_STATS !!!  "<<kcut<<"  "<<kmax<<std::endl;
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
    breakuphandler.DoBreakup(Fk.data(),coefs.data()); //Fill array of coefficients.
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
//       std::cout<<"Filling FK :"<<std::endl;
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      mRealType k=std::pow(KList.ksq[ki],0.5);
      mRealType uk=evalFk(k);
      Fk_symm[ks]=uk;
//         std::cout<<uk<<std::endl;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
    }
    //for(int ki=0; ki<KList.kpts_cart.size(); ki++){
    //  mRealType k=dot(KList.kpts_cart[ki],KList.kpts_cart[ki]);
    //  k=std::sqrt(k);
    //  Fk[ki] = evalFk(k); //Call derived fn.
    //}
  }
};
}
#endif
