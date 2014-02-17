//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file LRHandlerSRCoulomb.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRHANLDERSRCOULOMBTEMP_H
#define QMCPLUSPLUS_LRHANLDERSRCOULOMBTEMP_H

#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHISRCoulombBasis.h"
#include "LongRange/LRBreakup.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{

/* Templated LRHandler class
 *
 * LRHandlerSRCoulomb<Func,BreakupBasis> is a modification of LRHandler
 * and a derived class from LRHanlderBase. Implements the LR breakup http://dx.doi.org/10.1006/jcph.1995.1054 .
 * The first template parameter Func is a generic functor, e.g., CoulombFunctor.
 * The second template parameter is a BreakupBasis and the default is set to LPQHIBasis.
 * LRHandlerBase is introduced to enable run-time options. See RPAContstraints.h
 */
template<class Func, class BreakupBasis=LPQHISRCoulombBasis>
class LRHandlerSRCoulomb: public LRHandlerBase
{

public:
  //Typedef for the lattice-type.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  typedef BreakupBasis BreakupBasisType;

  bool FirstTime;
  RealType rs;
  BreakupBasisType Basis; //This needs a Lattice for the constructor...
  Func myFunc;

  //Constructor
  LRHandlerSRCoulomb(ParticleSet& ref, RealType kc_in=-1.0):
    LRHandlerBase(kc_in),FirstTime(true), Basis(ref.LRBox)
  {
    myFunc.reset(ref);
  }

  //LRHandlerSRCoulomb(ParticleSet& ref, RealType rs, RealType kc=-1.0): LRHandlerBase(kc), Basis(ref.LRBox)
  //{
  //  myFunc.reset(ref,rs);
  //}

  /** "copy" constructor
   * @param aLR LRHandlerSRCoulomb
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  LRHandlerSRCoulomb(const LRHandlerSRCoulomb& aLR, ParticleSet& ref):
    LRHandlerBase(aLR), FirstTime(true), Basis(aLR.Basis, ref.LRBox)
  {
    myFunc.reset(ref);
    fillYk(ref.SK->KLists);
    fillYkg(ref.SK->KLists);
  }

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new LRHandlerSRCoulomb<Func,BreakupBasis>(*this,ref);
  }

  void initBreakup(ParticleSet& ref)
  {
    InitBreakup(ref.LRBox,1);
    fillYk(ref.SK->KLists);
    fillYkg(ref.SK->KLists);
    LR_rc=Basis.get_rc();
  }

  void Breakup(ParticleSet& ref, RealType rs_ext)
  {
    //ref.LRBox.Volume=ref.getTotalNum()*4.0*M_PI/3.0*rs*rs*rs;
    rs=rs_ext;
    myFunc.reset(ref,rs);
    InitBreakup(ref.LRBox,1);
    fillYk(ref.SK->KLists);
    fillYkg(ref.SK->KLists);
    LR_rc=Basis.get_rc();
  }

  void resetTargetParticleSet(ParticleSet& ref)
  {
    myFunc.reset(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref, RealType rs)
  {
    myFunc.reset(ref,rs);
  }

  inline RealType evaluate(RealType r, RealType rinv)
  {
    RealType v=0.0;
    for(int n=0; n<coefs.size(); n++)
      v += coefs[n]*Basis.h(n,r);
    return v;
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline RealType srDf(RealType r, RealType rinv)
  {
    RealType df = 0.0;
    //RealType df = myFunc.df(r, rinv);
    for(int n=0; n<coefs.size(); n++)
      df += gcoefs[n]*Basis.df(n,r);
    return df;
  }


  /** evaluate the contribution from the long-range part for for spline
   */
  inline RealType evaluateLR(RealType r)
  {
    RealType v=0.0;
    for(int n=0; n<coefs.size(); n++)
      v -= coefs[n]*Basis.h(n,r);
    return v;
  }

  inline RealType evaluateSR_k0()
  {
    RealType v0=myFunc.integrate_r2(Basis.get_rc());
    for(int n=0; n<coefs.size(); n++)
      v0 -= coefs[n]*Basis.hintr2(n);
    return v0*2.0*TWOPI/Basis.get_CellVolume();
  }

  inline RealType evaluateLR_r0()
  {
    RealType v0=0.0;
    for(int n=0; n<coefs.size(); n++)
      v0 += coefs[n]*Basis.h(n,0.0);
    return v0;
  }

private:

  inline RealType evalYk(RealType k)
  {
    //FatK = 4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
    RealType FatK=myFunc.Vk(k);
    for(int n=0; n<Basis.NumBasisElem(); n++)
      FatK -= coefs[n]*Basis.c(n,k);
    return FatK;
  }
  inline RealType evalYkg(RealType k)
  {
    RealType FatK=myFunc.Vk(k);
    for(int n=0; n<Basis.NumBasisElem(); n++)
      FatK -= gcoefs[n]*Basis.c(n,k);
    return FatK;
    
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
    RealType chisqr_f=0.0;
    RealType chisqr_df=0.0;
    
    //First we send the new Lattice to the Basis, in case it has been updated.
    Basis.set_Lattice(ref);
    //Compute RC from box-size - in constructor?
    //No here...need update if box changes
    int NumKnots(21);
    Basis.set_NumKnots(NumKnots);
    Basis.set_rc(ref.LR_rc);
    //Initialise the breakup - pass in basis.
    LRBreakup<BreakupBasis> breakuphandler(Basis);
    //Find size of basis from cutoffs
    RealType kc = (LR_kc<0)?ref.LR_kc:LR_kc;
    //RealType kc(ref.LR_kc); //User cutoff parameter...
    //kcut is the cutoff for switching to approximate k-point degeneracies for
    //better performance in making the breakup. A good bet is 30*K-spacing so that
    //there are 30 "boxes" in each direction that are treated with exact degeneracies.
    //Assume orthorhombic cell just for deriving this cutoff - should be insensitive.
    //K-Spacing = (kpt_vol)**1/3 = 2*pi/(cellvol**1/3)
    RealType kcut = 60*M_PI*std::pow(Basis.get_CellVolume(),-1.0/3.0);
    //Use 3000/LMax here...==6000/rc for non-ortho cells
    RealType kmax(6000.0/ref.LR_rc);
    MaxKshell = static_cast<int>(breakuphandler.SetupKVecs(kc,kcut,kmax));
    if(FirstTime)
    {
      app_log() <<" finding kc:  "<<ref.LR_kc<<" , "<<LR_kc<<endl;
      app_log() << "  LRBreakp parameter Kc =" << kc << endl;
      app_log() << "    Continuum approximation in k = [" << kcut << "," << kmax << ")" << endl;
      FirstTime=false;
    }
    //Set up x_k
    //This is the FT of -V(r) from r_c to infinity.
    //This is the only data that the breakup handler needs to do the breakup.
    //We temporarily store it in Fk, which is replaced with the full FT (0->inf)
    //of V_l(r) after the breakup has been done.
    fillVk(breakuphandler.KList);
    //Allocate the space for the coefficients.
    IndexType Nbasis=Basis.NumBasisElem();
    coefs.resize(Nbasis); //This must be after SetupKVecs.
    gcoefs.resize(Nbasis);

    //Going to implement a smooth real space cutoff.  This means that alpha=0,1,2 for the LPQHI basis at knot r_c
    //all equal the 0, 1st, and 2nd derivatives of our bare function.  
    //These three functions are the last three basis elements in our set.  


   
    Vector<RealType> constraints;
    constraints.resize(Nbasis);
    for (int i=0; i < Nbasis; i++) constraints[i]=1;
    

    RealType rc=Basis.get_rc();
    
    ///This is to make sure there's no cusp in the LR part.  
    gcoefs[0]=coefs[0] = 1.0;
    constraints[0]=0;
    gcoefs[1] = coefs[1] = 0.0;
    constraints[1]=0;

    gcoefs[2] = coefs[2] = 0.0; 
    constraints[2]=0.0;

    //2nd derivative of SR will go to zero by setting LR at rc equal to the bare function.
   // gcoefs[Nbasis-1]=myFunc.df2(rc);


    gcoefs[Nbasis-1]=coefs[Nbasis-1]=0.0;
    constraints[Nbasis-1]=0;
   
    //1st derivative
    
    gcoefs[Nbasis-2]=coefs[Nbasis-2]=0.0;
    constraints[Nbasis-2]=0;

    //Function value 
    gcoefs[Nbasis-3]=coefs[Nbasis-3]=0.0;
    constraints[Nbasis-3]=0;
    //And now to impose the constraints
    

    
    chisqr_f=breakuphandler.DoBreakup(Fk.data(),coefs.data(),constraints.data()); //Fill array of coefficients.
    chisqr_df=breakuphandler.DoGradBreakup(Fkg.data(), gcoefs.data(), constraints.data());
    
    app_log()<<"LR function chi^2 = "<<chisqr_f<<endl;
    app_log()<<"LR grad function chi^2 = "<<chisqr_df<<endl;
    app_log()<<"  n  tn   gtn h(n)\n";

    for (int i=0; i<Nbasis; i++)
    {
       app_log()<<"  "<<i<<" "<<coefs[i]<<" "<<gcoefs[i]<<" "<<Basis.rh(i,i*2.0/(14.0))<<endl;
    }
    RealType delta=rc/RealType(NumKnots-1);

    for (int i=0;i<NumKnots; i++)
	app_log()<<i*delta<<" "<<evaluate(i*delta, 1.0/(delta*i))<<" "<<srDf(i*delta, 1.0/(delta*i))<<endl;
    

  }


  void fillVk(vector<TinyVector<RealType,2> >& KList)
  {
    Fk.resize(KList.size());
    Fkg.resize(KList.size());
    for(int ki=0; ki<KList.size(); ki++)
    {
      RealType k=KList[ki][0];
      Fk[ki] = myFunc.Vk(k); //Call derived fn.
      Fkg[ki]= myFunc.Vk(k);
    }
  }

  void fillYk(KContainer& KList)
  {
    Fk.resize(KList.kpts_cart.size());
    const vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    Fk_symm.resize(MaxKshell);
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      RealType uk=evalYk(std::sqrt(KList.ksq[ki]));
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
  void fillYkg(KContainer& KList)
  {
    Fkg.resize(KList.kpts_cart.size());
    const vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    Fk_symmg.resize(MaxKshell);
    for(int ks=0,ki=0; ks<Fk_symmg.size(); ks++)
    {
      RealType uk=evalYkg(std::sqrt(KList.ksq[ki]));
      Fk_symmg[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fkg.size())
        Fkg[ki++]=uk;
    }
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5981 $   $Date: 2013-09-17 14:47:32 -0500 (Tue, 17 Sep 2013) $
 * $Id: LRHandlerSRCoulomb.h 5981 2013-09-17 19:47:32Z jnkim $
 ***************************************************************************/
