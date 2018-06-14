//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerSRCoulomb.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRHANLDERSRCOULOMBTEMP_H
#define QMCPLUSPLUS_LRHANLDERSRCOULOMBTEMP_H

#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHISRCoulombBasis.h"
#include "LongRange/LRBreakup.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimLinearSpline.h"

#include <sstream>
#include <string>

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
  typedef LinearGrid<mRealType>                               GridType;
  typedef OneDimCubicSpline<mRealType>                       RadFunctorType;

  bool FirstTime;
  mRealType rs;
  BreakupBasisType Basis; //This needs a Lattice for the constructor...
  Func myFunc;

  std::vector<mRealType> Fk_copy;
  
  GridType* aGrid;
  
  RadFunctorType* rV_energy;
  RadFunctorType* rV_force;
  RadFunctorType* drV_force;
  RadFunctorType* rV_stress;
  RadFunctorType* drV_stress;
  

  //Constructor
  LRHandlerSRCoulomb(ParticleSet& ref, mRealType kc_in=-1.0):
    LRHandlerBase(kc_in),FirstTime(true), Basis(ref.LRBox), aGrid(0), 
    rV_force(0), rV_energy(0), drV_force(0), rV_stress(0), drV_stress(0)

  {
    myFunc.reset(ref);
  }
   
  ~LRHandlerSRCoulomb()
  {
//	 delete aGrid;
//	 delete rV_energy;
//     delete rV_force;
//     delete drV_force;
//     delete rV_stress;
//     delete drV_stress;  
  }
  //LRHandlerSRCoulomb(ParticleSet& ref, mRealType rs, mRealType kc=-1.0): LRHandlerBase(kc), Basis(ref.LRBox)
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
//    myFunc.reset(ref);
//    fillYk(ref.SK->KLists);
//    fillYkg(ref.SK->KLists);
//    filldFk_dk(ref.SK->KLists);
   // app_log()<<"copy constructor called.  thread #"<<omp_get_num_threads()<< std::endl;
    aGrid = new GridType(*(aLR.aGrid));
   // new GridType(aLR.aGrid)
    rV_energy = aLR.rV_energy->makeClone();
    rV_force= aLR.rV_force->makeClone();
    drV_force= aLR.drV_force->makeClone();
    rV_stress= aLR.rV_stress->makeClone();
    drV_stress= aLR.drV_stress->makeClone();
//    rV_force = new RadFunctorType(*(aLR.rV_force));
//    drV_force = new RadFunctorType(*(aLR.drV_force));
//    rV_stress = new RadFunctorType(*(aLR.rV_stress));
//    drV_stress = new RadFunctorType(*(aLR.drV_stress));  
    

  }

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    LRHandlerSRCoulomb* tmp= new LRHandlerSRCoulomb<Func,BreakupBasis>(*this,ref);
//    tmp->makeSplines(1001);
    return tmp;
  }

  void initBreakup(ParticleSet& ref)
  {
    InitBreakup(ref.LRBox,1);
    fillYk(ref.SK->KLists);
    fillYkg(ref.SK->KLists);
    filldFk_dk(ref.SK->KLists);
    LR_rc=Basis.get_rc();
    //makeSplines(1000);
  }

  void Breakup(ParticleSet& ref, mRealType rs_ext)
  {
    //ref.LRBox.Volume=ref.getTotalNum()*4.0*M_PI/3.0*rs*rs*rs;
    rs=rs_ext;
    myFunc.reset(ref,rs);
    InitBreakup(ref.LRBox,1);
    fillYk(ref.SK->KLists);
    fillYkg(ref.SK->KLists);
    filldFk_dk(ref.SK->KLists);
    LR_rc=Basis.get_rc();
   // makeSplines(1000);
  }
  void makeSplines(int ngrid)
  {
     if(aGrid == 0)
     {
       aGrid = new GridType;
       aGrid->set(0.0,Basis.get_rc(),ngrid);
     }
     
     std::vector<mRealType> vE(ngrid);
     std::vector<mRealType> vF(ngrid);
     std::vector<mRealType> dvF(ngrid);
     std::vector<mRealType> vS(ngrid);
     std::vector<mRealType> dvS(ngrid);
     
     for( int i=1; i<ngrid; i++)
     {
		mRealType r=(*aGrid)[i];
		
		vE[i]=r*Basis.f(r,coefs);
		vF[i]=r*Basis.f(r,gcoefs);
		dvF[i]= r*r*Basis.df_dr(r,gcoefs);
		vS[i]=r*Basis.f(r,gstraincoefs);
		dvS[i]= r*r*Basis.df_dr(r,gstraincoefs);
	 }
	 
	 vE[0]=1.0;
	 vF[0]=1.0;
	 dvF[0]=1.0;
	 vS[0]=1.0;
	 dvS[0]=1.0;
	 

     rV_energy=new RadFunctorType(aGrid,vE);
     rV_force=new RadFunctorType(aGrid,vF);
     drV_force=new RadFunctorType(aGrid,dvF);
     rV_stress=new RadFunctorType(aGrid,vS);
     drV_stress=new RadFunctorType(aGrid,dvS);
     
     rV_energy->spline(0,vE[0],ngrid-1,0);
     rV_force->spline(0,vF[0],ngrid-1,0);
     drV_force->spline(0,dvF[0],ngrid-1,0);
     rV_stress->spline(0,vS[0],ngrid-1,0);
     drV_stress->spline(0,dvS[0],ngrid-1,0);	  
	  
  }
 /* void makeSplines(int ngrid)
  {
     if(aGrid == 0)
     {
       aGrid = new GridType;
       aGrid->set(0.0,Basis.get_rc(),ngrid);
     }
     
     std::vector<mRealType> vE(ngrid);
     std::vector<mRealType> vF(ngrid);
     std::vector<mRealType> dvF(ngrid);
     std::vector<mRealType> vS(ngrid);
     std::vector<mRealType> dvS(ngrid);
     
     for( int i=1; i<ngrid; i++)
     {
		mRealType r=(*aGrid)[i];
		
		vE[i]=r*Basis.f(r,coefs);
		vF[i]=r*Basis.f(r,gcoefs);
		dvF[i]=r*Basis.df_dr(r,gcoefs)+Basis.f(r,gcoefs);
		vS[i]=r*Basis.f(r,gstraincoefs);
		dvS[i]= r*Basis.df_dr(r,gstraincoefs)+Basis.f(r,gstraincoefs);
	 }
	 
	 vE[0]=1.0;
	 vF[0]=1.0;
	 dvF[0]=0.0;
	 vS[0]=1.0;
	 dvS[0]=1.0;
	 

     rV_energy=new RadFunctorType(aGrid,vE);
     rV_force=new RadFunctorType(aGrid,vF);
     drV_force=new RadFunctorType(aGrid,dvF);
     rV_stress=new RadFunctorType(aGrid,vS);
     drV_stress=new RadFunctorType(aGrid,dvS);
     
     rV_energy->spline(0,vE[0],ngrid-1,vE[ngrid-1]);
     rV_force->spline(0,vF[0],ngrid-1,vF[ngrid-1]);
     drV_force->spline(0,dvF[0],ngrid-1,dvF[ngrid-1]);
     rV_stress->spline(0,vS[0],ngrid-1,vS[ngrid-1]);
     drV_stress->spline(0,dvS[0],ngrid-1,dvS[ngrid-1]);	  
	  
  }*/

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
    mRealType v = Basis.f(r, coefs);
    //app_log()<<"evaluate() #"<<omp_get_num_threads()<<" rmax="<<aGrid->rmax()<<" size="<<aGrid->size()<< std::endl;
   
  //  return df;
    return v;
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv)
  {
   // mRealType df = Basis.df_dr(r, gcoefs);
     //app_log()<<"evaluate() #"<<omp_get_thread_num()<<" rmax="<<aGrid->rmax()<<" size="<<aGrid->size()<< std::endl;
 //    return df; 
//    std::stringstream wee;
//    wee<<"srDf() #"<<omp_get_thread_num()<<" dspl= "<<rinv*rinv*du-df<<" ref= "<<df<<" r= "<<r<< std::endl;
//   app_log()<<wee.str();  
    return drV_force->splint(r)/mRealType(r*r) ; 
  }

  inline mRealType srDf_strain(mRealType r, mRealType rinv)
  {
  //  mRealType df = Basis.df_dr(r, gstraincoefs);
  //  return df;
    
    mRealType du=drV_stress->splint(r);
    return rinv*rinv*du; 
  }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r)
  {
    mRealType v=0.0;
   
  // for(int n=0; n<coefs.size(); n++)
  //    v -= coefs[n]*Basis.h(n,r);
    return v;
  }

  inline mRealType evaluateSR_k0()
  {
    mRealType v0=0.0;
    for(int n=0; n<coefs.size(); n++)
      v0 += coefs[n]*Basis.hintr2(n);
    return v0*2.0*TWOPI/Basis.get_CellVolume();
  }


  inline mRealType evaluateLR_r0()
  {
	//this is because the constraint v(r)=sigma(r) as r-->0.  
	// so v(r)-sigma(r)="0".  Divergence prevents me from coding this.
    mRealType v0=0.0;
//    for(int n=0; n<coefs.size(); n++)
//      v0 += coefs[n]*Basis.h(n,0.0);
    return v0;
  }
  
    //This returns the stress derivative of Fk, except for the explicit volume dependence.  The explicit volume dependence is factored away into V.
  inline SymTensor<mRealType, OHMMS_DIM> evaluateLR_dstrain(TinyVector<mRealType, OHMMS_DIM> k, mRealType kmag)
  {
	  SymTensor<mRealType, OHMMS_DIM> deriv_tensor = 0;
	 // mRealType derivconst = Basis.fk(kmag, dcoefs);
	//  app_log()<<"squoo "<<derivconst<< std::endl;
	  
	  for (int dim1=0; dim1<OHMMS_DIM; dim1++)
		for(int dim2=dim1; dim2<OHMMS_DIM; dim2++)
		{
		  mRealType v=0.0;
          deriv_tensor(dim1,dim2)=- evaldYkgstrain(kmag)*k[dim1]*k[dim2]/kmag; //- evaldFk_dk(kmag)*k[dim1]*k[dim2]/kmag ;
          
          if (dim1==dim2) deriv_tensor(dim1,dim2)-= evalYkgstrain(kmag); //+ derivconst;
         // app_log()<<"squoo "<<Basis.fk(kmag, dcoefs(dim1,dim2))<< std::endl;
		}
	  	
		
	  return deriv_tensor;
  }
  
  
  inline SymTensor<mRealType, OHMMS_DIM> evaluateSR_dstrain(TinyVector<mRealType, OHMMS_DIM> r, mRealType rmag)
  {
    SymTensor<mRealType, OHMMS_DIM> deriv_tensor=0;

    mRealType Sr_r=srDf_strain(rmag, 1.0/mRealType(rmag))/mRealType(rmag);

    for (int dim1=0; dim1<OHMMS_DIM; dim1++)
    {
		for(int dim2=dim1; dim2<OHMMS_DIM; dim2++)
		{
	       mRealType v=0.0;

	       deriv_tensor(dim1,dim2)=r[dim1]*r[dim2]*Sr_r;

	    }
	}

     	
	return deriv_tensor;
  }
  
/*  inline mRealType evaluateSR_k0_dstrain()
  {
    mRealType v0=0.0;
    mRealType norm=2.0*TWOPI/Basis.get_CellVolume();
   
    for(int n=0; n<coefs.size(); n++)
      v0 += gstraincoefs[n]*Basis.hintr2(n);
    
    v0*=-norm
    SymTensor stress(v0,0.0,v0,0.0,0.0,v0);
    return stress;
 }
   */
  inline SymTensor<mRealType, OHMMS_DIM> evaluateSR_k0_dstrain()
  {
    mRealType v0=0.0;
    mRealType norm=2.0*TWOPI/Basis.get_CellVolume();
   
    for(int n=0; n<coefs.size(); n++)
      v0 += gstraincoefs[n]*Basis.hintr2(n);
    
    v0*=-norm;
    SymTensor<mRealType, OHMMS_DIM> stress;
    for (int i=0; i<OHMMS_DIM; i++) stress(i,i)=v0;
    
    return stress;
  }
  
  inline mRealType evaluateLR_r0_dstrain(int i, int j)
  {
	//the t derivative for the relevant basis elements are all zero because of constraints.
    return 0.0; //Basis.f(0,dcoefs(i,j));
  }
  
  inline SymTensor<mRealType, OHMMS_DIM> evaluateLR_r0_dstrain()
  {
    SymTensor<mRealType, OHMMS_DIM> stress;
	return stress;
  }

private:


  inline mRealType evalYk(mRealType k)
  {
    //FatK = 4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
    mRealType FatK=myFunc.Vk(k)- Basis.fk(k,coefs);
  //  for(int n=0; n<Basis.NumBasisElem(); n++)
  //    FatK -= coefs[n]*Basis.c(n,k);
    return FatK;
  }
  inline mRealType evalYkg(mRealType k)
  {
    mRealType FatK=myFunc.Vk(k) - Basis.fk(k,gcoefs);
    //for(int n=0; n<Basis.NumBasisElem(); n++)
   //   FatK -= gcoefs[n]*Basis.c(n,k);
    return FatK;
    
  }
  inline mRealType evalYkgstrain(mRealType k)
  {
    mRealType FatK=myFunc.Vk(k) - Basis.fk(k,gstraincoefs);
    //for(int n=0; n<Basis.NumBasisElem(); n++)
   //   FatK -= gcoefs[n]*Basis.c(n,k);
    return FatK;
    
  }
  
  inline mRealType evaldYkgstrain(mRealType k)
  {
    mRealType dFk_dk = myFunc.dVk_dk(k) - Basis.dfk_dk(k,gstraincoefs);
  //  mRealType dFk_dk = myFunc.dVk_dk(k,Basis.get_rc()) - Basis.dfk_dk(k,coefs);
    return dFk_dk;
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
    mRealType chisqr_f=0.0;
    mRealType chisqr_df=0.0;
    mRealType chisqr_strain=0.0; 
    //First we send the new Lattice to the Basis, in case it has been updated.
    Basis.set_Lattice(ref);
    //Compute RC from box-size - in constructor?
    //No here...need update if box changes
    int NumKnots(17);
    Basis.set_NumKnots(NumKnots);
    Basis.set_rc(ref.LR_rc);
    //Initialise the breakup - pass in basis.
    LRBreakup<BreakupBasis> breakuphandler(Basis);
    //Find size of basis from cutoffs
    mRealType kc = (LR_kc<0)?ref.LR_kc:LR_kc;
    //mRealType kc(ref.LR_kc); //User cutoff parameter...
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
	  app_log() <<"\nPerforming Optimized Breakup with Short Range Coulomb Basis\n";	
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
    fillVk(breakuphandler.KList);
    //Allocate the space for the coefficients.
    int Nbasis=Basis.NumBasisElem();
    coefs.resize(Nbasis); //This must be after SetupKVecs.
    gcoefs.resize(Nbasis);
    gstraincoefs.resize(Nbasis);

    //Going to implement a smooth real space cutoff.  This means that alpha=0,1,2 for the LPQHI basis at knot r_c
    //all equal the 0, 1st, and 2nd derivatives of our bare function.  
    //These three functions are the last three basis elements in our set.  


   
    Vector<mRealType> constraints;
    //Vector<mRealType> strainconstraints;

    constraints.resize(Nbasis);
  //  strainconstraints.resize(Nbasis);
    for (int i=0; i < Nbasis; i++) constraints[i]=1;
    

   // mRealType rc=Basis.get_rc();
    
    ///This is to make sure there's no cusp in the LR part.  
    gstraincoefs[0]=gcoefs[0]=coefs[0] = 1.0;
    constraints[0]=0;
   
    gstraincoefs[1]=gcoefs[1] = coefs[1] = 0.0;
    constraints[1]=0;
    
    gstraincoefs[2]=gcoefs[2] = coefs[2] = 0.0; 
    constraints[2]=0.0;
   

    gstraincoefs[Nbasis-1]= gcoefs[Nbasis-1]=coefs[Nbasis-1]=0.0;
    constraints[Nbasis-1]=0;
   
    //1st derivative
    
    gstraincoefs[Nbasis-2]=gcoefs[Nbasis-2]=coefs[Nbasis-2]=0.0;
    constraints[Nbasis-2]=0;

    //Function value 
    gstraincoefs[Nbasis-3]=gcoefs[Nbasis-3]=coefs[Nbasis-3]=0.0;
    constraints[Nbasis-3]=0;
    //And now to impose the constraints
    


    Vector<mRealType> chisqr(3);
  //  chisqr_f=breakuphandler.DoBreakup(Fk.data(),coefs.data(),constraints.data()); //Fill array of coefficients.
  //  chisqr_df=breakuphandler.DoGradBreakup(Fkg.data(), gcoefs.data(), constraints.data());
  //  chisqr_strain=breakuphandler.DoStrainBreakup(Fk.data(), Fkgstrain.data(), gstraincoefs.data(), strainconstraints.data());   
    breakuphandler.DoAllBreakup(chisqr.data(), Fk.data(), Fkgstrain.data(), coefs.data(), gcoefs.data(), gstraincoefs.data(), constraints.data());
    app_log()<<"         LR function chi^2 = "<<chisqr[0]<< std::endl;
    app_log()<<"    LR grad function chi^2 = "<<chisqr[1]<< std::endl;
    app_log()<<"  LR strain function chi^2 = "<<chisqr[2]<< std::endl;
   // app_log()<<"  n  tn   gtn h(n)\n";
     
  //  myFunc.reset(ref);
  //  fillYk(ref.SK->KLists);
  //  fillYkg(ref.SK->KLists);
  //  filldFk_dk(ref.SK->KLists);

    makeSplines(10001);
  }



  void fillVk(std::vector<TinyVector<mRealType,2> >& KList)
  {
    Fk.resize(KList.size());
    Fkg.resize(KList.size());
    Fkgstrain.resize(KList.size());
   // Fk_copy.resize(KList.size());
    for(int ki=0; ki<KList.size(); ki++)
    {
      mRealType k=KList[ki][0];
      Fk[ki] = myFunc.Vk(k); //Call derived fn.
      Fkg[ki]= myFunc.Vk(k);
      Fkgstrain[ki] = myFunc.dVk_dk(k);
     // Fk_copy[ki]=myFunc.Vk(k);
    }
  }

  void fillYk(KContainer& KList)
  {
    Fk.resize(KList.kpts_cart.size());
    const std::vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    Fk_symm.resize(MaxKshell);
    for(int ks=0,ki=0; ks<Fk_symm.size(); ks++)
    {
      mRealType uk=evalYk(std::sqrt(KList.ksq[ki]));
      Fk_symm[ks]=uk;
      while(ki<KList.kshell[ks+1] && ki<Fk.size())
        Fk[ki++]=uk;
    }
    //for(int ki=0; ki<KList.kpts_cart.size(); ki++){
    //  mRealType k=dot(KList.kpts_cart[ki],KList.kpts_cart[ki]);
    //  k=std::sqrt(k);
    //  Fk[ki] = evalFk(k); //Call derived fn.
    //}
  }
  void fillYkg(KContainer& KList)
  {
    Fkg.resize(KList.kpts_cart.size());
    const std::vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;

    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType uk=evalYkg(std::sqrt(KList.ksq[ki]));

      while(ki<KList.kshell[ks+1] && ki<Fkg.size())
        Fkg[ki++]=uk;
    }
  }
  
  void fillYkgstrain(KContainer& KList)
  {
    Fkgstrain.resize(KList.kpts_cart.size());
    const std::vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      mRealType uk=evalYkgstrain(std::sqrt(KList.ksq[ki]));
      while(ki<KList.kshell[ks+1] && ki<Fkgstrain.size())
        Fkgstrain[ki++]=uk;
    }
  }
  void filldFk_dk(KContainer& KList)
  {
    dFk_dstrain.resize(KList.kpts_cart.size());
    

    for (int ki=0; ki<dFk_dstrain.size(); ki++)
    {
	  dFk_dstrain[ki] = evaluateLR_dstrain(KList.kpts_cart[ki], std::sqrt(KList.ksq[ki]));
	}

  }
};
}
#endif
