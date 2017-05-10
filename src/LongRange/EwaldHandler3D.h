//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_EWALD_HANLDER3D_H
#define QMCPLUSPLUS_EWALD_HANLDER3D_H

#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{

/* LR breakup for the standard Ewald method
 *
 * Quasi-2D Ewald method : J. Phys.: Condens. Matter 16, 891 (2004)
 * http://iopscience.iop.org/0953-8984/16/6/017/
 * Note that \f$ \simga \rightarrow 1/\sigma\f$
 * It is possible to use 3D Ewald but for the bulk system, the optimal breakup method
 * is used.
 */
class EwaldHandler3D: public LRHandlerBase
{

public:
  ///type of supercell
  int SuperCellEnum;
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  RealType Sigma;
  ///Volume of the supercell
  RealType Volume;
  ///Area of the supercell: always z is the slab direction
  RealType Area;
  /** Define prefactors for the mixed boundary conditions
   *
   * For quasi-2D (see Appendix A of JPC)
   * PreFactors[0] = \f$ \frac{2\pi}{A}\f$
   * PreFactors[1] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[2] = \f$ \frac{2\pi}{A}\frac{1}{\sigma\pi}\f$
   * PreFactors[3] = \f$ 2\frac{\sqrt{\pi}}{A*\sigma}-\frac{2\pi}{k*A} erfc(\frac{k}{2\sigma}\f$
   */
  TinyVector<RealType,4> PreFactors;
  ///store |k|
  vector<RealType> kMag;
  /// Constructor
  EwaldHandler3D(ParticleSet& ref, RealType kc_in=-1.0)
    : LRHandlerBase(kc_in)
  {
    Sigma=LR_kc=ref.Lattice.LR_kc;
  }

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  EwaldHandler3D(const EwaldHandler3D& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref)
  {
    return new EwaldHandler3D(*this,ref);
  }

  void initBreakup(ParticleSet& ref);

  void Breakup(ParticleSet& ref, RealType rs_in)
  {
    initBreakup(ref);
  }

  void resetTargetParticleSet(ParticleSet& ref) { }

  inline RealType evaluate(RealType r, RealType rinv)
  {
    return erfc(r*Sigma)*rinv;
  }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline RealType evaluateLR(RealType r)
  {
    return erf(r*Sigma)/r;
  }

  inline RealType evaluateSR_k0()
  {
    RealType v0=M_PI/Sigma/Sigma/Volume;
    return v0;
  }

  inline RealType evaluateLR_r0()
  {
    return 2.0*Sigma/std::sqrt(M_PI);
  }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline RealType srDf(RealType r, RealType rinv)
  {
    return -2.0*Sigma*std::exp(-Sigma*Sigma*r*r)/(std::sqrt(M_PI)*r) - erfc(Sigma*r)*rinv*rinv;
  }

  void fillFk(KContainer& KList);

  void fillYkgstrain(KContainer& KList)
  {
    Fkgstrain.resize(KList.kpts_cart.size());
    const vector<int>& kshell(KList.kshell);
    if(MaxKshell >= kshell.size())
      MaxKshell=kshell.size()-1;
    for(int ks=0,ki=0; ks<MaxKshell; ks++)
    {
      RealType uk=evalYkgstrain(std::sqrt(KList.ksq[ki]));
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
  
    //This returns the stress derivative of Fk, except for the explicit volume dependence.  The explicit volume dependence is factored away into V.
  inline SymTensor<RealType, OHMMS_DIM> evaluateLR_dstrain(TinyVector<RealType, OHMMS_DIM> k, RealType kmag)
  {
	  SymTensor<RealType, OHMMS_DIM> deriv_tensor = 0;
	 // RealType derivconst = Basis.fk(kmag, dcoefs);
	//  app_log()<<"squoo "<<derivconst<<endl;
	  
	  for (int dim1=0; dim1<OHMMS_DIM; dim1++)
		for(int dim2=dim1; dim2<OHMMS_DIM; dim2++)
		{
		  RealType v=0.0;
          deriv_tensor(dim1,dim2)=- evaldYkgstrain(kmag)*k[dim1]*k[dim2]/kmag; //- evaldFk_dk(kmag)*k[dim1]*k[dim2]/kmag ;
          
          if (dim1==dim2) deriv_tensor(dim1,dim2)-= evalYkgstrain(kmag); //+ derivconst;
         // app_log()<<"squoo "<<Basis.fk(kmag, dcoefs(dim1,dim2))<<endl;
		}
	  	
		
	  return deriv_tensor;
  }
  
  
  inline SymTensor<RealType, OHMMS_DIM> evaluateSR_dstrain(TinyVector<RealType, OHMMS_DIM> r, RealType rmag)
  {
    SymTensor<RealType, OHMMS_DIM> deriv_tensor=0;

    RealType Sr_r=srDf(rmag, 1.0/RealType(rmag))/RealType(rmag);

    for (int dim1=0; dim1<OHMMS_DIM; dim1++)
    {
		for(int dim2=dim1; dim2<OHMMS_DIM; dim2++)
		{
	       RealType v=0.0;

	       deriv_tensor(dim1,dim2)=r[dim1]*r[dim2]*Sr_r;

	    }
	}

     	
	return deriv_tensor;
  }
  
/*  inline RealType evaluateSR_k0_dstrain()
  {
    RealType v0=0.0;
    RealType norm=2.0*TWOPI/Basis.get_CellVolume();
   
    for(int n=0; n<coefs.size(); n++)
      v0 += gstraincoefs[n]*Basis.hintr2(n);
    
    v0*=-norm
    SymTensor stress(v0,0.0,v0,0.0,0.0,v0);
    return stress;
 }
   */
  inline SymTensor<RealType, OHMMS_DIM> evaluateSR_k0_dstrain()
  {
    RealType v0=-M_PI/Sigma/Sigma/Volume;
    SymTensor<RealType, OHMMS_DIM> stress;
    for (int i=0; i<OHMMS_DIM; i++) stress(i,i)=v0;
    
    return stress;
  }
  
  inline RealType evaluateLR_r0_dstrain(int i, int j)
  {
    return 0.0; 
  }
  
  inline SymTensor<RealType, OHMMS_DIM> evaluateLR_r0_dstrain()
  {
    SymTensor<RealType, OHMMS_DIM> stress;
	return stress;
  }

private:

  inline RealType evalYkgstrain(RealType k)
  {
    RealType norm=4.0*M_PI/Volume;
    RealType denom=4.0*Sigma*Sigma;
    RealType k2=k*k;
    return norm*std::exp(-k2/denom)/k2; 
    
  }
  
  inline RealType evaldYkgstrain(RealType k)
  {
    RealType norm=4.0*M_PI/Volume;
    RealType denom=4.0*Sigma*Sigma;
    RealType sigma2=Sigma*Sigma;
    RealType k2=k*k;
    return -norm*std::exp(-k2/denom)*(denom+k2)/(k*k2*2.0*sigma2);
    
  }

};
}
#endif
