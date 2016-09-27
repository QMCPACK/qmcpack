//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                     Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LRHANDLER_H
#define QMCPLUSPLUS_LRHANDLER_H

#include "Particle/ParticleSet.h"
#include "LongRange/KContainer.h"
#include "LongRange/LRBasis.h"
#include "LongRange/LRBreakup.h"

//Template-class declaration

namespace qmcplusplus
{

/** @ingroup longrange
 *\brief Base-class for long-range evaluations of functions.
 * Override for each functional form, such as Coulomb and Pade Jastrows.
 */

template<class BreakupBasis>
class LRHandler: public QMCTraits
{
protected:
  //Typedef for the lattice-type.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  //The basis that has been used for the breakup.
  BreakupBasis Basis; //This needs a Lattice for the constructor...

  //This is how many different functions the derived class will decompose
  int NumFns;

  //m_consts: the breakup has 3 parts: short, long & const. Const doesn't
  //change unless cell varies => store as member
  RealType m_consts;

public:
  //Constructor
  LRHandler(ParticleLayout_t& ref) : Basis(ref) //Pass lattice to constructor of Basis.
  {
    /* DO NOT CALL INITBREAKUP IN THIS BASECLASS CONSTRUCTOR!!
       IT CALLS PURE VIRTUAL FUNCTIONS FOLLOWED TO THE DERIVED CLASS...
       Instead, we call it from the derived class constructor.
    */
  }

  //Initialise the basis&breakup.
  //This can be re-called when Lattice is changed.
  //This can be overriden if required.
  virtual void InitBreakup(ParticleLayout_t& ref, int NumFunctions);

  //Override these to perform evaluations of terms.
  virtual RealType evalLR() = 0;
  virtual RealType evalSR() = 0;
  virtual RealType evalConsts() = 0;
  //Override IF needed.
  //Constants are not particle position dependent - compute once only...
  virtual RealType evalTotal()
  {
    RealType LR = evalLR();
    RealType SR = evalSR();
    //cout << "Constant   terms = " << m_consts<< std::endl;
    //cout << "LongRange  terms = " << LR << std::endl;
    //cout << "ShortRange terms = " << SR<< std::endl;
    //cout << "Total            = " << LR+SR+m_consts<< std::endl<< std::endl;
    return (LR+SR+m_consts);
  }

protected:
  //Override this to provide Fk or Xk for a |k|.
  //Fk is the Fourier transform of V_l(r).
  //Xk is the Fourier transform of -V(r) from rc->infinity.
  //Function index is for evaluating multiple functional forms with the same
  //derived class.
  virtual RealType evalFk(RealType k,int FunctionIndex) = 0;
  virtual RealType evalXk(RealType k,int FunctionIndex) = 0;
  //Now we have two functions to fill the array with above values.
  //fillFk needs a full KList which matches (in order and size) the list
  //that will be used from the structure-factor.
  //fillXk requires a 2-dimensional array with index1 = |k|,
  //                                           index2 = degeneracy
  void fillFk(KContainer& KList);
  void fillXk(std::vector<TinyVector<RealType,2> >& KList);

  //These both have 2 indices. The first is the function index, so
  //multiple functions can be broken up by the same derived class (examples:
  // species index on potential, different Jastrow functions...). The
  // second index is the basis index for coefs and the k-vector index for
  // Fk.
  //Breakup coefficients.
  Matrix<RealType> coefs;
  //k-space representation of function for breakup.
  Matrix<RealType> Fk;

};

//Template-class Definitions
template<class BreakupBasis>
void
LRHandler<BreakupBasis>::InitBreakup(ParticleLayout_t& ref,int NumFunctions)
{
  //Here we initialise the basis and coefficients for the long-range
  //beakup. We loocally create a breakup handler and pass in the basis
  //that has been initialised here. We then discard the handler, leaving
  //basis and coefs in a usable state.
  //This method can be re-called later if lattice changes shape.
  NumFns = NumFunctions; //How many functions the derived class will breakup
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
  RealType kc(ref.LR_kc); //User cutoff parameter...
  //kcut is the cutoff for switching to approximate k-point degeneracies for
  //better performance in making the breakup. A good bet is 30*K-spacing so that
  //there are 30 "boxes" in each direction that are treated with exact degeneracies.
  //Assume orthorhombic cell just for deriving this cutoff - should be insensitive.
  //K-Spacing = (kpt_vol)**1/3 = 2*pi/(cellvol**1/3)
  RealType kcut = 60*M_PI*std::pow(Basis.get_CellVolume(),-1.0/3.0);
  //Use 3000/LMax here...==6000/rc for non-ortho cells
  RealType kmax(6000.0/ref.LR_rc);
  breakuphandler.SetupKVecs(kc,kcut,kmax);
  //Set up x_k
  //This is the FT of -V(r) from r_c to infinity.
  //This is the only data that the breakup handler needs to do the breakup.
  //We temporarily store it in Fk, which is replaced with the full FT (0->inf)
  //of V_l(r) after the breakup has been done.
  this->fillXk(breakuphandler.KList);
  //Allocate the space for the coefficients.
  coefs.resize(NumFns,Basis.NumBasisElem()); //This must be after SetupKVecs.
  for(int fn=0; fn<NumFns; fn++)
    breakuphandler.DoBreakup(Fk[fn],coefs[fn]); //Fill array of coefficients.
  //Here we store the constant part of the potential. It only changes
  //when cell-size or shape changes, in which case this function is
  //called again anyway.
  m_consts = this->evalConsts();
}

template<class BreakupBasis>
void
LRHandler<BreakupBasis>::fillXk(std::vector<TinyVector<RealType,2> >& KList)
{
  Fk.resize(NumFns,KList.size());
  for(int fn=0; fn<NumFns; fn++)
    for(int ki=0; ki<KList.size(); ki++)
    {
      RealType k=KList[ki][0];
      Fk[fn][ki] = evalXk(k,fn); //Call derived fn.
    }
}

template<class BreakupBasis>
void
LRHandler<BreakupBasis>::fillFk(KContainer& KList)
{
  Fk.resize(NumFns,KList.kpts_cart.size());
  for(int fn=0; fn<NumFns; fn++)
    for(int ki=0; ki<KList.kpts_cart.size(); ki++)
    {
      RealType k=dot(KList.kpts_cart[ki],KList.kpts_cart[ki]);
      k=std::sqrt(k);
      Fk[fn][ki] = evalFk(k,fn); //Call derived fn.
    }
}
}
#endif
