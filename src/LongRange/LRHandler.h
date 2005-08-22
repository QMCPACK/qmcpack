#ifndef OHMMS_QMC_LRHANDLER_H
#define OHMMS_QMC_LRHANDLER_H

#include "Particle/ParticleSet.h"
#include "LongRange/KContainer.h"
#include "LongRange/LRBasis.h"
#include "LongRange/LRBreakup.h"

//Template-class declaration

namespace ohmmsqmc {
  
  /** @ingroup longrange
   *\brief Base-class for long-range evaluations of functions.
   * Override for each functional form, such as Coulomb and Pade Jastrows.
   */

  template<class BreakupBasis>
    class LRHandler: public QMCTraits {
  protected:
    //Typedef for the lattice-type.
    typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
    //The basis that has been used for the breakup.
    BreakupBasis Basis; //This needs a Lattice for the constructor...

    //This is how many different functions the derived class will decompose
    int NumFns;

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
    virtual RealType evalTotal() { return (evalLR()+evalSR()+evalConsts()); }

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
    void fillXk(vector<TinyVector<RealType,2> >& KList);

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

}


using namespace ohmmsqmc;

//Template-class Definitions
template<class BreakupBasis>
void
LRHandler<BreakupBasis>::InitBreakup(ParticleLayout_t& ref,int NumFunctions) {
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
  RealType rc(ref.a(0)[0]*0.5);
  int NumKnots(15);
  Basis.set_NumKnots(NumKnots);
  Basis.set_rc(rc);

  //Initialise the breakup - pass in basis.
  LRBreakup<BreakupBasis> breakuphandler(Basis);

  //Find size of basis from cutoffs
  RealType kc(3.2); //User parameter...
  RealType kcut = max(25.0,kc); //
  RealType kmax(300.0); //Use 3000/L here...
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
}

template<class BreakupBasis>
void
LRHandler<BreakupBasis>::fillXk(vector<TinyVector<RealType,2> >& KList) {
  Fk.resize(NumFns,KList.size());

  for(int fn=0; fn<NumFns; fn++)
    for(int ki=0; ki<KList.size(); ki++) {
      RealType k=KList[ki][0];
      Fk[fn][ki] = evalXk(k,fn); //Call derived fn.
    }
}

template<class BreakupBasis>
void
LRHandler<BreakupBasis>::fillFk(KContainer& KList) {
  Fk.resize(NumFns,KList.kpts_cart.size());

  for(int fn=0; fn<NumFns; fn++)
    for(int ki=0; ki<KList.kpts_cart.size(); ki++){
      RealType k=dot(KList.kpts_cart[ki],KList.kpts_cart[ki]);
      k=sqrt(k);
      Fk[fn][ki] = evalFk(k,fn); //Call derived fn.
    }
}
#endif
