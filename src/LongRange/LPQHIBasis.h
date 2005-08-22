#ifndef OHMMS_QMC_LPQHIBASIS_H
#define OHMMS_QMC_LPQHIBASIS_H

#include "LongRange/LRBasis.h"
#include <complex>

namespace ohmmsqmc {

  /** @ingroup longrange
   *\brief A derivative of LRBasis class to provide the functionality 
   * of the LPQHI basis. Based on code by Ken Esler from PIMC++.
   */

  class LPQHIBasis: public LRBasis {
  private:
    int NumKnots; //Number of knots for basis.
    RealType delta, deltainv;
    Matrix<RealType> S; //Coefficients for LPQHI

    //Helper functions for computing FT of basis functions (used in c(n,k))
    inline complex<RealType> Eplus(int i, RealType k, int n);
    inline complex<RealType> Eminus(int i, RealType k, int n);
    inline RealType Dplus(int i, RealType k, int n);
    inline RealType Dminus(int i, RealType k, int n);

    vector<RealType> tvec; //Coefficients
      

  public:
    inline RealType get_delta() { return delta; }
    void set_NumKnots(int n); // n >= 2 required
    void set_rc(RealType rc);
    inline int NumBasisElem() {return 3*NumKnots;}
    RealType h(int n, RealType r);
    RealType hintr2(int n);
    RealType c(int n, RealType k);
    //Constructor...fill S matrix...call correct base-class constructor
    LPQHIBasis(ParticleLayout_t& ref) : LRBasis(ref), NumKnots(0), delta(0.0) {
      S.resize(3,6);
      S(0,0)=1.0;  S(0,1)=0.0; S(0,2)=0.0; S(0,3)=-10.0; 
      S(0,4)=15.0; S(0,5)=-6.0; 

      S(1,0)=0.0; S(1,1)=1.0; S(1,2)=0.0; S(1,3)=-6.0;  
      S(1,4)=8.0; S(1,5)=-3.0; 
      
      S(2,0)=0.0; S(2,1)=0.0; S(2,2)=0.5; S(2,3)=-1.5;  
      S(2,4)=1.5; S(2,5)=-0.5; 
    }
  };

}

#endif
