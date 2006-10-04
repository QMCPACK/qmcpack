#ifndef QMCPLUSPLUS_LRCOULOMBAA_H
#define QMCPLUSPLUS_LRCOULOMBAA_H

#include "LongRange/LRHandler.h"
#include "Particle/ParticleSet.h"
#include "LongRange/StructFact.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

  template<class T=double>
  struct CoulombFunctor {
    T NormFactor;
    inline CoulombFunctor(){}
    void reset(T volume) {
      NormFactor=4.0*M_PI/volume;
    }
    inline T operator()(T r, T rinv) { return rinv;}
    inline T Fk(T k, T rc) {
      return NormFactor/(k*k)* std::cos(k*rc);
    }
    inline T Xk(T k, T rc) {
      return -NormFactor/(k*k)* std::cos(k*rc);
    }
  };

  template<class Func, class BreakupBasis=LPQHIBasis>
  class LRHandlerTemp: public LRHandler<BreakupBasis> {
  private:
    //Typedef for the lattice-type.
    typedef ParticleSet::ParticleLayout_t ParticleLayout_t;

    //Import types from the base-class
    typedef LRHandler<BreakupBasis> base_type;
    typedef typename base_type::RealType RealType;
    //Import members from the base-class
    using base_type::Basis;
    using base_type::coefs;
    using base_type::Fk;

    //Private members
    const StructFact* PtclRhoK;
    Func myFunc;

  public:

    //Constructor
    LRHandlerTemp(ParticleSet& ref): 
      //non-default base-class construct. We use 1 function for each species.
      LRHandler<BreakupBasis>(ref.Lattice), PtclRhoK(ref.SK) {
      myFunc.reset(ref.Lattice.Volume);
      LRHandler<BreakupBasis>::InitBreakup(ref.Lattice,1); 
      LRHandler<BreakupBasis>::fillFk(PtclRhoK.KLists);
    }
      
    RealType evalLR() {return 0.0;}
    RealType evalSR() {return 0.0;}
    RealType evalConsts() {return 0.0;}

    void resetTargetParticleSet(ParticleSet& ref) {
      PtclRhok=ref.SK;
      myFunc.reset(ref.Lattice.Volume);
    }

    inline RealType evaluate(RealType r, RealType rinv) {
      RealType v=myFunc(r,rinv);
      for(int n=0; n<coefs.size(); n++) v -= coefs[0][n]*Basis.h(n,sep);
      return v;
    }

    inline RealType evaluate(const vector<int>& minusk, 
        const ComplexType* restrict rk1, const ComplexType* restrict rk2) {
      RealType vk=0.0;
      for(int ki=0; ki<PtclRhoK.KLists.kpts_cart.size(); ki++) {
	vk += (rk1[ki]*rk2[minusk[ki]]).real()*Fk[0][ki];
      } //ki
      return vk;
    }

  private:

    inline RealType evalFk(RealType k,int FunctionIndex) {
      //FatK = 4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
      RealType FatK=myFunc.Fk(k,Basis.get_rc());
      for(int n=0; n<Basis.NumBasisElem(); n++)
        FatK += coefs[0][n]*Basis.c(n,k);
      return FatK;
    }
    inline RealType evalXk(RealType k,int FunctionIndex) {
      //RealType FatK;
      //FatK = -4.0*M_PI/(Basis.get_CellVolume()*k*k)* std::cos(k*Basis.get_rc());
      //return (FatK);
      return myFunc.Xk(k,Basis.get_rc());
    }

  };
}
#endif
