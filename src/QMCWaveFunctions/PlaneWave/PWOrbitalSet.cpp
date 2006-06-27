//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Message/Communicate.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {


  void 
  PWOrbitalSet::evaluate(const ParticleSet& P, int iat, std::vector<ValueType>& psi) {
    //Evaluate every orbital for particle iat.
    //Evaluate the basis-set at these coordinates:
    BasisSet->evaluate(P,iat);
    MatrixOperators::product(Coefs,BasisSet->Zv,&psi[0]);
  }

  void 
  PWOrbitalSet::evaluate(const ParticleSet& P, int iat, std::vector<ValueType>& psi, std::vector<GradType>& dpsi, 
      std::vector<ValueType>& d2psi) {
    //Evaluate the orbitals and derivatives for particle iat only.
    BasisSet->evaluateAll(P,iat);
    MatrixOperators::product(Coefs,BasisSet->Z,Temp);
    const ValueType* restrict tptr=Temp.data();
    for(int j=0; j< NumBands; j++, tptr+=PW_MAXINDEX) {
      psi[j]    =tptr[PW_VALUE];
      d2psi[j]  =tptr[PW_LAP];
      dpsi[j]=GradType(tptr[PW_GRADX],tptr[PW_GRADY],tptr[PW_GRADZ]);
    }
  }
    
  void 
  PWOrbitalSet::evaluate(const ParticleSet& P, int first, int last,
      Matrix<ValueType>& logdet, Matrix<GradType>& dlogdet, Matrix<ValueType>& d2logdet) {
    for(int iat=first,i=0; iat<last; iat++,i++){
      BasisSet->evaluateAll(P,iat);
      MatrixOperators::product(Coefs,BasisSet->Z,Temp);
      const ValueType* restrict tptr=Temp.data();
      for(int j=0; j< NumBands; j++,tptr+=PW_MAXINDEX) {
        logdet(j,i)= tptr[PW_VALUE];
        d2logdet(i,j)= tptr[PW_LAP];
        dlogdet(i,j)=GradType(tptr[PW_GRADX],tptr[PW_GRADY],tptr[PW_GRADZ]);
      }
    }
  }
}
