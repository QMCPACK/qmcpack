//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/DeterminantOperators.h"
#include <algorithm>

namespace qmcplusplus {

  class PWOrbitalSet: public QMCTraits {

  public:

    typedef PWBasis BasisSet_t;

    /** inherit the enum of BasisSet_t */
    enum {PW_VALUE=   BasisSet_t::PW_VALUE, 
          PW_LAP=     BasisSet_t::PW_LAP, 
          PW_GRADX=   BasisSet_t::PW_GRADX, 
          PW_GRADY=   BasisSet_t::PW_GRADY, 
          PW_GRADZ=   BasisSet_t::PW_GRADZ,
          PW_MAXINDEX=BasisSet_t::PW_MAXINDEX
    };
    
    /** default constructor
     */
    PWOrbitalSet(): NumBands(0), BasisSize(0), BasisSet(0), OwnBasisSet(false) {
    }

    /** delete BasisSet only it owns this
     *
     * Builder takes care of who owns what
     */
    ~PWOrbitalSet() {
      if(OwnBasisSet&&BasisSet) delete BasisSet;
    }

    void resize(BasisSet_t* bset, int nbands, bool cleanup=false) {
      BasisSet=bset;
      NumBands=nbands;
      OwnBasisSet=cleanup;
      BasisSize=BasisSet->NumPlaneWaves;
      Coefs.resize(NumBands,BasisSize);
      Temp.resize(NumBands,PW_MAXINDEX);
    }


    /** Builder class takes care of the assertion
     */
    inline void addVector(const std::vector<ValueType>& coefs,int jorb) {
      std::copy(coefs.begin(),coefs.end(),Coefs[jorb]);
    }

    ///return the number of single particle orbitals
    inline int numOrbitals() const { return NumBands;}

    ///return the number of basis functions
    inline int numBasis() const { return BasisSize;}

    inline void reset() { }

    void resetTargetParticleSet(ParticleSet& P) {
      cout << "resetTargetParticleSet not yet coded." << endl;
      OHMMS::Controller->abort();
    }

    inline RealType
      evaluate(const ParticleSet& P, int iat, int jorb) {
      LOGMSG("PWOSet: this should not be used");
      OHMMS::Controller->abort();
      return 0.0;
    }

    void 
      evaluate(const ParticleSet& P, int iat, std::vector<ValueType>& psi);

    void 
      evaluate(const ParticleSet& P, int iat, 
          std::vector<ValueType>& psi, 
          std::vector<GradType>& dpsi, 
          std::vector<ValueType>& d2psi);

    void 
      evaluate(const ParticleSet& P, int first, int last,
          Matrix<ValueType>& logdet, Matrix<GradType>& dlogdet, Matrix<ValueType>& d2logdet);
    
    /** boolean
     *
     * If true, this has to delete the BasisSet
     */
    bool OwnBasisSet;
    ///Number of bands of this PWOrbitalSet
    int NumBands;
    ///Size of the basis set
    int BasisSize;
    ///TwistAngle of this PWOrbitalSet
    PosType TwistAngle;
    ///My basis set
    BasisSet_t* BasisSet;
    ///Plane-wave coefficients: (iband,g-vector)
    Matrix<ValueType> Coefs;
    /** temporary array to perform gemm operation */
    Matrix<ValueType> Temp;
    };
}
#endif
