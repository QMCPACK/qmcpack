//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALBUILD_BLAS_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
namespace qmcplusplus {

  /** OrbitalBuilder for Slater determinants in PW basis
  */
  class PWOrbitalBuilder: public OrbitalBuilderBase {

  private:
    ///Number of up and down particles from ParticleSet
    int nup, ndown, upindx;
    ///Index of spin data from HDF5 file to use for updet and downdet
    int updetspinindex, downdetspinindex;

    ///Read routine for HDF wavefunction file version 0.10
    void ReadHDFWavefunction010(hid_t hfile,double& ecut);

    //Storage for the orbitals and basis is created in PWOSet.
    std::map<std::string,PWOrbitalSet*> PWOSet;
    //the last PWOrbitalSet created by the builder
    PWOrbitalSet* myPWOSUp;
    PWOrbitalSet* myPWOSDown;
    //will do something for twist
    PWBasis* myBasisSet;

  public:

    ///constructor
    PWOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);
    ~PWOrbitalBuilder();

    ///implement vritual function
    bool put(xmlNodePtr cur);
  };
}
#endif
