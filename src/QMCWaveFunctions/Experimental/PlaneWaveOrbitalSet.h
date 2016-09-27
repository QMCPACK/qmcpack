//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_H

#include "QMCWaveFunctions/PlaneWaveBasis.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{

class PlaneWaveOrbitalSet: public OhmmsElementBase
{
private:
  ParticleSet* ptargetPtcl;
public:
  typedef PlaneWaveBasis BasisSet_t;
  //Typedefs from PlaneWaveBasis's parent
  typedef PlaneWaveBasis::RealType  RealType;
  typedef PlaneWaveBasis::ValueType ValueType;

  //My basis set
  BasisSet_t* BasisSet;

  //Wavefunction parameters
  int numspins, numbands;

  PlaneWaveOrbitalSet(ParticleSet* ref, int nspins,int nbands,
                      TinyVector<RealType,3> twist_angle):
    ptargetPtcl(ref), numbands(nbands), numspins(nspins)
  {
    //Allocate the basisset.
    BasisSet = new PlaneWaveBasis(twist_angle,ptargetPtcl->getTotalNum());
  }

  ~PlaneWaveOrbitalSet()
  {
    //Clean up the basis set
    delete BasisSet;
  }

  inline void resizeCoefs()
  {
    //Resize Matrix based on number of bands and number of PWs
    //Need +1 to discard elements of basis functions outside ecut
    Coefs.resize(ptargetPtcl->getTotalNum(),BasisSet->NumPlaneWaves+1);
  }

  inline void addVector(std::vector<ValueType> coefs,int jorb)
  {
    //Jorb is the orbital index.
    //This is the same as the grouping order in ParticleSet.
    assert(jorb < ptargetPtcl->getTotalNum());
    for(int ig=0; ig<coefs.size(); ig++)
      Coefs(jorb,BasisSet->inputmap[ig]) = coefs[ig];
  }

  ///return the number of single particle orbitals
  inline int numOrbitals() const
  {
    return Coefs.rows();
  }

  ///return the number of basis functions
  inline int numBasis() const
  {
    return BasisSet->NumPlaneWaves;
  }

  inline void reset() { }

  ///put from XML
  bool put(xmlNodePtr cur)
  {
    std::cout << "PUT FROM XML CALLED: NOT VALID HERE" << std::endl;
    OHMMS::Controller->abort();
    return true;
  }

  ///write to a std::ostream
  bool get(std::ostream& ) const
  {
    std::cout << "GET TO OSTREAM CALLED. Not yet coded." << std::endl;
    std::exit(0);
    return true;
  }

  ///read from std::istream
  bool put(std::istream& )
  {
    std::cout << "PUT FROM ISTREAM CALLED. Not yet coded." << std::endl;
    std::exit(0);
    return true;
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    std::cout << "resetTargetParticleSet not yet coded." << std::endl;
    OHMMS::Controller->abort();
    //      change in basis and here. This involves resizing of storage.
    //      What happens to Coefs in number of particles changes? Must reread from HDF.
    //      LOGMSG("PlaneWaveOrbitalSet::resetTargetParticleSet with " << P.getName())
    //      BasisSet->resetTargetParticleSet(P);
  }

  ///resize the internal storage of BasisSet by the number of particles
  //keep same number of basis elements.
  inline void resize(int nptcl)
  {
    BasisSet->resize(nptcl);
  }

  //Plane-wave coefficients: (iband,g-vector)
  Matrix<ValueType> Coefs;

  inline RealType
  evaluate(const ParticleSet& P, int iat, int jorb)
  {
    LOGMSG("PWOSet: this should not be used");
    OHMMS::Controller->abort();
    return 0.0;
  }

  template<class VV>
  inline void
  evaluate(const ParticleSet& P, int iat, VV& psi)
  {
    //Evaluate every orbital for particle iat.
    //Evaluate the basis-set at these coordinates:
    BasisSet->evaluate(P,iat,iat+1);
    int NumPtcls = P.getTotalNum();
    for(int j=0; j<NumPtcls; j++)
    {
      psi[j] = dot(&Coefs(j,0),BasisSet->y(iat),BasisSet->NumPlaneWaves);
    }
  }

  template<class VV, class GV>
  inline void
  evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    //Evaluate the orbitals and derivatives for particle iat only.
    BasisSet->evaluateAll(P,iat,iat+1);
    TinyVector<ValueType,3> tmpvec;
    int NumPtcls = P.getTotalNum(), BasisSize=BasisSet->NumPlaneWaves;
    for(int j=0 ; j<NumPtcls; j++)
    {
      psi[j] = dot(&Coefs(j,0),BasisSet->y(iat),BasisSize);
      tmpvec = dot(&Coefs(j,0),BasisSet->dy(iat),BasisSize);
      for(int idim=0; idim<3; idim++)
        dpsi[j][idim]  = tmpvec[idim];
      d2psi[j] = dot(&Coefs(j,0),BasisSet->d2y(iat),BasisSize);
    }
  }

  template<class VM, class GM>
  inline void
  evaluate(const ParticleSet& P, int first, int last,
           VM& logdet, GM& dlogdet, VM& d2logdet)
  {
    //Check Ptclset has same size as storage allocated
    if(P.getTotalNum() != ptargetPtcl->getTotalNum())
    {
      LOGMSG("PlaneWaveOrbitalSet Error: Coefs storage not allocated correctly");
      OHMMS::Controller->abort();
    }
    //Evaluate for all particles in a group.
    //Order of "groups" in particleset is same as order of groups in Coefs,
    //so same the first/last is valid.
    //Parameters
    int NumPtcls=last-first, BasisSize = BasisSet->NumPlaneWaves;
    //Evaluate the basis-set at these coordinates:
    BasisSet->evaluateAll(P,first,last);
    TinyVector<ValueType,3> tmpvec;
    //Indices for slater matrices are zero-based. Indices for basis functions and coefficients
    //are "first"-based.
    for(int i=0, iat=first; i<NumPtcls; i++, iat++)
    {
      //Orbital index...
      for(int j=0, jat=first ; j<NumPtcls; j++, jat++)
      {
        logdet(j,i) = dot(&Coefs(jat,0),BasisSet->y(iat),BasisSize);
        tmpvec = dot(&Coefs(jat,0),BasisSet->dy(iat),BasisSize);
        for(int idim=0; idim<3; idim++)
          dlogdet(i,j)[idim] = tmpvec[idim];
        d2logdet(i,j) = dot(&Coefs(jat,0),BasisSet->d2y(iat),BasisSize);
      }
    }
  }
};
}
#endif
