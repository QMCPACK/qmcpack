//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Message/Communicate.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
PWOrbitalSet::~PWOrbitalSet()
{
#if !defined(ENABLE_SMARTPOINTER)
  if(OwnBasisSet&&myBasisSet)
    delete myBasisSet;
#endif
}

SPOSetBase* PWOrbitalSet::makeClone() const
{
  PWOrbitalSet *myclone=new PWOrbitalSet(*this);
  myclone->myBasisSet = new PWBasis(*myBasisSet);
  myclone->IsCloned=true;
  return myclone;
}

void PWOrbitalSet::resetParameters(const opt_variables_type& optVariables)
{
  //DO NOTHING FOR NOW
}

void PWOrbitalSet::setOrbitalSetSize(int norbs)
{
}

void PWOrbitalSet::resetTargetParticleSet(ParticleSet& P)
{
  // Not sure what to do here, if anything
  //app_error() << "PWOrbitalSet::resetTargetParticleSet not yet coded." << std::endl;
  //OHMMS::Controller->abort();
}

void PWOrbitalSet::resize(PWBasisPtr bset, int nbands, bool cleanup)
{
  myBasisSet=bset;
  OrbitalSetSize=nbands;
  OwnBasisSet=cleanup;
  BasisSetSize=myBasisSet->NumPlaneWaves;
  //C.resize(OrbitalSetSize,BasisSetSize);
  setIdentity(false);
  Temp.resize(OrbitalSetSize,PW_MAXINDEX);
  app_log() << "  PWOrbitalSet::resize OrbitalSetSize =" << OrbitalSetSize << " BasisSetSize = " << BasisSetSize << std::endl;
}

void PWOrbitalSet::addVector(const std::vector<ComplexType>& coefs,int jorb)
{
  int ng=myBasisSet->inputmap.size();
  if(ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int> &inputmap(myBasisSet->inputmap);
  for(int ig=0; ig<ng; ig++)
  {
    if(inputmap[ig]>-1)
      (*C)[jorb][inputmap[ig]]=coefs[ig];
  }
}

void PWOrbitalSet::addVector(const std::vector<RealType>& coefs,int jorb)
{
  int ng=myBasisSet->inputmap.size();
  if(ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int> &inputmap(myBasisSet->inputmap);
  for(int ig=0; ig<ng; ig++)
  {
    if(inputmap[ig]>-1)
      (*C)[jorb][inputmap[ig]]=coefs[ig];
  }
}

void
PWOrbitalSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  //Evaluate every orbital for particle iat.
  //Evaluate the basis-set at these coordinates:
  //myBasisSet->evaluate(P,iat);
  myBasisSet->evaluate(P.activeR(iat));
  MatrixOperators::product(*C,myBasisSet->Zv,&psi[0]);
}

void
PWOrbitalSet::evaluate(const ParticleSet& P, int iat,
                       ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  //Evaluate the orbitals and derivatives for particle iat only.
  myBasisSet->evaluateAll(P,iat);
  MatrixOperators::product(*C,myBasisSet->Z,Temp);
  const ValueType* restrict tptr=Temp.data();
  for(int j=0; j< OrbitalSetSize; j++, tptr+=PW_MAXINDEX)
  {
    psi[j]    =tptr[PW_VALUE];
    d2psi[j]  =tptr[PW_LAP];
    dpsi[j]=GradType(tptr[PW_GRADX],tptr[PW_GRADY],tptr[PW_GRADZ]);
  }
}

void
PWOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
                                   ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    myBasisSet->evaluateAll(P,iat);
    MatrixOperators::product(*C,myBasisSet->Z,Temp);
    const ValueType* restrict tptr=Temp.data();
    for(int j=0; j< OrbitalSetSize; j++,tptr+=PW_MAXINDEX)
    {
      logdet(i,j)= tptr[PW_VALUE];
      d2logdet(i,j)= tptr[PW_LAP];
      dlogdet(i,j)=GradType(tptr[PW_GRADX],tptr[PW_GRADY],tptr[PW_GRADZ]);
    }
  }
}
}
