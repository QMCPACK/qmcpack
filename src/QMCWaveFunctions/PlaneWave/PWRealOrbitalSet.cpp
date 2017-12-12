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
    
    
/** @file PWRealOrbitalSet.cpp
 * @brief declaration of the member functions of PWRealOrbitalSet
 *
 * Not the most optimized method to use wavefunctions in a plane-wave basis.
 */
#include "Message/Communicate.h"
#include "QMCWaveFunctions/PlaneWave/PWRealOrbitalSet.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
PWRealOrbitalSet::~PWRealOrbitalSet()
{
#if !defined(ENABLE_SMARTPOINTER)
  if(OwnBasisSet&&myBasisSet)
    delete myBasisSet;
#endif
}

SPOSetBase* PWRealOrbitalSet::makeClone() const
{
  PWRealOrbitalSet *myclone=new PWRealOrbitalSet(*this);
  myclone->myBasisSet = new PWBasis(*myBasisSet);
  return myclone;
}


void PWRealOrbitalSet::resetParameters(const opt_variables_type& active)
{
  //DO NOTHING FOR NOW
}

void PWRealOrbitalSet::setOrbitalSetSize(int norbs)
{
}

void PWRealOrbitalSet::resetTargetParticleSet(ParticleSet& P)
{
  //  Not sure what to do here, if anything
  //app_error() << "PWRealOrbitalSet::resetTargetParticleSet not yet coded." << std::endl;
  //OHMMS::Controller->abort();
}

void PWRealOrbitalSet::resize(PWBasisPtr bset, int nbands, bool cleanup)
{
  myBasisSet=bset;
  OrbitalSetSize=nbands;
  OwnBasisSet=cleanup;
  BasisSetSize=myBasisSet->NumPlaneWaves;
  CC.resize(OrbitalSetSize,BasisSetSize);
  Temp.resize(OrbitalSetSize,PW_MAXINDEX);
  tempPsi.resize(OrbitalSetSize);
  app_log() << "  PWRealOrbitalSet::resize OrbitalSetSize =" << OrbitalSetSize << " BasisSetSize = " << BasisSetSize << std::endl;
}

void PWRealOrbitalSet::addVector(const std::vector<RealType>& coefs,int jorb)
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
      CC[jorb][inputmap[ig]]=coefs[ig];
  }
}

void PWRealOrbitalSet::addVector(const std::vector<ComplexType>& coefs,int jorb)
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
      CC[jorb][inputmap[ig]]=coefs[ig];
  }
}

void
PWRealOrbitalSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  myBasisSet->evaluate(P.activeR(iat));
  MatrixOperators::product(CC,myBasisSet->Zv,tempPsi.data());
  for(int j=0; j<OrbitalSetSize; j++)
    psi[j]=tempPsi[j].real();
}

void
PWRealOrbitalSet::evaluate(const ParticleSet& P, int iat,
                           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  myBasisSet->evaluateAll(P,iat);
  MatrixOperators::product(CC,myBasisSet->Z,Temp);
  const ComplexType* restrict tptr=Temp.data();
  for(int j=0; j< OrbitalSetSize; j++, tptr+=PW_MAXINDEX)
  {
    psi[j]   =tptr[PW_VALUE].real();
    d2psi[j] =tptr[PW_LAP].real();
#if OHMMS_DIM==3
    dpsi[j]  =GradType(tptr[PW_GRADX].real(),tptr[PW_GRADY].real(),tptr[PW_GRADZ].real());
#elif OHMMS_DIM==2
    dpsi[j]  =GradType(tptr[PW_GRADX].real(),tptr[PW_GRADY].real());
#elif OHMMS_DIM==1
    dpsi[j]  =GradType(tptr[PW_GRADX].real());
#else
#error "Only physical dimensions 1/2/3 are supported."
#endif
  }
}

void
PWRealOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
                                       ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    myBasisSet->evaluateAll(P,iat);
    MatrixOperators::product(CC,myBasisSet->Z,Temp);
    const ComplexType* restrict tptr=Temp.data();
    for(int j=0; j< OrbitalSetSize; j++,tptr+=PW_MAXINDEX)
    {
      convert(tptr[PW_VALUE],logdet(i,j));
      convert(tptr[PW_LAP],d2logdet(i,j));
#if OHMMS_DIM==3
      dlogdet(i,j) = GradType(tptr[PW_GRADX].real(),tptr[PW_GRADY].real(),tptr[PW_GRADZ].real());
#elif OHMMS_DIM==2
      dlogdet(i,j) = GradType(tptr[PW_GRADX].real(),tptr[PW_GRADY].real());
#elif OHMMS_DIM==1
      dlogdet(i,j) = GradType(tptr[PW_GRADX].real());
#else
#error "Only physical dimensions 1/2/3 are supported."
#endif
    }
  }
}
}
