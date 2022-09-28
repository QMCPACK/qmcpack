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
#include "PWOrbitalSet.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
PWOrbitalSet::~PWOrbitalSet()
{
  if (OwnBasisSet && myBasisSet)
    delete myBasisSet;
  if (!IsCloned && C != nullptr)
    delete C;
}

std::unique_ptr<SPOSet> PWOrbitalSet::makeClone() const
{
  auto myclone        = std::make_unique<PWOrbitalSet>(*this);
  myclone->myBasisSet = new PWBasis(*myBasisSet);
  myclone->IsCloned   = true;
  return myclone;
}

void PWOrbitalSet::setOrbitalSetSize(int norbs) {}

void PWOrbitalSet::resize(PWBasisPtr bset, int nbands, bool cleanup)
{
  myBasisSet     = bset;
  OrbitalSetSize = nbands;
  OwnBasisSet    = cleanup;
  BasisSetSize   = myBasisSet->NumPlaneWaves;
  C              = new ValueMatrix(OrbitalSetSize, BasisSetSize);
  Temp.resize(OrbitalSetSize, PW_MAXINDEX);
  app_log() << "  PWOrbitalSet::resize OrbitalSetSize =" << OrbitalSetSize << " BasisSetSize = " << BasisSetSize
            << std::endl;
}

void PWOrbitalSet::addVector(const std::vector<ComplexType>& coefs, int jorb)
{
  int ng = myBasisSet->inputmap.size();
  if (ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int>& inputmap(myBasisSet->inputmap);
  for (int ig = 0; ig < ng; ig++)
  {
    if (inputmap[ig] > -1)
      (*C)[jorb][inputmap[ig]] = coefs[ig];
  }
}

void PWOrbitalSet::addVector(const std::vector<RealType>& coefs, int jorb)
{
  int ng = myBasisSet->inputmap.size();
  if (ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int>& inputmap(myBasisSet->inputmap);
  for (int ig = 0; ig < ng; ig++)
  {
    if (inputmap[ig] > -1)
      (*C)[jorb][inputmap[ig]] = coefs[ig];
  }
}

void PWOrbitalSet::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  //Evaluate every orbital for particle iat.
  //Evaluate the basis-set at these coordinates:
  //myBasisSet->evaluate(P,iat);
  myBasisSet->evaluate(P.activeR(iat));
  MatrixOperators::product(*C, myBasisSet->Zv, psi);
}

void PWOrbitalSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  //Evaluate the orbitals and derivatives for particle iat only.
  myBasisSet->evaluateAll(P, iat);
  MatrixOperators::product(*C, myBasisSet->Z, Temp);
  const ValueType* restrict tptr = Temp.data();
  for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
  {
    psi[j]   = tptr[PW_VALUE];
    d2psi[j] = tptr[PW_LAP];
    dpsi[j]  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
  }
}

void PWOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                        int first,
                                        int last,
                                        ValueMatrix& logdet,
                                        GradMatrix& dlogdet,
                                        ValueMatrix& d2logdet)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    myBasisSet->evaluateAll(P, iat);
    MatrixOperators::product(*C, myBasisSet->Z, Temp);
    const ValueType* restrict tptr = Temp.data();
    for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
    {
      logdet(i, j)   = tptr[PW_VALUE];
      d2logdet(i, j) = tptr[PW_LAP];
      dlogdet(i, j)  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
    }
  }
}
} // namespace qmcplusplus
