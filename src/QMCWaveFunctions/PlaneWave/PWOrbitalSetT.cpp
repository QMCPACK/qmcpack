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
#include "PWOrbitalSetT.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
template<class T>
PWOrbitalSetT<T>::~PWOrbitalSetT()
{
  if (OwnBasisSet && myBasisSet)
    delete myBasisSet;
  if (!IsCloned && this->C != nullptr)
    delete this->C;
}

template<class T>
std::unique_ptr<SPOSetT<T>> PWOrbitalSetT<T>::makeClone() const
{
  auto myclone        = std::make_unique<PWOrbitalSetT<T>>(*this);
  myclone->myBasisSet = new PWBasisT<T>(*myBasisSet);
  myclone->IsCloned   = true;
  return myclone;
}

template<class T>
void PWOrbitalSetT<T>::setOrbitalSetSize(int norbs) {}

template<class T>
void PWOrbitalSetT<T>::resize(PWBasisPtr bset, int nbands, bool cleanup)
{
  myBasisSet     = bset;
  this->OrbitalSetSize = nbands;
  OwnBasisSet    = cleanup;
  BasisSetSize   = myBasisSet->NumPlaneWaves;
  this->C              = new ValueMatrix(this->OrbitalSetSize, BasisSetSize);
  this->Temp.resize(this->OrbitalSetSize, PW_MAXINDEX);
  app_log() << "  PWOrbitalSetT<T>::resize OrbitalSetSize =" << this->OrbitalSetSize << " BasisSetSize = " << BasisSetSize
            << std::endl;
}

template<class T>
void PWOrbitalSetT<T>::addVector(const std::vector<ComplexType>& coefs, int jorb)
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
      (*(this->C))[jorb][inputmap[ig]] = coefs[ig];
  }
}

template<class T>
void PWOrbitalSetT<T>::addVector(const std::vector<RealType>& coefs, int jorb)
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
      (*(this->C))[jorb][inputmap[ig]] = coefs[ig];
  }
}

template<class T>
void PWOrbitalSetT<T>::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  //Evaluate every orbital for particle iat.
  //Evaluate the basis-set at these coordinates:
  //myBasisSet->evaluate(P,iat);
  myBasisSet->evaluate(P.activeR(iat));
  MatrixOperators::product<T>(*(this->C), myBasisSet->Zv, psi);
}

template<class T>
void PWOrbitalSetT<T>::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  //Evaluate the orbitals and derivatives for particle iat only.
  myBasisSet->evaluateAll(P, iat);
  MatrixOperators::product<T>(*(this->C), myBasisSet->Z, this->Temp);
  const T* restrict tptr = this->Temp.data();
  for (int j = 0; j < this->OrbitalSetSize; j++, tptr += PW_MAXINDEX)
  {
    psi[j]   = tptr[PW_VALUE];
    d2psi[j] = tptr[PW_LAP];
    dpsi[j]  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
  }
}

template<class T>
void PWOrbitalSetT<T>::evaluate_notranspose(const ParticleSet& P,
                                        int first,
                                        int last,
                                        ValueMatrix& logdet,
                                        GradMatrix& dlogdet,
                                        ValueMatrix& d2logdet)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    myBasisSet->evaluateAll(P, iat);
    MatrixOperators::product<T>(*(this->C), myBasisSet->Z, this->Temp);
    const T* restrict tptr = this->Temp.data();
    for (int j = 0; j < this->OrbitalSetSize; j++, tptr += PW_MAXINDEX)
    {
      logdet(i, j)   = tptr[PW_VALUE];
      d2logdet(i, j) = tptr[PW_LAP];
      dlogdet(i, j)  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
    }
  }
}

// Class concrete types from T
// NOTE: This class only gets compiled if QMC_COMPLEX is defined, thus it is inherently complex
// template class PWOrbitalSetT<double>;
// template class PWOrbitalSetT<float>;
template class PWOrbitalSetT<std::complex<double>>;
template class PWOrbitalSetT<std::complex<float>>;
} // namespace qmcplusplus
