//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CPU/e2iphi.h"
#include "EinsplineSet.h"
#include "einspline/multi_bspline.h"
#include "CPU/math.hpp"
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
template<typename StorageType>
inline void EinsplineSetExtended<StorageType>::computePhaseFactors(const TinyVector<RealType, OHMMS_DIM>& r)
{
  APP_ABORT("EinsplineSetExtended<StorageType>::computePhaseFactors called");
  for (int i = 0; i < kPoints.size(); i++)
    phase[i] = -dot(r, kPoints[i]);
  eval_e2iphi(kPoints.size(), phase.data(), eikr.data());
  //eval_e2iphi(phase,eikr);
  //#ifdef HAVE_MKL
  //    for (int i=0; i<kPoints.size(); i++)
  //      phase[i] = -dot(r, kPoints[i]);
  //    vzCIS(OrbitalSetSize, phase, (double*)eikr.data());
  //#else
  //    double s, c;
  //    for (int i=0; i<kPoints.size(); i++) {
  //      phase[i] = -dot(r, kPoints[i]);
  //      qmcplusplus::sincos (phase[i], &s, &c);
  //      eikr[i] = std::complex<double>(c,s);
  //    }
  //#endif
}


EinsplineSet::UnitCellType EinsplineSet::GetLattice() { return SuperLattice; }

void EinsplineSet::resetSourceParticleSet(ParticleSet& ions) {}

void EinsplineSet::setOrbitalSetSize(int norbs) { OrbitalSetSize = norbs; }

// Real evaluation functions
inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<double>& psi)
{
  eval_multi_UBspline_3d_d(spline, r[0], r[1], r[2], psi.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline, TinyVector<double, 3> r, std::vector<double>& psi)
{
  eval_multi_UBspline_3d_d(spline, r[0], r[1], r[2], &(psi[0]));
}

inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<double>& psi,
                               Vector<TinyVector<double, 3>>& grad)
{
  eval_multi_UBspline_3d_d_vg(spline, r[0], r[1], r[2], psi.data(), (double*)grad.data());
}


inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<double>& psi,
                               Vector<TinyVector<double, 3>>& grad,
                               Vector<Tensor<double, 3>>& hess)
{
  eval_multi_UBspline_3d_d_vgh(spline, r[0], r[1], r[2], psi.data(), (double*)grad.data(), (double*)hess.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_d* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<double>& psi,
                               Vector<TinyVector<double, 3>>& grad,
                               Vector<Tensor<double, 3>>& hess,
                               Vector<TinyVector<Tensor<double, 3>, 3>>& gradhess)
{
  eval_multi_UBspline_3d_d_vghgh(spline, r[0], r[1], r[2], psi.data(), (double*)grad.data(), (double*)hess.data(),
                                 (double*)gradhess.data());
}


//////////////////////////////////
// Complex evaluation functions //
//////////////////////////////////
inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<std::complex<double>>& psi)
{
  eval_multi_UBspline_3d_z(spline, r[0], r[1], r[2], psi.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<std::complex<double>>& psi,
                               Vector<TinyVector<std::complex<double>, 3>>& grad)
{
  eval_multi_UBspline_3d_z_vg(spline, r[0], r[1], r[2], psi.data(), (std::complex<double>*)grad.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<std::complex<double>>& psi,
                               Vector<TinyVector<std::complex<double>, 3>>& grad,
                               Vector<Tensor<std::complex<double>, 3>>& hess)
{
  eval_multi_UBspline_3d_z_vgh(spline, r[0], r[1], r[2], psi.data(), (std::complex<double>*)grad.data(),
                               (std::complex<double>*)hess.data());
}

inline void EinsplineMultiEval(multi_UBspline_3d_z* restrict spline,
                               const TinyVector<double, 3>& r,
                               Vector<std::complex<double>>& psi,
                               Vector<TinyVector<std::complex<double>, 3>>& grad,
                               Vector<Tensor<std::complex<double>, 3>>& hess,
                               Vector<TinyVector<Tensor<std::complex<double>, 3>, 3>>& gradhess)
{
  eval_multi_UBspline_3d_z_vghgh(spline, r[0], r[1], r[2], psi.data(), (std::complex<double>*)grad.data(),
                                 (std::complex<double>*)hess.data(), (std::complex<double>*)gradhess.data());
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::resetParameters(const opt_variables_type& active)
{}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::setOrbitalSetSize(int norbs)
{
  OrbitalSetSize = norbs;
}

#if !defined(QMC_COMPLEX)
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateValue(const ParticleSet& P, int iat, RealValueVector_t& psi)
{
  ValueTimer.start();
  const PosType& r(P.activeR(iat));
  // Check atomic orbitals
  bool inAtom = false;
  for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
  {
    inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector);
    if (inAtom)
      break;
  }
  StorageValueVector_t& valVec = StorageValueVector;
  if (!inAtom)
  {
    PosType ru(PrimLattice.toUnit(r));
    for (int i = 0; i < OHMMS_DIM; i++)
      ru[i] -= std::floor(ru[i]);
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, valVec);
    EinsplineTimer.stop();
    // Add e^ikr phase to B-spline orbitals
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      PosType k = kPoints[j];
      double s, c;
      double phase = -dot(r, k);
      qmcplusplus::sincos(phase, &s, &c);
      std::complex<double> e_mikr(c, s);
      valVec[j] *= e_mikr;
    }
  }
  const int N  = StorageValueVector.size();
  int psiIndex = 0;
  for (int j = 0; j < N; j++)
  {
    std::complex<double> psi_val = StorageValueVector[j];
    psi[psiIndex]                = real(psi_val);
    psiIndex++;
    if (MakeTwoCopies[j])
    {
      psi[psiIndex] = imag(psi_val);
      psiIndex++;
    }
  }
  ValueTimer.stop();
}


// This is an explicit specialization of the above for real orbitals
// with a real return value, i.e. simulations at the gamma or L
// point.
template<>
void EinsplineSetExtended<double>::evaluateValue(const ParticleSet& P, int iat, RealValueVector_t& psi)
{
  ValueTimer.start();
  const PosType& r(P.activeR(iat));
  bool inAtom = false;
  for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
  {
    inAtom = AtomicOrbitals[jat].evaluate(r, psi);
    if (inAtom)
      break;
  }
  if (!inAtom)
  {
    PosType ru(PrimLattice.toUnit(r));
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      sign += HalfG[i] * (int)img;
    }
    // Check atomic orbitals
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, psi);
    EinsplineTimer.stop();
    if (sign & 1)
      for (int j = 0; j < psi.size(); j++)
        psi[j] *= -1.0;
  }
  ValueTimer.stop();
}


// Value, gradient, and laplacian
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateVGL(const ParticleSet& P,
                                                    int iat,
                                                    RealValueVector_t& psi,
                                                    RealGradVector_t& dpsi,
                                                    RealValueVector_t& d2psi)
{
  VGLTimer.start();
  const PosType& r(P.activeR(iat));
  std::complex<double> eye(0.0, 1.0);
  bool inAtom = false;
  for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
  {
    inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector, StorageGradVector, StorageLaplVector);
    if (inAtom)
      break;
  }
  StorageValueVector_t& valVec  = StorageValueVector;
  StorageGradVector_t& gradVec  = StorageGradVector;
  StorageValueVector_t& laplVec = StorageLaplVector;
  // Finally, copy into output vectors
  int psiIndex = 0;
  const int N  = StorageValueVector.size();
  for (int j = 0; j < N; j++)
  {
    std::complex<double> psi_val, psi_lapl;
    TinyVector<std::complex<double>, OHMMS_DIM> psi_grad;
    psi_val       = StorageValueVector[j];
    psi_grad      = StorageGradVector[j];
    psi_lapl      = StorageLaplVector[j];
    psi[psiIndex] = real(psi_val);
    for (int n = 0; n < OHMMS_DIM; n++)
      dpsi[psiIndex][n] = real(psi_grad[n]);
    d2psi[psiIndex] = real(psi_lapl);
    psiIndex++;
    if (MakeTwoCopies[j])
    {
      psi[psiIndex] = imag(psi_val);
      for (int n = 0; n < OHMMS_DIM; n++)
        dpsi[psiIndex][n] = imag(psi_grad[n]);
      d2psi[psiIndex] = imag(psi_lapl);
      psiIndex++;
    }
  }
  VGLTimer.stop();
}

template<>
void EinsplineSetExtended<double>::evaluateVGL(const ParticleSet& P,
                                               int iat,
                                               RealValueVector_t& psi,
                                               RealGradVector_t& dpsi,
                                               RealValueVector_t& d2psi)
{
  VGLTimer.start();
  const PosType& r(P.activeR(iat));
  bool inAtom = false;
  for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
  {
    inAtom = AtomicOrbitals[jat].evaluate(r, psi, dpsi, d2psi);
    if (inAtom)
      break;
  }
  if (!inAtom)
  {
    PosType ru(PrimLattice.toUnit(r));
    int sign = 0;
    for (int i = 0; i < OHMMS_DIM; i++)
    {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      sign += HalfG[i] * (int)img;
    }
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, psi, StorageGradVector, StorageHessVector);
    EinsplineTimer.stop();
    if (sign & 1)
      for (int j = 0; j < psi.size(); j++)
      {
        psi[j] *= -1.0;
        StorageGradVector[j] *= -1.0;
        StorageHessVector[j] *= -1.0;
      }
    for (int i = 0; i < psi.size(); i++)
    {
      dpsi[i]  = dot(PrimLattice.G, StorageGradVector[i]);
      d2psi[i] = trace(StorageHessVector[i], GGt);
    }
  }
  VGLTimer.stop();
}


template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             RealValueMatrix_t& psi,
                                                             RealGradMatrix_t& dpsi,
                                                             RealValueMatrix_t& d2psi)
{
  std::complex<double> eye(0.0, 1.0);
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    bool inAtom = false;
    for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
    {
      inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector, StorageGradVector, StorageLaplVector);
      if (inAtom)
        break;
    }
    StorageValueVector_t& valVec  = StorageValueVector;
    StorageGradVector_t& gradVec  = StorageGradVector;
    StorageValueVector_t& laplVec = StorageLaplVector;
    // Finally, copy into output vectors
    int psiIndex = 0;
    const int N  = StorageValueVector.size();
    for (int j = 0; j < N; j++)
    {
      std::complex<double> psi_val, psi_lapl;
      TinyVector<std::complex<double>, OHMMS_DIM> psi_grad;
      psi_val          = StorageValueVector[j];
      psi_grad         = StorageGradVector[j];
      psi_lapl         = StorageLaplVector[j];
      psi(i, psiIndex) = real(psi_val);
      for (int n = 0; n < OHMMS_DIM; n++)
        dpsi(i, psiIndex)[n] = real(psi_grad[n]);
      d2psi(i, psiIndex) = real(psi_lapl);
      psiIndex++;
      // if (psiIndex >= dpsi.cols()) {
      //   std::cerr << "Error:  out of bounds writing in EinsplineSet::evalate.\n"
      // 	 << "psiIndex = " << psiIndex << "  dpsi.cols() = " << dpsi.cols() << std::endl;
      // }
      if (MakeTwoCopies[j])
      {
        psi(i, psiIndex) = imag(psi_val);
        for (int n = 0; n < OHMMS_DIM; n++)
          dpsi(i, psiIndex)[n] = imag(psi_grad[n]);
        d2psi(i, psiIndex) = imag(psi_lapl);
        psiIndex++;
      }
    }
  }
  VGLMatTimer.stop();
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             RealValueMatrix_t& psi,
                                                             RealGradMatrix_t& dpsi,
                                                             RealHessMatrix_t& grad_grad_psi)
{
  std::complex<double> eye(0.0, 1.0);
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    bool inAtom = false;
    for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
    {
      inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector, StorageGradVector, StorageHessVector);
      if (inAtom)
        break;
    }
    StorageValueVector_t& valVec = StorageValueVector;
    StorageGradVector_t& gradVec = StorageGradVector;
    StorageHessVector_t& hessVec = StorageHessVector;
    Tensor<std::complex<double>, OHMMS_DIM> tmphs;
    // Finally, copy into output vectors
    int psiIndex = 0;
    const int N  = StorageValueVector.size();
    for (int j = 0; j < N; j++)
    {
      std::complex<double> psi_val;
      TinyVector<std::complex<double>, OHMMS_DIM> psi_grad;
      psi_val          = StorageValueVector[j];
      psi_grad         = StorageGradVector[j];
      tmphs            = StorageHessVector[j];
      psi(i, psiIndex) = real(psi_val);
      for (int n = 0; n < OHMMS_DIM; n++)
        dpsi(i, psiIndex)[n] = real(psi_grad[n]);
      //d2psi(i,psiIndex) = real(psi_lapl);
      // FIX FIX FIX
      for (int n = 0; n < OHMMS_DIM * OHMMS_DIM; n++)
        grad_grad_psi(i, psiIndex)[n] = real(tmphs(n));
      psiIndex++;
      // if (psiIndex >= dpsi.cols()) {
      //   std::cerr << "Error:  out of bounds writing in EinsplineSet::evalate.\n"
      //     << "psiIndex = " << psiIndex << "  dpsi.cols() = " << dpsi.cols() << std::endl;
      // }
      if (MakeTwoCopies[j])
      {
        psi(i, psiIndex) = imag(psi_val);
        for (int n = 0; n < OHMMS_DIM; n++)
          dpsi(i, psiIndex)[n] = imag(psi_grad[n]);
        //d2psi(i,psiIndex) = imag(psi_lapl);
        for (int n = 0; n < OHMMS_DIM * OHMMS_DIM; n++)
          grad_grad_psi(i, psiIndex)[n] = imag(tmphs(n));
        psiIndex++;
      }
    }
  }
  VGLMatTimer.stop();
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateGradSource(const ParticleSet& P,
                                                           int first,
                                                           int last,
                                                           const ParticleSet& source,
                                                           int iat,
                                                           RealGradMatrix_t& dpsi)
{
  if (ionDerivs)
  {
    // Loop over dimensions
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      // Loop over electrons
      for (int iel = first, i = 0; iel < last; iel++, i++)
      {
        const PosType& r(P.activeR(iel));
        PosType ru(PrimLattice.toUnit(r));
        assert(FirstOrderSplines[iat][dim]);
        EinsplineMultiEval(FirstOrderSplines[iat][dim], ru, StorageValueVector);
        int dpsiIndex = 0;
        for (int j = 0; j < NumValenceOrbs; j++)
        {
          PosType k = kPoints[j];
          double s, c;
          double phase = -dot(r, k);
          qmcplusplus::sincos(phase, &s, &c);
          std::complex<double> e_mikr(c, s);
          StorageValueVector[j] *= e_mikr;
          dpsi(i, dpsiIndex)[dim] = real(StorageValueVector[j]);
          dpsiIndex++;
          if (MakeTwoCopies[j])
          {
            dpsi(i, dpsiIndex)[dim] = imag(StorageValueVector[j]);
            dpsiIndex++;
          }
        }
      }
    }
    for (int i = 0; i < (last - first); i++)
      for (int j = 0; j < (last - first); j++)
        dpsi(i, j) = dot(PrimLattice.G, dpsi(i, j));
  }
}


// Evaluate the gradient w.r.t. to ion iat of the gradient and
// laplacian of the orbitals w.r.t. the electrons
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateGradSource(const ParticleSet& P,
                                                           int first,
                                                           int last,
                                                           const ParticleSet& source,
                                                           int iat_src,
                                                           RealGradMatrix_t& dphi,
                                                           RealHessMatrix_t& dgrad_phi,
                                                           RealGradMatrix_t& dlapl_phi)
{
  if (ionDerivs)
  {
    std::complex<double> eye(0.0, 1.0);
    // Loop over dimensions
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      // Loop over electrons
      for (int iel = first, i = 0; iel < last; iel++, i++)
      {
        const PosType& r(P.activeR(iel));
        PosType ru(PrimLattice.toUnit(r));
        assert(FirstOrderSplines[iat_src][dim]);
        EinsplineMultiEval(FirstOrderSplines[iat_src][dim], ru, StorageValueVector, StorageGradVector,
                           StorageHessVector);
        int dphiIndex = 0;
        for (int j = 0; j < NumValenceOrbs; j++)
        {
          StorageGradVector[j]                              = dot(PrimLattice.G, StorageGradVector[j]);
          StorageLaplVector[j]                              = trace(StorageHessVector[j], GGt);
          std::complex<double> u                            = StorageValueVector[j];
          TinyVector<std::complex<double>, OHMMS_DIM> gradu = StorageGradVector[j];
          std::complex<double> laplu                        = StorageLaplVector[j];
          PosType k                                         = kPoints[j];
          TinyVector<std::complex<double>, OHMMS_DIM> ck;
          for (int n = 0; n < OHMMS_DIM; n++)
            ck[n] = k[n];
          double s, c;
          double phase = -dot(r, k);
          qmcplusplus::sincos(phase, &s, &c);
          std::complex<double> e_mikr(c, s);
          StorageValueVector[j]   = e_mikr * u;
          StorageGradVector[j]    = e_mikr * (-eye * u * ck + gradu);
          StorageLaplVector[j]    = e_mikr * (-dot(k, k) * u - 2.0 * eye * dot(ck, gradu) + laplu);
          dphi(i, dphiIndex)[dim] = real(StorageValueVector[j]);
          for (int k = 0; k < OHMMS_DIM; k++)
            dgrad_phi(dphiIndex)[dim] = real(StorageGradVector[j][k]);
          dlapl_phi(dphiIndex)[dim] = real(StorageLaplVector[j]);
          dphiIndex++;
          if (MakeTwoCopies[j])
          {
            dphi(i, dphiIndex)[dim] = imag(StorageValueVector[j]);
            for (int k = 0; k < OHMMS_DIM; k++)
              dgrad_phi(i, dphiIndex)(dim, k) = imag(StorageGradVector[j][k]);
            dlapl_phi(i, dphiIndex)[dim] = imag(StorageLaplVector[j]);
            dphiIndex++;
          }
        }
      }
    }
    for (int i = 0; i < (last - first); i++)
      for (int j = 0; j < (last - first); j++)
      {
        dphi(i, j) = dot(PrimLattice.G, dphi(i, j));
        // Check this one!
        dgrad_phi(i, j) = dot(PrimLattice.G, dgrad_phi(i, j));
        dlapl_phi(i, j) = dot(PrimLattice.G, dlapl_phi(i, j));
      }
  }
}


template<>
void EinsplineSetExtended<double>::evaluateGradSource(const ParticleSet& P,
                                                      int first,
                                                      int last,
                                                      const ParticleSet& source,
                                                      int iat_src,
                                                      RealGradMatrix_t& dphi,
                                                      RealHessMatrix_t& dgrad_phi,
                                                      RealGradMatrix_t& dlapl_phi)
{
  if (ionDerivs)
  {
    // Loop over dimensions
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      assert(FirstOrderSplines[iat_src][dim]);
      // Loop over electrons
      for (int iel = first, i = 0; iel < last; iel++, i++)
      {
        const PosType& r(P.activeR(iel));
        PosType ru(PrimLattice.toUnit(r));
        int sign = 0;
        for (int n = 0; n < OHMMS_DIM; n++)
        {
          RealType img = std::floor(ru[n]);
          ru[n] -= img;
          sign += HalfG[n] * (int)img;
        }
        for (int n = 0; n < OHMMS_DIM; n++)
          ru[n] -= std::floor(ru[n]);
        EinsplineMultiEval(FirstOrderSplines[iat_src][dim], ru, StorageValueVector, StorageGradVector,
                           StorageHessVector);
        if (sign & 1)
          for (int j = 0; j < OrbitalSetSize; j++)
          {
            dphi(i, j)[dim] = -1.0 * StorageValueVector[j];
            PosType g       = -1.0 * dot(PrimLattice.G, StorageGradVector[j]);
            for (int k = 0; k < OHMMS_DIM; k++)
              dgrad_phi(i, j)(dim, k) = g[k];
            dlapl_phi(i, j)[dim] = -1.0 * trace(StorageHessVector[j], GGt);
          }
        else
          for (int j = 0; j < OrbitalSetSize; j++)
          {
            dphi(i, j)[dim] = StorageValueVector[j];
            PosType g       = dot(PrimLattice.G, StorageGradVector[j]);
            for (int k = 0; k < OHMMS_DIM; k++)
              dgrad_phi(i, j)(dim, k) = g[k];
            dlapl_phi(i, j)[dim] = trace(StorageHessVector[j], GGt);
          }
      }
    }
    for (int i = 0; i < (last - first); i++)
      for (int j = 0; j < (last - first); j++)
      {
        dphi(i, j) = dot(PrimLattice.G, dphi(i, j));
        // Check this one!
        dgrad_phi(i, j) = dot(PrimLattice.G, dgrad_phi(i, j));
        dlapl_phi(i, j) = dot(PrimLattice.G, dlapl_phi(i, j));
      }
  }
}
template<>
void EinsplineSetExtended<double>::evaluateGradSource(const ParticleSet& P,
                                                      int first,
                                                      int last,
                                                      const ParticleSet& source,
                                                      int iat,
                                                      RealGradMatrix_t& dpsi)
{
  if (ionDerivs)
  {
    // Loop over dimensions
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      assert(FirstOrderSplines[iat][dim]);
      // Loop over electrons
      for (int iel = first, i = 0; iel < last; iel++, i++)
      {
        const PosType& r(P.activeR(iel));
        PosType ru(PrimLattice.toUnit(r));
        int sign = 0;
        for (int n = 0; n < OHMMS_DIM; n++)
        {
          RealType img = std::floor(ru[n]);
          ru[n] -= img;
          sign += HalfG[n] * (int)img;
        }
        for (int n = 0; n < OHMMS_DIM; n++)
          ru[n] -= std::floor(ru[n]);
        EinsplineMultiEval(FirstOrderSplines[iat][dim], ru, StorageValueVector);
        if (sign & 1)
          for (int j = 0; j < OrbitalSetSize; j++)
            dpsi(i, j)[dim] = -1.0 * StorageValueVector[j];
        else
          for (int j = 0; j < OrbitalSetSize; j++)
            dpsi(i, j)[dim] = StorageValueVector[j];
      }
    }
    for (int i = 0; i < (last - first); i++)
      for (int j = 0; j < (last - first); j++)
      {
        dpsi(i, j) = dot(PrimLattice.G, dpsi(i, j));
      }
  }
}

#else

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateValue(const ParticleSet& P, int iat, ComplexValueVector_t& psi)
{
  ValueTimer.start();
  const PosType& r(P.activeR(iat));
  PosType ru(PrimLattice.toUnit(r));
  for (int i = 0; i < OHMMS_DIM; i++)
    ru[i] -= std::floor(ru[i]);
  EinsplineTimer.start();
  EinsplineMultiEval(MultiSpline, ru, StorageValueVector);
  EinsplineTimer.stop();
  //computePhaseFactors(r);
  for (int i = 0; i < psi.size(); i++)
  {
    PosType k = kPoints[i];
    double s, c;
    double phase = -dot(r, k);
    qmcplusplus::sincos(phase, &s, &c);
    std::complex<double> e_mikr(c, s);
    psi[i] = e_mikr * StorageValueVector[i];
  }
  ValueTimer.stop();
}

// Value, gradient, and laplacian
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateVGL(const ParticleSet& P,
                                                    int iat,
                                                    ComplexValueVector_t& psi,
                                                    ComplexGradVector_t& dpsi,
                                                    ComplexValueVector_t& d2psi)
{
  VGLTimer.start();
  const PosType& r(P.activeR(iat));
  PosType ru(PrimLattice.toUnit(r));
  for (int i = 0; i < OHMMS_DIM; i++)
    ru[i] -= std::floor(ru[i]);
  EinsplineTimer.start();
  EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
  EinsplineTimer.stop();
  //computePhaseFactors(r);
  std::complex<double> eye(0.0, 1.0);
  for (int j = 0; j < psi.size(); j++)
  {
    std::complex<double> u, laplu;
    TinyVector<std::complex<double>, OHMMS_DIM> gradu;
    u         = StorageValueVector[j];
    gradu     = dot(PrimLattice.G, StorageGradVector[j]);
    laplu     = trace(StorageHessVector[j], GGt);
    PosType k = kPoints[j];
    TinyVector<std::complex<double>, OHMMS_DIM> ck;
    for (int n = 0; n < OHMMS_DIM; n++)
      ck[n] = k[n];
    double s, c;
    double phase = -dot(r, k);
    qmcplusplus::sincos(phase, &s, &c);
    std::complex<double> e_mikr(c, s);
    psi[j]  = e_mikr * u;
    dpsi[j] = e_mikr * (-eye * u * ck + gradu);
    //convertVec(e_mikr*(-eye*u*ck + gradu), dpsi[j]);
    d2psi[j] = e_mikr * (-dot(k, k) * u - 2.0 * eye * dot(ck, gradu) + laplu);
  }
  VGLTimer.stop();
}

// Value, gradient, and laplacian
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluateVGH(const ParticleSet& P,
                                                    int iat,
                                                    ComplexValueVector_t& psi,
                                                    ComplexGradVector_t& dpsi,
                                                    ComplexHessVector_t& grad_grad_psi)
{
  VGLTimer.start();
  const PosType& r(P.activeR(iat));
  PosType ru(PrimLattice.toUnit(r));
  for (int i = 0; i < OHMMS_DIM; i++)
    ru[i] -= std::floor(ru[i]);
  EinsplineTimer.start();
  EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
  EinsplineTimer.stop();
  //computePhaseFactors(r);
  std::complex<double> eye(0.0, 1.0);
  for (int j = 0; j < psi.size(); j++)
  {
    std::complex<double> u;
    TinyVector<std::complex<double>, OHMMS_DIM> gradu;
    Tensor<std::complex<double>, OHMMS_DIM> hs, tmphs;
    u     = StorageValueVector[j];
    gradu = dot(PrimLattice.G, StorageGradVector[j]);
    ////laplu = trace(StorageHessVector[j], GGt);
    tmphs = dot(PrimLattice.G, StorageHessVector[j]);
    //hs = dot(tmphs,PrimLattice.G);
    hs        = dot(tmphs, PrimLattice.Gt);
    PosType k = kPoints[j];
    TinyVector<std::complex<double>, OHMMS_DIM> ck;
    for (int n = 0; n < OHMMS_DIM; n++)
      ck[n] = k[n];
    double s, c;
    double phase = -dot(r, k);
    qmcplusplus::sincos(phase, &s, &c);
    std::complex<double> e_mikr(c, s);
    psi[j]  = e_mikr * u;
    dpsi[j] = e_mikr * (-eye * u * ck + gradu);
    //convertVec(e_mikr*(-eye*u*ck + gradu), dpsi[j]);
    //d2psi[j] = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
    grad_grad_psi[j] =
        e_mikr * (hs - u * outerProduct(ck, ck) - eye * outerProduct(ck, gradu) - eye * outerProduct(gradu, ck));
  }
  VGLTimer.stop();
}
#endif

#if !defined(QMC_COMPLEX)
template<>
void EinsplineSetExtended<double>::evaluate_notranspose(const ParticleSet& P,
                                                        int first,
                                                        int last,
                                                        RealValueMatrix_t& psi,
                                                        RealGradMatrix_t& dpsi,
                                                        RealValueMatrix_t& d2psi)
{
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    bool inAtom = false;
    for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
    {
      inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector, StorageGradVector, StorageLaplVector);
      if (inAtom)
      {
        for (int j = 0; j < OrbitalSetSize; j++)
        {
          psi(i, j)   = StorageValueVector[j];
          dpsi(i, j)  = StorageGradVector[j];
          d2psi(i, j) = StorageLaplVector[j];
        }
        break;
      }
    }
    if (!inAtom)
    {
      PosType ru(PrimLattice.toUnit(r));
      int sign = 0;
      for (int n = 0; n < OHMMS_DIM; n++)
      {
        RealType img = std::floor(ru[n]);
        ru[n] -= img;
        sign += HalfG[n] * (int)img;
      }
      for (int n = 0; n < OHMMS_DIM; n++)
        ru[n] -= std::floor(ru[n]);
      EinsplineTimer.start();
      EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
      EinsplineTimer.stop();
      if (sign & 1)
        for (int j = 0; j < OrbitalSetSize; j++)
        {
          StorageValueVector[j] *= -1.0;
          StorageGradVector[j] *= -1.0;
          StorageHessVector[j] *= -1.0;
        }
      for (int j = 0; j < OrbitalSetSize; j++)
      {
        psi(i, j)   = StorageValueVector[j];
        dpsi(i, j)  = dot(PrimLattice.G, StorageGradVector[j]);
        d2psi(i, j) = trace(StorageHessVector[j], GGt);
      }
    }
  }
  VGLMatTimer.stop();
}

template<>
void EinsplineSetExtended<double>::evaluate_notranspose(const ParticleSet& P,
                                                        int first,
                                                        int last,
                                                        RealValueMatrix_t& psi,
                                                        RealGradMatrix_t& dpsi,
                                                        RealHessMatrix_t& grad_grad_psi)
{
  //APP_ABORT("evaluate_notranspose:  Check Hessian, then remove this error message.\n")
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    bool inAtom = false;
    for (int jat = 0; jat < AtomicOrbitals.size(); jat++)
    {
      inAtom = AtomicOrbitals[jat].evaluate(r, StorageValueVector, StorageGradVector, StorageHessVector);
      if (inAtom)
      {
        for (int j = 0; j < OrbitalSetSize; j++)
        {
          psi(i, j)           = StorageValueVector[j];
          dpsi(i, j)          = StorageGradVector[j];
          grad_grad_psi(i, j) = StorageHessVector[j];
        }
        break;
      }
    }
    if (!inAtom)
    {
      PosType ru(PrimLattice.toUnit(r));
      int sign = 0;
      for (int n = 0; n < OHMMS_DIM; n++)
      {
        RealType img = std::floor(ru[n]);
        ru[n] -= img;
        sign += HalfG[n] * (int)img;
      }
      for (int n = 0; n < OHMMS_DIM; n++)
        ru[n] -= std::floor(ru[n]);
      EinsplineTimer.start();
      EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
      EinsplineTimer.stop();
      if (sign & 1)
        for (int j = 0; j < OrbitalSetSize; j++)
        {
          StorageValueVector[j] *= -1.0;
          StorageGradVector[j] *= -1.0;
          StorageHessVector[j] *= -1.0;
        }
      for (int j = 0; j < OrbitalSetSize; j++)
      {
        psi(i, j)           = StorageValueVector[j];
        dpsi(i, j)          = dot(PrimLattice.G, StorageGradVector[j]);
        grad_grad_psi(i, j) = dot(PrimLattice.G, dot(StorageHessVector[j], PrimLattice.Gt));
      }
    }
  }
  VGLMatTimer.stop();
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             RealValueMatrix_t& psi,
                                                             RealGradMatrix_t& dpsi,
                                                             RealHessMatrix_t& grad_grad_psi,
                                                             RealGGGMatrix_t& grad_grad_grad_logdet)
{
  //      APP_ABORT(" EinsplineSetExtended<StorageType>::evaluate_notranspose not implemented for grad_grad_grad_logdet yet. \n");
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    PosType ru(PrimLattice.toUnit(r));
    for (int n = 0; n < OHMMS_DIM; n++)
      ru[n] -= std::floor(ru[n]);
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector,
                       StorageGradHessVector);
    EinsplineTimer.stop();
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      TinyVector<std::complex<double>, OHMMS_DIM> tmpg;
      Tensor<std::complex<double>, OHMMS_DIM> tmphs;
      TinyVector<Tensor<std::complex<double>, OHMMS_DIM>, OHMMS_DIM> tmpghs;
      tmpg                 = dot(PrimLattice.G, StorageGradVector[j]);
      StorageGradVector[j] = tmpg;
      tmphs                = dot(PrimLattice.G, StorageHessVector[j]);
      StorageHessVector[j] = dot(tmphs, PrimLattice.Gt);
      for (int n = 0; n < OHMMS_DIM; n++)
      {
        tmphs     = dot(PrimLattice.G, StorageGradHessVector[j][n]);
        tmpghs[n] = dot(tmphs, PrimLattice.Gt);
      }
      StorageGradHessVector[j] = dot(PrimLattice.G, tmpghs);
    }
    std::complex<double> eye(0.0, 1.0);
    //    StorageValueVector_t &valVec =
    //      StorageValueVector;
    //    StorageGradVector_t &gradVec =
    //      StorageGradVector;
    //    StorageHessVector_t &hessVec =
    //      StorageHessVector;
    //    Tensor<std::complex<double>,OHMMS_DIM> tmphs;
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      //            std::complex<double> u = valVec[j];
      //            TinyVector<std::complex<double>,OHMMS_DIM> gradu = gradVec[j];
      //            tmphs = hessVec[j];
      //            PosType k = kPoints[j];
      //            TinyVector<std::complex<double>,OHMMS_DIM> ck;
      //            for (int n=0; n<OHMMS_DIM; n++)       ck[n] = k[n];
      //            double s,c;
      //            double phase = -dot(r, k);
      //            qmcplusplus::sincos (phase, &s, &c);
      //            std::complex<double> e_mikr (c,s);
      //            valVec[j]   = e_mikr*u;
      //            gradVec[j]  = e_mikr*(-eye*u*ck + gradu);
      //            hessVec[j]  = e_mikr*(tmphs -u*outerProduct(ck,ck) - eye*outerProduct(ck,gradu) - eye*outerProduct(gradu,ck));
      std::complex<double> u                            = (StorageValueVector[j]);
      TinyVector<std::complex<double>, OHMMS_DIM> gradu = (StorageGradVector[j]);
      Tensor<std::complex<double>, OHMMS_DIM> tmphs     = (StorageHessVector[j]);
      //        TinyVector<Tensor<std::complex<double>,OHMMS_DIM>,OHMMS_DIM> tmpghs=(StorageGradHessVector[j]);
      PosType k = kPoints[j];
      TinyVector<std::complex<double>, OHMMS_DIM> ck;
      for (int n = 0; n < OHMMS_DIM; n++)
        ck[n] = k[n];
      double s, c;
      double phase = -dot(r, k);
      qmcplusplus::sincos(phase, &s, &c);
      std::complex<double> e_mikr(c, s);
      StorageValueVector[j] = e_mikr * u;
      StorageGradVector[j]  = e_mikr * (-eye * u * ck + gradu);
      StorageHessVector[j] =
          e_mikr * (tmphs - u * outerProduct(ck, ck) - eye * outerProduct(ck, gradu) - eye * outerProduct(gradu, ck));
      //Is this right?
      StorageGradHessVector[j] *= e_mikr;
      for (unsigned a0(0); a0 < OHMMS_DIM; a0++)
        for (unsigned a1(0); a1 < OHMMS_DIM; a1++)
          for (unsigned a2(0); a2 < OHMMS_DIM; a2++)
            StorageGradHessVector[j][a0](a1, a2) += e_mikr *
                (-1.0 * eye * (ck[a0] * tmphs(a1, a2) + ck[a1] * tmphs(a0, a2) + ck[a2] * tmphs(a0, a1)) -
                 (ck[a0] * ck[a1] * gradu[a2] + ck[a0] * ck[a2] * gradu[a1] + ck[a1] * ck[a2] * gradu[a0]) +
                 eye * ck[a0] * ck[a1] * ck[a2] * u);
    }
    int psiIndex(0);
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      if (MakeTwoCopies[j])
      {
        psi(i, psiIndex) = imag(StorageValueVector[j]);
        for (int n = 0; n < OHMMS_DIM; n++)
          dpsi(i, psiIndex)[n] = imag(StorageGradVector[j][n]);
        for (int n = 0; n < OHMMS_DIM * OHMMS_DIM; n++)
          grad_grad_psi(i, psiIndex)[n] = imag(StorageHessVector[j](n));
        for (int n = 0; n < OHMMS_DIM; n++)
          for (int m = 0; m < OHMMS_DIM * OHMMS_DIM; m++)
            grad_grad_grad_logdet(i, psiIndex)[n][m] = imag(StorageGradHessVector[j][n](m));
        psiIndex++;
        psi(i, psiIndex) = real(StorageValueVector[j]);
        for (int n = 0; n < OHMMS_DIM; n++)
          dpsi(i, psiIndex)[n] = real(StorageGradVector[j][n]);
        for (int n = 0; n < OHMMS_DIM * OHMMS_DIM; n++)
          grad_grad_psi(i, psiIndex)[n] = real(StorageHessVector[j](n));
        for (int n = 0; n < OHMMS_DIM; n++)
          for (int m = 0; m < OHMMS_DIM * OHMMS_DIM; m++)
            grad_grad_grad_logdet(i, psiIndex)[n][m] = real(StorageGradHessVector[j][n](m));
        psiIndex++;
      }
      else
      {
        psi(i, psiIndex) = real(StorageValueVector[j]);
        for (int n = 0; n < OHMMS_DIM; n++)
          dpsi(i, psiIndex)[n] = real(StorageGradVector[j][n]);
        for (int n = 0; n < OHMMS_DIM * OHMMS_DIM; n++)
          grad_grad_psi(i, psiIndex)[n] = real(StorageHessVector[j](n));
        for (int n = 0; n < OHMMS_DIM; n++)
          for (int m = 0; m < OHMMS_DIM * OHMMS_DIM; m++)
            grad_grad_grad_logdet(i, psiIndex)[n][m] = real(StorageGradHessVector[j][n](m));
        psiIndex++;
      }
    }
  }
  VGLMatTimer.stop();
}

template<>
void EinsplineSetExtended<double>::evaluate_notranspose(const ParticleSet& P,
                                                        int first,
                                                        int last,
                                                        RealValueMatrix_t& psi,
                                                        RealGradMatrix_t& dpsi,
                                                        RealHessMatrix_t& grad_grad_psi,
                                                        RealGGGMatrix_t& grad_grad_grad_logdet)
{
  //      APP_ABORT(" EinsplineSetExtended<StorageType>::evaluate_notranspose not implemented for grad_grad_grad_logdet yet. \n");
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    bool inAtom  = false;
    int psiIndex = 0;
    const int N  = StorageValueVector.size();
    for (int j = 0; j < N; j++)
    {
      psi(i, psiIndex)                   = StorageValueVector[j];
      dpsi(i, psiIndex)                  = dot(StorageGradVector[j], PrimLattice.G);
      grad_grad_psi(i, psiIndex)         = StorageHessVector[j];
      grad_grad_grad_logdet(i, psiIndex) = dot(StorageGradHessVector[j], PrimLattice.G);
      psiIndex++;
    }
  }
  VGLMatTimer.stop();
}


#else
template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             ComplexValueMatrix_t& psi,
                                                             ComplexGradMatrix_t& dpsi,
                                                             ComplexValueMatrix_t& d2psi)
{
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    PosType ru(PrimLattice.toUnit(r));
    for (int n = 0; n < OHMMS_DIM; n++)
      ru[n] -= std::floor(ru[n]);
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
    EinsplineTimer.stop();
    //computePhaseFactors(r);
    std::complex<double> eye(0.0, 1.0);
    for (int j = 0; j < OrbitalSetSize; j++)
    {
      std::complex<double> u, laplu;
      TinyVector<std::complex<double>, OHMMS_DIM> gradu;
      u         = StorageValueVector[j];
      gradu     = dot(PrimLattice.G, StorageGradVector[j]);
      laplu     = trace(StorageHessVector[j], GGt);
      PosType k = kPoints[j];
      TinyVector<std::complex<double>, OHMMS_DIM> ck;
      for (int n = 0; n < OHMMS_DIM; n++)
        ck[n] = k[n];
      double s, c;
      double phase = -dot(r, k);
      qmcplusplus::sincos(phase, &s, &c);
      std::complex<double> e_mikr(c, s);
      psi(i, j) = e_mikr * u;
      //psi(j,i) = e_mikr * u;
      dpsi(i, j) = e_mikr * (-eye * u * ck + gradu);
      //convertVec(e_mikr*(-eye*u*ck + gradu), dpsi(i,j));
      d2psi(i, j) = e_mikr * (-dot(k, k) * u - 2.0 * eye * dot(ck, gradu) + laplu);
    }
  }
  VGLMatTimer.stop();
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             ComplexValueMatrix_t& psi,
                                                             ComplexGradMatrix_t& dpsi,
                                                             ComplexHessMatrix_t& grad_grad_psi)
{
  VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    PosType ru(PrimLattice.toUnit(r));
    for (int n = 0; n < OHMMS_DIM; n++)
      ru[n] -= std::floor(ru[n]);
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector);
    EinsplineTimer.stop();
    //computePhaseFactors(r);
    std::complex<double> eye(0.0, 1.0);
    for (int j = 0; j < OrbitalSetSize; j++)
    {
      std::complex<double> u;
      TinyVector<std::complex<double>, OHMMS_DIM> gradu;
      Tensor<std::complex<double>, OHMMS_DIM> hs, tmphs;
      u     = StorageValueVector[j];
      gradu = dot(PrimLattice.G, StorageGradVector[j]);
      // tmphs = dot(transpose(PrimLattice.G),StorageHessVector[j]);
      tmphs = dot(PrimLattice.G, StorageHessVector[j]);
      hs    = dot(tmphs, PrimLattice.Gt);
      //laplu = trace(StorageHessVector[j], GGt);
      PosType k = kPoints[j];
      TinyVector<std::complex<double>, OHMMS_DIM> ck;
      for (int n = 0; n < OHMMS_DIM; n++)
        ck[n] = k[n];
      double s, c;
      double phase = -dot(r, k);
      qmcplusplus::sincos(phase, &s, &c);
      std::complex<double> e_mikr(c, s);
      psi(i, j) = e_mikr * u;
      //psi(j,i) = e_mikr * u;
      dpsi(i, j) = e_mikr * (-eye * u * ck + gradu);
      //convertVec(e_mikr*(-eye*u*ck + gradu), dpsi(i,j));
      grad_grad_psi(i, j) =
          e_mikr * (hs - u * outerProduct(ck, ck) - eye * outerProduct(ck, gradu) - eye * outerProduct(gradu, ck));
    }
  }
  VGLMatTimer.stop();
}

template<>
void EinsplineSetExtended<double>::evaluate_notranspose(const ParticleSet& P,
                                                        int first,
                                                        int last,
                                                        ComplexValueMatrix_t& psi,
                                                        ComplexGradMatrix_t& dpsi,
                                                        ComplexHessMatrix_t& grad_grad_psi,
                                                        ComplexGGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT(
      " EinsplineSetExtended<StorageType>::evaluate_notranspose not implemented for grad_grad_grad_logdet yet. \n");
}

template<typename StorageType>
void EinsplineSetExtended<StorageType>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             ComplexValueMatrix_t& psi,
                                                             ComplexGradMatrix_t& dpsi,
                                                             ComplexHessMatrix_t& grad_grad_psi,
                                                             ComplexGGGMatrix_t& grad_grad_grad_logdet)
{
  //      APP_ABORT(" EinsplineSetExtended<StorageType>::evaluate_notranspose not implemented for grad_grad_grad_logdet yet. \n");
  //    VGLMatTimer.start();
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    const PosType& r(P.activeR(iat));
    PosType ru(PrimLattice.toUnit(r));
    for (int n = 0; n < OHMMS_DIM; n++)
      ru[n] -= std::floor(ru[n]);
    EinsplineTimer.start();
    EinsplineMultiEval(MultiSpline, ru, StorageValueVector, StorageGradVector, StorageHessVector,
                       StorageGradHessVector);
    EinsplineTimer.stop();
    Tensor<double, OHMMS_DIM> PG;
    PG = PrimLattice.G;
    Tensor<double, OHMMS_DIM> TPG;
    TPG = transpose(PrimLattice.G);
    Tensor<std::complex<double>, OHMMS_DIM> hs, tmphs;
    TinyVector<Tensor<std::complex<double>, OHMMS_DIM>, OHMMS_DIM> tmpghs, hvdot;
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      StorageGradVector[j] = dot(PG, StorageGradVector[j]);
      tmphs                = dot(PG, StorageHessVector[j]);
      StorageHessVector[j] = dot(tmphs, TPG);
      for (int n = 0; n < OHMMS_DIM; n++)
      {
        tmpghs[n]                   = dot(PG, StorageGradHessVector[j][n]);
        StorageGradHessVector[j][n] = dot(tmpghs[n], TPG);
      }
      StorageGradHessVector[j] = dot(PG, StorageGradHessVector[j]);
      //              grad_grad_grad_logdet(i,j)=StorageGradHessVector[j];
      //              grad_grad_psi(i,j)=StorageHessVector[j];
      //              dpsi(i,j)=StorageGradVector[j];
      //              psi(i,j)=StorageValueVector[j];
    }
    const std::complex<double> eye(0.0, 1.0);
    const std::complex<double> meye(0.0, -1.0);
    for (int j = 0; j < NumValenceOrbs; j++)
    {
      std::complex<double> u(StorageValueVector[j]);
      TinyVector<std::complex<double>, OHMMS_DIM> gradu(StorageGradVector[j]);
      tmphs     = StorageHessVector[j];
      tmpghs    = StorageGradHessVector[j];
      PosType k = kPoints[j];
      TinyVector<double, OHMMS_DIM> ck;
      for (int n = 0; n < OHMMS_DIM; n++)
        ck[n] = k[n];
      double s, c;
      double phase = -dot(r, k);
      qmcplusplus::sincos(phase, &s, &c);
      std::complex<double> e_mikr(c, s);
      psi(i, j)  = e_mikr * u;
      dpsi(i, j) = e_mikr * (-eye * u * ck + gradu);
      grad_grad_psi(i, j) =
          e_mikr * (tmphs - u * outerProduct(ck, ck) - eye * outerProduct(ck, gradu) - eye * outerProduct(gradu, ck));
      //Is this right?
      StorageGradHessVector[j] *= e_mikr;
      for (unsigned a0(0); a0 < OHMMS_DIM; a0++)
        for (unsigned a1(0); a1 < OHMMS_DIM; a1++)
          for (unsigned a2(0); a2 < OHMMS_DIM; a2++)
            StorageGradHessVector[j][a0](a1, a2) += e_mikr *
                (meye * (ck[a0] * tmphs(a1, a2) + ck[a1] * tmphs(a0, a2) + ck[a2] * tmphs(a0, a1)) -
                 (ck[a0] * ck[a1] * gradu[a2] + ck[a0] * ck[a2] * gradu[a1] + ck[a1] * ck[a2] * gradu[a0]) +
                 eye * ck[a0] * ck[a1] * ck[a2] * u);
      grad_grad_grad_logdet(i, j) = StorageGradHessVector[j];
    }
  }
}

#endif

template<typename StorageType>
std::string EinsplineSetExtended<StorageType>::Type()
{
  return "EinsplineSetExtended";
}


template<typename StorageType>
void EinsplineSetExtended<StorageType>::registerTimers()
{
  ValueTimer.reset();
  VGLTimer.reset();
  VGLMatTimer.reset();
  EinsplineTimer.reset();
}


template<typename StorageType>
std::unique_ptr<SPOSet> EinsplineSetExtended<StorageType>::makeClone() const
{
  auto clone = std::make_unique<EinsplineSetExtended<StorageType>>(*this);
  clone->registerTimers();
  for (int iat = 0; iat < clone->AtomicOrbitals.size(); iat++)
    clone->AtomicOrbitals[iat].registerTimers();
  return clone;
}

template class EinsplineSetExtended<std::complex<double>>;
template class EinsplineSetExtended<double>;


#ifdef QMC_CUDA
///////////////////////////////
// Real StorageType versions //
///////////////////////////////
template<>
std::string EinsplineSetHybrid<double>::Type()
{
  return "EinsplineSetHybrid<double>";
}


template<typename StorageType>
std::unique_ptr<SPOSet> EinsplineSetHybrid<StorageType>::makeClone() const
{
  auto clone = std::make_unique<EinsplineSetHybrid<StorageType>>(*this);
  clone->registerTimers();
  return clone;
}


//////////////////////////////////
// Complex StorageType versions //
//////////////////////////////////


template<>
std::string EinsplineSetHybrid<std::complex<double>>::Type()
{
  return "EinsplineSetHybrid<std::complex<double> >";
}

template<>
EinsplineSetHybrid<double>::EinsplineSetHybrid() : CurrentWalkers(0)
{
  className = "EinsplineSeHybrid";
  for (int i = 0; i < OHMMS_DIM; i++)
    HalfG[i] = 0;
}

template<>
EinsplineSetHybrid<std::complex<double>>::EinsplineSetHybrid() : CurrentWalkers(0)
{
  className = "EinsplineSeHybrid";
  for (int i = 0; i < OHMMS_DIM; i++)
    HalfG[i] = 0;
}

template class EinsplineSetHybrid<std::complex<double>>;
template class EinsplineSetHybrid<double>;
#endif
} // namespace qmcplusplus
