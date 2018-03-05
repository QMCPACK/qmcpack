//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"
#include "LongRange/StructFact.h"
#include "Numerics/e2iphi.h"
#include <sstream>
#include <algorithm>


namespace qmcplusplus
{

void
kSpaceJastrow::StructureFactor(PosType G, std::vector<ComplexType> &rho_G)
{
  for (int i=0; i<NumIonSpecies; i++)
    rho_G[i] = ComplexType();
  for (int iat=0; iat<Ions.getTotalNum(); iat++)
  {
    PosType r(Ions.R[iat]);
    RealType phase = dot(r,G);
    int id = Ions.GroupID[iat];
    rho_G[id] += ComplexType(std::cos(phase), std::sin(phase));
  }
}

inline bool
Include (int i, int j, int k)
{
  if (i > 0)
    return true;
  else
    if (i==0)
    {
      if (j>0)
        return true;
      else
        if ((j==0) && (k>0))
          return true;
    }
  return false;
}

inline bool
Include (int i, int j)
{
  if (i > 0)
    return true;
  else
    if (i==0)
    {
      if (j>0)
        return true;
    }
  return false;
}

void
kSpaceJastrow::setupGvecs(RealType kc, std::vector<PosType> &gvecs,
                          bool useStructFact)
{
  gvecs.clear();
  int maxIndex[OHMMS_DIM];
  for (int i=0; i<OHMMS_DIM; i++)
    maxIndex[i] = 2 + (int)std::floor
                  (std::sqrt(dot(Ions.Lattice.a(i), Ions.Lattice.a(i)))*kc / (2.0*M_PI));
  std::vector<ComplexType> rho_G(NumIonSpecies);
#if OHMMS_DIM==3
  for (int i=0; i<=maxIndex[0]; i++)
    for (int j=-maxIndex[1]; j<=maxIndex[1]; j++)
      for (int k=-maxIndex[2]; k<=maxIndex[2]; k++)
      {
        // Omit half the G-vectors because of time-reversal symmetry
        if (Include(i,j,k))
        {
          PosType G = 2.0*M_PI*((RealType)i*Ions.Lattice.Gv[0] +
                                (RealType)j*Ions.Lattice.Gv[1] +
                                (RealType)k*Ions.Lattice.Gv[2]);
          if (dot(G,G) <= (kc*kc))
          {
            bool notZero(false);
            StructureFactor(G, rho_G);
            for (int sp=0; sp<NumIonSpecies; sp++)
              notZero = notZero || (norm(rho_G[sp]) > 1.0e-12);
            if (notZero || !useStructFact)
              gvecs.push_back(G);
          }
        }
      }
#elif OHMMS_DIM==2
  for (int i=0; i<=maxIndex[0]; i++)
    for (int j=-maxIndex[1]; j<=maxIndex[1]; j++)
    {
      // Omit half the G-vectors because of time-reversal symmetry
      if (Include(i,j))
      {
        PosType G = 2.0*M_PI*((RealType)i*Ions.Lattice.Gv[0] +
                              (RealType)j*Ions.Lattice.Gv[1]);
        if (dot(G,G) <= (kc*kc))
        {
          bool notZero(false);
          StructureFactor(G, rho_G);
          for (int sp=0; sp<NumIonSpecies; sp++)
            notZero = notZero || (norm(rho_G[sp]) > 1.0e-12);
          if (notZero || !useStructFact)
            gvecs.push_back(G);
        }
      }
    }
#endif
}

struct magLess
{
  inline bool operator()(kSpaceJastrow::PosType a, kSpaceJastrow::PosType b)
  {
    return dot(a,a) < dot(b,b);
  }
};

bool
kSpaceJastrow::operator()(PosType G1, PosType G2)
{
  if (std::abs(dot(G1,G1) - dot(G2,G2)) > 1.0e-8)
    return dot(G1,G1) < dot(G2,G2);
  // if (Equivalent(G1,G2)) return false;
  std::vector<ComplexType> rho_G1(NumIonSpecies), rho_G2(NumIonSpecies);
  StructureFactor(G1, rho_G1);
  StructureFactor(G2, rho_G2);
  for (int i=0; i<NumIonSpecies; i++ )
    for (int j=i+1; j<NumIonSpecies; j++)
    {
      ComplexType zG1 = rho_G1[i]*qmcplusplus::conj(rho_G1[j]);
      ComplexType zG2 = rho_G2[i]*qmcplusplus::conj(rho_G2[j]);
      double SG1  = std::real(zG1);
      double SG2  = std::real(zG2);
      if (std::abs(SG1 - SG2) > 1.0e-8)
        return SG1 < SG2;
    }
  return false;
}

bool
kSpaceJastrow::Equivalent(PosType G1, PosType G2)
{
  return (!(*this)(G1,G2) && !(*this)(G2,G1));
  if (std::abs(dot(G1,G1) - dot(G2,G2)) > 1.0e-8)
    return false;
  std::vector<ComplexType> rho_G1(NumIonSpecies), rho_G2(NumIonSpecies);
  StructureFactor(G1, rho_G1);
  StructureFactor(G2, rho_G2);
  for (int j=0; j<NumIonSpecies; j++)
    for (int i=j; i<NumIonSpecies; i++ )
    {
      ComplexType zG1 = rho_G1[i]*qmcplusplus::conj(rho_G1[j]);
      ComplexType zG2 = rho_G2[i]*qmcplusplus::conj(rho_G2[j]);
      double SG1  = std::real(zG1);
      double SG2  = std::real(zG2);
      if (std::abs(SG1 - SG2) > 1.0e-8)
        return false;
    }
  return true;
}

template<typename T> void
kSpaceJastrow::sortGvecs(std::vector<PosType> &gvecs,
                         std::vector<kSpaceCoef<T> > &coefs,
                         SymmetryType symm)
{
  if (!gvecs.size())
    return;
  if (symm == CRYSTAL)
  {
    // First pass:  sort all the G-vectors into equivalent groups
    std::sort(gvecs.begin(), gvecs.end(), *this);
    // Now, look through the sorted G-vectors and group them
    kSpaceCoef<T> coef;
    coef.cG = T();
    coef.firstIndex = coef.lastIndex = 0;
    for (int i=1; i<gvecs.size(); i++)
      if (Equivalent(gvecs[i], gvecs[i-1]))
        coef.lastIndex=i;
      else
      {
        coefs.push_back(coef);
        coef.firstIndex = coef.lastIndex = i;
      }
    coefs.push_back(coef);
  }
  else
    if (symm == ISOTROPIC)
    {
      magLess comparator;
      std::sort(gvecs.begin(), gvecs.end(), comparator);
      double curMag2 = dot(gvecs[0], gvecs[0]);
      kSpaceCoef<T> coef;
      coef.cG = T();
      coef.firstIndex = coef.lastIndex = 0;
      for (int i=1; i<gvecs.size(); i++)
      {
        double mag2 = dot(gvecs[i],gvecs[i]);
        if (std::abs(mag2-curMag2) < 1.0e-10)
          coef.lastIndex = i;
        else
        {
          coefs.push_back(coef);
          coef.firstIndex = coef.lastIndex = i;
          curMag2 = mag2;
        }
      }
      coefs.push_back(coef);
    }
    else
      if (symm == NOSYMM)
      {
        coefs.resize(gvecs.size());
        for (int i=0; i<gvecs.size(); i++)
        {
          coefs[i].cG = T();
          coefs[i].firstIndex = coefs[i].lastIndex = i;
        }
      }
  app_log() << "Using a total of " << gvecs.size() << " G-vectors in " << coefs.size()
            << " symmetry groups.\n";
  app_log() << "kSpace coefficent groups:\n";
  for (int i=0; i<coefs.size(); i++)
  {
    app_log() << "  Group " << i << ":\n";
    for (int j=coefs[i].firstIndex; j<=coefs[i].lastIndex; j++)
      app_log() << "    " << gvecs[j] << "    " << std::sqrt(dot(gvecs[j], gvecs[j])) << std::endl;
  }
}

kSpaceJastrow::kSpaceJastrow(ParticleSet& ions, ParticleSet& elecs,
                             SymmetryType oneBodySymm, RealType oneBodyCutoff,
                             std::string onebodyid, bool oneBodySpin,
                             SymmetryType twoBodySymm, RealType twoBodyCutoff,
                             std::string twobodyid, bool twoBodySpin)
  : Ions(ions), Elecs(elecs),OneBodyID(onebodyid),TwoBodyID(twobodyid)
{
  Optimizable=true;
  Prefactor = 1.0/elecs.Lattice.Volume;
  NumIonSpecies = 0;
  NumElecs = elecs.getTotalNum();
  for (int iat=0; iat<ions.getTotalNum(); iat++)
    NumIonSpecies = std::max(NumIonSpecies, ions.GroupID[iat]+1);
  if (oneBodyCutoff > 0.0)
  {
    setupGvecs(oneBodyCutoff, OneBodyGvecs, true);
    sortGvecs (OneBodyGvecs, OneBodySymmCoefs, oneBodySymm);
    for (int i=0; i<OneBodySymmCoefs.size(); i++)
    {
      std::stringstream name_real, name_imag;
      name_real << OneBodyID << "_" << 2*i;
      name_imag << OneBodyID << "_" << 2*i+1;
      myVars.insert(name_real.str(),OneBodySymmCoefs[i].cG.real(),true,optimize::LOGLINEAR_P);
      myVars.insert(name_imag.str(),OneBodySymmCoefs[i].cG.imag(),true,optimize::LOGLINEAR_P);
      //VarMap[name_real.str()] = &(OneBodySymmCoefs[i].cG.real());
      //VarMap[name_imag.str()] = &(OneBodySymmCoefs[i].cG.imag());
    }
  }
  if (twoBodyCutoff > 0.0)
  {
    setupGvecs(twoBodyCutoff, TwoBodyGvecs, false);
    sortGvecs (TwoBodyGvecs, TwoBodySymmCoefs, twoBodySymm);
    for (int i=0; i<TwoBodySymmCoefs.size(); i++)
    {
      std::stringstream name;
      name << TwoBodyID << "_" << i;
      myVars.insert(name.str(),TwoBodySymmCoefs[i].cG,true,optimize::LOGLINEAR_P);
      //VarMap[name.str()] = &(TwoBodySymmCoefs[i].cG);
    }
  }
  if (oneBodySpin)
    app_log() << "One-body k-space Jastrow is spin-dependent.\n";
  if (twoBodySpin)
    app_log() << "Two-body k-space Jastrow is spin-dependent.\n";
  // Now resize all buffers
  int nOne = OneBodyGvecs.size();
  OneBodyCoefs.resize(nOne);
  OneBody_rhoG.resize(nOne);
  Ion_rhoG.resize(nOne);
  OneBodyPhase.resize(nOne);
  OneBody_e2iGr.resize(nOne);
  int nTwo = TwoBodyGvecs.size();
  int nElecs = Elecs.getTotalNum();
  TwoBodyCoefs.resize(nTwo);
  TwoBody_rhoG.resize(nTwo);
  TwoBodyPhase.resize(nTwo);
  TwoBody_e2iGr_new.resize(nTwo);
  TwoBody_e2iGr_old.resize(nTwo);
  Delta_e2iGr.resize(nElecs,nTwo);
  for (int iat=0; iat<nElecs; iat++)
    for (int i=0; i<nTwo; i++)
      Delta_e2iGr(iat,i) = ComplexType();
  // Set Ion_rhoG
  for (int i=0; i<OneBodyGvecs.size(); i++)
  {
    Ion_rhoG[0] = ComplexType();
    for (int iat=0; iat<ions.getTotalNum(); iat++)
    {
      double phase = dot(OneBodyGvecs[i],ions.R[iat]);
      Ion_rhoG[i] += ComplexType(std::cos(phase), std::sin(phase));
    }
  }
}

void
kSpaceJastrow::setCoefficients(std::vector<RealType> &oneBodyCoefs,
                               std::vector<RealType> &twoBodyCoefs)
{
  int kk(0);
//     for (int i=0; i<oneBodyCoefs.size(); i++)
//       std::cerr << "oneBodyCoefs[" << i << "] = " << oneBodyCoefs[i] << std::endl;
  if (oneBodyCoefs.size() != 2*OneBodySymmCoefs.size())
  {
    app_warning() << "Warning!  Wrong number of coefficients specified in "
                  << "kSpaceJastrow's one-body coefficients.\n"
                  << oneBodyCoefs.size() << " were specified.  Should have been "
                  << 2*OneBodySymmCoefs.size() << std::endl;
    for (int i=0; i<OneBodySymmCoefs.size(); i++)
    {
      OneBodySymmCoefs[i].cG = ComplexType();
      myVars[kk++]=0.0;
      myVars[kk++]=0.0;
      OneBodySymmCoefs[i].set (OneBodyCoefs);
    }
  }
  else
  {
    for (int i=0; i<OneBodySymmCoefs.size(); i++)
    {
      OneBodySymmCoefs[i].cG = ComplexType (oneBodyCoefs[2*i+0],oneBodyCoefs[2*i+1]);
      myVars[kk++]=oneBodyCoefs[2*i+0];
      myVars[kk++]=oneBodyCoefs[2*i+1];
      OneBodySymmCoefs[i].set (OneBodyCoefs);
    }
  }
  if (twoBodyCoefs.size() != TwoBodySymmCoefs.size())
  {
    app_warning() << "Warning!  Wrong number of coefficients specified in "
                  << "kSpaceJastrow's two-body coefficients.\n"
                  << twoBodyCoefs.size() << " were specified.  Should have been "
                  << TwoBodySymmCoefs.size() << std::endl;
    for (int i=0; i<TwoBodySymmCoefs.size(); i++)
    {
      TwoBodySymmCoefs[i].cG = 0.0;
      myVars[kk++]=0.0;
      TwoBodySymmCoefs[i].set (TwoBodyCoefs);
    }
  }
  else
  {
    for (int i=0; i<TwoBodySymmCoefs.size(); i++)
    {
      TwoBodySymmCoefs[i].cG = twoBodyCoefs[i];
      myVars[kk++]=twoBodyCoefs[i];
      TwoBodySymmCoefs[i].set (TwoBodyCoefs);
    }
  }
}

void
kSpaceJastrow::resetTargetParticleSet(ParticleSet& P)
{
  for (int i=0; i<TwoBodyGvecs.size(); i++)
  {
    TwoBody_rhoG[i] = ComplexType();
    for (int iat=0; iat<NumElecs; iat++)
    {
      double phase, s, c;
      phase = dot (TwoBodyGvecs[i], P.R[iat]);
      sincos(phase, &s, &c);
      TwoBody_rhoG[i] += ComplexType(c,s);
    }
  }
}

///////////////////////////////////////////////////////////////
//                  Evaluation functions                     //
///////////////////////////////////////////////////////////////

kSpaceJastrow::RealType
kSpaceJastrow::evaluateLog(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L)
{
  RealType J1(0.0), J2(0.0);
  int N = P.getTotalNum();
  ComplexType eye(0.0, 1.0);
  int nOne = OneBodyGvecs.size();
  for (int iat=0; iat<N; iat++)
  {
    PosType r(P.R[iat]);
    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], r);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
    for (int i=0; i<nOne; i++)
    {
      ComplexType z = OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]);
      J1 += Prefactor*real(z);
      G[iat] += -Prefactor*real(z*eye)*OneBodyGvecs[i];
      L[iat] += -Prefactor*dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z);
    }
  }
  // Do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i=0; i<nTwo; i++)
    TwoBody_rhoG[i] = ComplexType();
  for (int iat=0; iat<N; iat++)
  {
    PosType r(P.R[iat]);
    for (int iG=0; iG<nTwo; iG++)
      TwoBodyPhase[iG] = dot(TwoBodyGvecs[iG], r);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
    for (int iG=0; iG<nTwo; iG++)
      TwoBody_rhoG[iG] += TwoBody_e2iGr_new[iG];
  }
//     std::cerr << "TwoBody_rhoG = ";
//     for (int i=0; i<nTwo; i++)
//       std::cerr << TwoBody_rhoG[i]  << "  ";
//     std::cerr << std::endl;
  for (int i=0; i<nTwo; i++)
    J2 += Prefactor*TwoBodyCoefs[i]*norm(TwoBody_rhoG[i]);
  for (int iat=0; iat<N; iat++)
  {
    PosType r(P.R[iat]);
    for (int i=0; i<nTwo; i++)
      TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
    for (int i=0; i<nTwo; i++)
    {
      PosType Gvec(TwoBodyGvecs[i]);
      ComplexType z = TwoBody_e2iGr_new[i];
      G[iat] += -Prefactor*2.0*Gvec*TwoBodyCoefs[i]*imag(qmcplusplus::conj(TwoBody_rhoG[i])*z);
      L[iat] += Prefactor*2.0*TwoBodyCoefs[i]*dot(Gvec,Gvec)*(-real(z*qmcplusplus::conj(TwoBody_rhoG[i])) + 1.0);
    }
  }
  return J1 + J2;
}


kSpaceJastrow::GradType kSpaceJastrow::evalGrad(ParticleSet& P, int iat)
{
//     RealType J1(0.0), J2(0.0);
  kSpaceJastrow::GradType G;
  int N = P.getTotalNum();
  ComplexType eye(0.0, 1.0);
  int nOne = OneBodyGvecs.size();
//     for (int iat=0; iat<N; iat++) {
  PosType r(P.R[iat]);
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], r);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
  for (int i=0; i<nOne; i++)
  {
    ComplexType z = OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]);
//    J1 += Prefactor*real(z);
    G += -Prefactor*real(z*eye)*OneBodyGvecs[i];
//    L[iat] += -Prefactor*dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z);
  }
//     }
  // Do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i=0; i<nTwo; i++)
    TwoBody_rhoG[i] = ComplexType();
  for (int iat2=0; iat2<N; iat2++)
  {
    PosType rp(P.R[iat2]);
    for (int iG=0; iG<nTwo; iG++)
      TwoBodyPhase[iG] = dot(TwoBodyGvecs[iG], rp);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
    for (int iG=0; iG<nTwo; iG++)
      TwoBody_rhoG[iG] += TwoBody_e2iGr_new[iG];
  }
//     std::cerr << "TwoBody_rhoG = ";
//     for (int i=0; i<nTwo; i++)
//       std::cerr << TwoBody_rhoG[i]  << "  ";
//     std::cerr << std::endl;
//     for (int i=0; i<nTwo; i++)
//       J2 += Prefactor*TwoBodyCoefs[i]*norm(TwoBody_rhoG[i]);
//     for (int iat=0; iat<N; iat++) {
//       PosType r(P.R[iat]);
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
  for (int i=0; i<nTwo; i++)
  {
    PosType Gvec(TwoBodyGvecs[i]);
    ComplexType z = TwoBody_e2iGr_new[i];
    G += -Prefactor*2.0*Gvec*TwoBodyCoefs[i]*imag(qmcplusplus::conj(TwoBody_rhoG[i])*z);
//    L[iat] += Prefactor*2.0*TwoBodyCoefs[i]*dot(Gvec,Gvec)*(-real(z*qmcplusplus::conj(TwoBody_rhoG[i])) + 1.0);
  }
//     }
  return G;
}

kSpaceJastrow::ValueType
kSpaceJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  ComplexType eye(0.0, 1.0);
  RealType J1new(0.0), J1old(0.0), J2new(0.0), J2old(0.0);
  const PosType &rnew(P.activePos), &rold(P.R[iat]);
  // Compute one-body contribution
  int nOne = OneBodyGvecs.size();
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
  for (int i=0; i<nOne; i++)
  {
    ComplexType z = OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]);
    J1new += Prefactor*real(z);
    grad_iat += -Prefactor*real(z*eye)*OneBodyGvecs[i];
  }
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
  for (int i=0; i<nOne; i++)
    J1old += Prefactor*real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  // Now, do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_old);
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
  for (int i=0; i<nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i];
    J2old += Prefactor*TwoBodyCoefs[i]*std::norm(rho_G);
  }
  for (int i=0; i<nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
    J2new += Prefactor*TwoBodyCoefs[i]*std::norm(rho_G);
  }
  for (int i=0; i<nTwo; i++)
  {
    PosType Gvec(TwoBodyGvecs[i]);
    ComplexType z = TwoBody_e2iGr_new[i];
    grad_iat += -Prefactor*2.0*TwoBodyGvecs[i]*TwoBodyCoefs[i]*imag(qmcplusplus::conj(TwoBody_rhoG[i])*TwoBody_e2iGr_new[i]);
  }
  return std::exp(J1new+J2new - (J1old + J2old));
}

/* evaluate the ratio with P.R[iat]
 *
 */
kSpaceJastrow::ValueType
kSpaceJastrow::ratio(ParticleSet& P, int iat)
{
  RealType J1new(0.0), J1old(0.0), J2new(0.0), J2old(0.0);
  const PosType &rnew(P.activePos), &rold(P.R[iat]);
  // Compute one-body contribution
  int nOne = OneBodyGvecs.size();
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
  for (int i=0; i<nOne; i++)
    J1new += Prefactor*real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
  for (int i=0; i<nOne; i++)
    J1old += Prefactor*real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  // Now, do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_old);
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
  for (int i=0; i<nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i];
    J2old += Prefactor*TwoBodyCoefs[i]*std::norm(rho_G);
  }
  for (int i=0; i<nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
    J2new += Prefactor*TwoBodyCoefs[i]*std::norm(rho_G);
  }
  return std::exp(J1new+J2new - (J1old + J2old));
}

/** evaluate the ratio
*/
inline void kSpaceJastrow::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<kSpaceJastrow::ValueType>& ratios)
{
  RealType J1new(0.0);
  const PosType &rnew(P.activePos);
//     Compute one-body contribution
  int nOne = OneBodyGvecs.size();
  for (int i=0; i<nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
  eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
// //
  for (int i=0; i<nOne; i++)
    J1new += Prefactor*real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  // Now, do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i=0; i<nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
  eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
  int N = P.getTotalNum();
  for (int n=0; n<N; n++)
  {
    RealType J1old(0.0), J2Rat(0.0);
    const PosType &rold(P.R[n]);
    for (int i=0; i<nOne; i++)
      OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
    eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
    for (int i=0; i<nOne; i++)
      J1old += Prefactor*real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
    for (int i=0; i<nTwo; i++)
      TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
    eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_old);
    for (int i=0; i<nTwo; i++)
    {
      ComplexType rho_G_new = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
      ComplexType rho_G_old = TwoBody_rhoG[i];
      J2Rat += Prefactor*TwoBodyCoefs[i]*(std::norm(rho_G_new) - std::norm(rho_G_old));
    }
    ratios[n]=std::exp(J1new-J1old + J2Rat);
  }
}


void
kSpaceJastrow::restore(int iat)
{
  //substract the addition in logRatio
  //if(NeedToRestore) Rhok -= delta_eikr;
}

void
kSpaceJastrow::acceptMove(ParticleSet& P, int iat)
{
  for (int i=0; i<TwoBody_e2iGr_new.size(); i++)
    TwoBody_rhoG[i] += Delta_e2iGr(iat,i);
  //TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
  // copy(eikr_new.data(),eikr_new.data()+MaxK,eikr[iat]);
  // U += offU;
  // dU += offdU;
  // d2U += offd2U;
}

void
kSpaceJastrow::registerData(ParticleSet& P, WFBufferType& buf)
{
  LogValue=evaluateLog(P,P.G,P.L);
  // eikr.resize(NumPtcls,MaxK);
  // eikr_new.resize(MaxK);
  // delta_eikr.resize(MaxK);
  // for(int iat=0; iat<NumPtcls; iat++)
  //   copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
  // buf.add(Rhok.first_address(), Rhok.last_address());
  // buf.add(U.first_address(), U.last_address());
  // buf.add(d2U.first_address(), d2U.last_address());
  // buf.add(FirstAddressOfdU,LastAddressOfdU);
  // return LogValue;
}

kSpaceJastrow::RealType
kSpaceJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf,
                            bool fromscratch)
{
  LogValue=evaluateLog(P,P.G,P.L);
  // for(int iat=0; iat<NumPtcls; iat++)
  //   copy(P.SK->eikr[iat],P.SK->eikr[iat]+MaxK,eikr[iat]);
  // buf.put(Rhok.first_address(), Rhok.last_address());
  // buf.put(U.first_address(), U.last_address());
  // buf.put(d2U.first_address(), d2U.last_address());
  // buf.put(FirstAddressOfdU,LastAddressOfdU);
  // return LogValue;
  return LogValue;
}

void
kSpaceJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  for (int i=0; i<TwoBodyCoefs.size(); i++)
    TwoBody_rhoG[i] = ComplexType();
  for (int iat=0; iat<NumElecs; iat++)
  {
    for (int i=0; i<TwoBodyCoefs.size(); i++)
      TwoBodyPhase[i] = dot(TwoBodyGvecs[i],P.R[iat]);
    eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_new);
    for (int i=0; i<TwoBodyCoefs.size(); i++)
      TwoBody_rhoG[i] += TwoBody_e2iGr_new[i];
  }
}

void kSpaceJastrow::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
  int nOne = OneBodyGvecs.size();
  int obi=0;
  if (nOne)
  {
    OneBodyVarMap.resize(nOne);
    for (int i=0; i<nOne; i++)
    {
      //two coeffs for each of these points, imaginary coefficients.
      OneBodyVarMap[i]=obi;
      if (i==OneBodySymmCoefs[obi/2].lastIndex)
        obi+=2;
    }
  }
  int nTwo = TwoBodyGvecs.size();
  TwoBodyVarMap.resize(nTwo);
  int tbi=0;
  for (int i=0; i<nTwo; i++)
  {
    //one coeff for each of these points, real coefficients.
    TwoBodyVarMap[i]=obi+tbi;
    if (i==TwoBodySymmCoefs[tbi].lastIndex)
      tbi+=1;
  }
}

void kSpaceJastrow::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
  //Optimizable=myVars.is_optimizable();
}

void kSpaceJastrow::reportStatus(std::ostream& os)
{
}

void
kSpaceJastrow::resetParameters(const opt_variables_type& active)
{
  int ii=0;
  
  //Update the one body coefficients.
  for (int i=0; i<OneBodySymmCoefs.size(); i++)
  {
    //Reset and update the real part of one body term.

    int loc_r=myVars.where(ii);
    if(loc_r>=0)
    {
      myVars[ii]=active[loc_r];    //update the optimization parameter
      //lvalue error with LLVM
      OneBodySymmCoefs[i].cG=ComplexType(myVars[ii],OneBodySymmCoefs[i].cG.imag());  
      //OneBodySymmCoefs[i].cG.real()=myVars[ii];  //update the coefficient from local opt parametr
      ii++;  
    }
    //The imaginary part...
    int loc_i=myVars.where(ii);
    if(loc_i>=0)
    {
      myVars[ii]=active[loc_i];
      //lvalue error with LLVM 
      OneBodySymmCoefs[i].cG=ComplexType(OneBodySymmCoefs[i].cG.real(),myVars[ii]);  
      //OneBodySymmCoefs[i].cG.imag()=myVars[ii];
      ii++;
    }
  }

  //Update the two-body coefficients
  for (int i=0; i<TwoBodySymmCoefs.size(); i++)
  {
    int loc=myVars.where(ii);
    if(loc>=0)
    {
      myVars[ii]=active[loc];
      TwoBodySymmCoefs[i].cG=myVars[ii];
    }
    ii++;
  }


  for (int i=0; i<OneBodySymmCoefs.size(); i++)
    OneBodySymmCoefs[i].set(OneBodyCoefs);
  for (int i=0; i<TwoBodySymmCoefs.size(); i++)
    TwoBodySymmCoefs[i].set(TwoBodyCoefs);
}

bool
kSpaceJastrow::put(xmlNodePtr cur)
{
  return true;
}

OrbitalBasePtr kSpaceJastrow::makeClone(ParticleSet& tqp) const
{
  kSpaceJastrow *kj =new kSpaceJastrow(Ions,tqp);
  kj->copyFrom(*this);
  // kSpaceJastrow *kj = new kSpaceJastrow(*this);
  // kj->VarMap.clear();
  // for (int i=0; i<OneBodySymmCoefs.size(); i++) {
  //   std::stringstream name_real, name_imag;
  //   name_real << OneBodyID << "_" << 2*i;
  //   name_imag << OneBodyID << "_" << 2*i+1;
  //   kj->VarMap[name_real.str()] = &(kj->OneBodySymmCoefs[i].cG.real());
  //   kj->VarMap[name_imag.str()] = &(kj->OneBodySymmCoefs[i].cG.imag());
  // }
  // for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
  //   std::stringstream name;
  //   name << TwoBodyID << "_" << i;
  //   kj->VarMap[name.str()] = &(kj->TwoBodySymmCoefs[i].cG);
  // }
  return kj;
}

/** constructor to initialize Ions and Elecs
 */
kSpaceJastrow::kSpaceJastrow(const ParticleSet& ions, ParticleSet& els):
  Ions(ions), Elecs(els)
{
}

void kSpaceJastrow::copyFrom(const kSpaceJastrow& old)
{
  CellVolume=old.CellVolume;
  NormConstant=old.NormConstant;
  NumElecs=old.NumElecs;
  NumSpins=old.NumSpins;
  NumIons=old.NumIons;
  NumIonSpecies=old.NumIonSpecies;
  OneBodyGvecs=old.OneBodyGvecs;
  TwoBodyGvecs=old.TwoBodyGvecs;
  OneBodySymmCoefs=old.OneBodySymmCoefs;
  TwoBodySymmCoefs=old.TwoBodySymmCoefs;
  OneBodyCoefs=old.OneBodyCoefs;
  TwoBodyCoefs=old.TwoBodyCoefs;
  OneBodySymmType=old.OneBodySymmType;
  TwoBodySymmType=old.TwoBodySymmType;
  Ion_rhoG=old.Ion_rhoG;
  OneBody_rhoG=old.OneBody_rhoG;
  TwoBody_rhoG=old.TwoBody_rhoG;
  OneBodyPhase=old.OneBodyPhase;
  TwoBodyPhase=old.TwoBodyPhase;
  OneBody_e2iGr=old.OneBody_e2iGr;
  TwoBody_e2iGr_new=old.TwoBody_e2iGr_new;
  TwoBody_e2iGr_old=old.TwoBody_e2iGr_old;
  Delta_e2iGr=old.Delta_e2iGr;
  OneBodyID=old.OneBodyID;
  TwoBodyID=old.TwoBodyID;
  //copy the variable map
  myVars=old.myVars;
  Optimizable=old.Optimizable;
  TwoBodyVarMap=old.TwoBodyVarMap;
  OneBodyVarMap=old.OneBodyVarMap;
  Prefactor=old.Prefactor;
  //for (int i=0; i<OneBodySymmCoefs.size(); i++) {
  //  std::stringstream name_real, name_imag;
  //  name_real << OneBodyID << "_" << 2*i;
  //  name_imag << OneBodyID << "_" << 2*i+1;
  //  VarMap[name_real.str()] = &(OneBodySymmCoefs[i].cG.real());
  //  VarMap[name_imag.str()] = &(OneBodySymmCoefs[i].cG.imag());
  //}
  //for (int i=0; i<TwoBodySymmCoefs.size(); i++) {
  //  std::stringstream name;
  //  name << TwoBodyID << "_" << i;
  //  VarMap[name.str()] = &(TwoBodySymmCoefs[i].cG);
  //}
}

void kSpaceJastrow::evaluateDerivatives(ParticleSet& P,
                                        const opt_variables_type& active,
                                        std::vector<RealType>& dlogpsi,
                                        std::vector<RealType>& dhpsioverpsi)
{
  bool recalculate(false);
  for (int k=0; k<myVars.size(); ++k)
  {
    int kk=myVars.where(k);
    if (kk<0)
      continue;
    if (active.recompute(kk))
      recalculate=true;
  }
  if (recalculate)
  {
    int N = P.getTotalNum();
    ComplexType eye(0.0, 1.0);
    RealType tmp_dot;
    int nOne = OneBodyGvecs.size();
    if (nOne)
    {
      for (int iat=0; iat<N; iat++)
      {
        PosType r(P.R[iat]);
        for (int i=0; i<nOne; i++)
          OneBodyPhase[i] = dot(OneBodyGvecs[i], r);
        eval_e2iphi (OneBodyPhase, OneBody_e2iGr);
        for (int i=0; i<nOne; i++)
        {
          ComplexType z =  qmcplusplus::conj(OneBody_e2iGr[i]);
          int kk=myVars.where(OneBodyVarMap[i]);
          if (kk>=0)
          {
            //real part of coeff
            dlogpsi[kk] += Prefactor*real(z);
            //convert(dot(OneBodyGvecs[i],P.G[iat]),tmp_dot);
            convert(dot(P.G[iat],OneBodyGvecs[i]),tmp_dot);
            dhpsioverpsi[kk] +=  0.5*Prefactor*dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(z) + Prefactor*real(z*eye)*tmp_dot;
            //	+ Prefactor*real(z*eye)*real(dot(OneBodyGvecs[i],P.G[iat]));
            //imaginary part of coeff,
            dlogpsi[kk+1] += Prefactor*real(eye*z);
            //mius here due to i*i term
            //dhpsioverpsi[kk+1] += 0.5*Prefactor*dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(eye*z) - Prefactor*real(z)*real(dot(OneBodyGvecs[i],P.G[iat]));
            dhpsioverpsi[kk+1] += 0.5*Prefactor*dot(OneBodyGvecs[i],OneBodyGvecs[i])*real(eye*z) - Prefactor*real(z)*tmp_dot;
          }
        }
      }
    }
    // Do two-body part
    int nTwo = TwoBodyGvecs.size();
    for (int i=0; i<nTwo; i++)
      TwoBody_rhoG[i] = ComplexType();
    for (int iat=0; iat<N; iat++)
    {
      PosType r(P.R[iat]);
      for (int iG=0; iG<nTwo; iG++)
        TwoBodyPhase[iG] = dot(TwoBodyGvecs[iG], r);
      eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
      for (int iG=0; iG<nTwo; iG++)
        TwoBody_rhoG[iG] += TwoBody_e2iGr_new[iG];
    }
    for (int i=0; i<nTwo; i++)
    {
      int kk=myVars.where(TwoBodyVarMap[i]);
      if (kk>=0)
      {
        dlogpsi[kk] += Prefactor*norm(TwoBody_rhoG[i]);
      }
    }
    for (int iat=0; iat<N; iat++)
    {
      PosType r(P.R[iat]);
      for (int i=0; i<nTwo; i++)
        TwoBodyPhase[i] = dot(TwoBodyGvecs[i], r);
      eval_e2iphi (TwoBodyPhase, TwoBody_e2iGr_new);
      for (int i=0; i<nTwo; i++)
      {
        PosType Gvec(TwoBodyGvecs[i]);
        ComplexType z = TwoBody_e2iGr_new[i];
        int kk=myVars.where(TwoBodyVarMap[i]);
        if (kk>0)
        {
          convert(dot(P.G[iat],Gvec),tmp_dot);
          //dhpsioverpsi[kk] -= Prefactor*dot(Gvec,Gvec)*(-real(z*qmcplusplus::conj(TwoBody_rhoG[i])) + 1.0) - Prefactor*2.0*real(dot(P.G[iat],Gvec))*imag(qmcplusplus::conj(TwoBody_rhoG[i])*z);
          dhpsioverpsi[kk] -= Prefactor*dot(Gvec,Gvec)*(-real(z*qmcplusplus::conj(TwoBody_rhoG[i])) + 1.0) - Prefactor*2.0*tmp_dot*imag(qmcplusplus::conj(TwoBody_rhoG[i])*z);
        }
      }
    }
  }
}

}

