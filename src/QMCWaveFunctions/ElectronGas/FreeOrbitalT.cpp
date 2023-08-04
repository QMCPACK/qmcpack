//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//                    William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "FreeOrbitalT.h"

namespace qmcplusplus
{

template<class T>
FreeOrbitalT<T>::FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart) : SPOSetT<T>(my_name)
{}

//Explicit template specialization
template<>
FreeOrbitalT<float>::FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart)
    : SPOSetT<float>(my_name),
      kvecs(kpts_cart),
      mink(1), // treat k=0 as special case
      maxk(kpts_cart.size()),
      k2neg(maxk)
{
  this->OrbitalSetSize = 2 * maxk - 1; // k=0 has no (cos, sin) split, SPOSet member
  for (int ik = 0; ik < maxk; ik++)
    k2neg[ik] = -dot(kvecs[ik], kvecs[ik]);
}

template<>
FreeOrbitalT<double>::FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart)
    : SPOSetT<double>(my_name),
      kvecs(kpts_cart),
      mink(1), // treat k=0 as special case
      maxk(kpts_cart.size()),
      k2neg(maxk)
{
  this->OrbitalSetSize = 2 * maxk - 1; // k=0 has no (cos, sin) split, SPOSet member
  for (int ik = 0; ik < maxk; ik++)
    k2neg[ik] = -dot(kvecs[ik], kvecs[ik]);
}

template<>
FreeOrbitalT<std::complex<float>>::FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart)
    : SPOSetT<std::complex<float>>(my_name),
      kvecs(kpts_cart),
      mink(0), // treat k=0 as special case
      maxk(kpts_cart.size()),
      k2neg(maxk)
{
  this->OrbitalSetSize = maxk; // SPOSet member
  for (int ik = 0; ik < maxk; ik++)
    k2neg[ik] = -dot(kvecs[ik], kvecs[ik]);
}

template<>
FreeOrbitalT<std::complex<double>>::FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart)
    : SPOSetT<std::complex<double>>(my_name),
      kvecs(kpts_cart),
      mink(0), // treat k=0 as special case
      maxk(kpts_cart.size()),
      k2neg(maxk)
{
  this->OrbitalSetSize = maxk; // SPOSet member
  for (int ik = 0; ik < maxk; ik++)
    k2neg[ik] = -dot(kvecs[ik], kvecs[ik]);
}


template<class T>
void FreeOrbitalT<T>::evaluateVGL(const ParticleSet& P,
                                  int iat,
                                  ValueVector& pvec,
                                  GradVector& dpvec,
                                  ValueVector& d2pvec)
{}

template<>
void FreeOrbitalT<float>::evaluateVGL(const ParticleSet& P,
                                      int iat,
                                      ValueVector& pvec,
                                      GradVector& dpvec,
                                      ValueVector& d2pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
    dpvec[j1]    = -sinkr * kvecs[ik];
    dpvec[j2]    = coskr * kvecs[ik];
    d2pvec[j1]   = k2neg[ik] * coskr;
    d2pvec[j2]   = k2neg[ik] * sinkr;
  }
  pvec[0]   = 1.0;
  dpvec[0]  = 0.0;
  d2pvec[0] = 0.0;
}

template<>
void FreeOrbitalT<double>::evaluateVGL(const ParticleSet& P,
                                       int iat,
                                       ValueVector& pvec,
                                       GradVector& dpvec,
                                       ValueVector& d2pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
    dpvec[j1]    = -sinkr * kvecs[ik];
    dpvec[j2]    = coskr * kvecs[ik];
    d2pvec[j1]   = k2neg[ik] * coskr;
    d2pvec[j2]   = k2neg[ik] * sinkr;
  }
  pvec[0]   = 1.0;
  dpvec[0]  = 0.0;
  d2pvec[0] = 0.0;
}


template<>
void FreeOrbitalT<std::complex<float>>::evaluateVGL(const ParticleSet& P,
                                                    int iat,
                                                    ValueVector& pvec,
                                                    GradVector& dpvec,
                                                    ValueVector& d2pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);

    pvec[ik]   = ValueType(coskr, sinkr);
    dpvec[ik]  = ValueType(-sinkr, coskr) * kvecs[ik];
    d2pvec[ik] = ValueType(k2neg[ik] * coskr, k2neg[ik] * sinkr);
  }
}

template<>
void FreeOrbitalT<std::complex<double>>::evaluateVGL(const ParticleSet& P,
                                                     int iat,
                                                     ValueVector& pvec,
                                                     GradVector& dpvec,
                                                     ValueVector& d2pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);

    pvec[ik]   = ValueType(coskr, sinkr);
    dpvec[ik]  = ValueType(-sinkr, coskr) * kvecs[ik];
    d2pvec[ik] = ValueType(k2neg[ik] * coskr, k2neg[ik] * sinkr);
  }
}


template<>
void FreeOrbitalT<float>::evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
  }
  pvec[0] = 1.0;
}

template<>
void FreeOrbitalT<double>::evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
  }
  pvec[0] = 1.0;
}

template<>
void FreeOrbitalT<std::complex<float>>::evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);

    pvec[ik]     = std::complex<float>(coskr, sinkr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
  }
}

template<>
void FreeOrbitalT<std::complex<double>>::evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec)
{
  const PosType& r = P.activeR(iat);
  RealType sinkr, coskr;
  for (int ik = mink; ik < maxk; ik++)
  {
    sincos(dot(kvecs[ik], r), &sinkr, &coskr);

    pvec[ik]     = std::complex<double>(coskr, sinkr);
    const int j2 = 2 * ik;
    const int j1 = j2 - 1;
    pvec[j1]     = coskr;
    pvec[j2]     = sinkr;
  }
}

template<class T>
void FreeOrbitalT<T>::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& phi,
                                           GradMatrix& dphi,
                                           HessMatrix& d2phi_mat)
{}


template<>
void FreeOrbitalT<float>::evaluate_notranspose(const ParticleSet& P,
                                               int first,
                                               int last,
                                               ValueMatrix& phi,
                                               GradMatrix& dphi,
                                               HessMatrix& d2phi_mat)
{
  RealType sinkr, coskr;
  float phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], this->OrbitalSetSize);
    GradVector dp(dphi[i], this->OrbitalSetSize);
    HessVector hess(d2phi_mat[i], this->OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const int j2 = 2 * ik;
      const int j1 = j2 - 1;
      p[j1]        = coskr;
      p[j2]        = sinkr;
      dp[j1]       = -sinkr * kvecs[ik];
      dp[j2]       = coskr * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[j1](la, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la];
        hess[j2](la, la) = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[j1](la, lb) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j2](la, lb) = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j1](lb, la) = hess[j1](la, lb);
          hess[j2](lb, la) = hess[j2](la, lb);
        }
      }
    }
    p[0]    = 1.0;
    dp[0]   = 0.0;
    hess[0] = 0.0;
  }
}

template<>
void FreeOrbitalT<double>::evaluate_notranspose(const ParticleSet& P,
                                                int first,
                                                int last,
                                                ValueMatrix& phi,
                                                GradMatrix& dphi,
                                                HessMatrix& d2phi_mat)
{
  RealType sinkr, coskr;
  double phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], this->OrbitalSetSize);
    GradVector dp(dphi[i], this->OrbitalSetSize);
    HessVector hess(d2phi_mat[i], this->OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const int j2 = 2 * ik;
      const int j1 = j2 - 1;
      p[j1]        = coskr;
      p[j2]        = sinkr;
      dp[j1]       = -sinkr * kvecs[ik];
      dp[j2]       = coskr * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[j1](la, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la];
        hess[j2](la, la) = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[j1](la, lb) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j2](la, lb) = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j1](lb, la) = hess[j1](la, lb);
          hess[j2](lb, la) = hess[j2](la, lb);
        }
      }
    }
    p[0]    = 1.0;
    dp[0]   = 0.0;
    hess[0] = 0.0;
  }
}


template<>
void FreeOrbitalT<std::complex<float>>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             ValueMatrix& phi,
                                                             GradMatrix& dphi,
                                                             HessMatrix& d2phi_mat)
{
  RealType sinkr, coskr;
  std::complex<float> phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], this->OrbitalSetSize);
    GradVector dp(dphi[i], this->OrbitalSetSize);
    HessVector hess(d2phi_mat[i], this->OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);

      phi_of_r = std::complex<float>(coskr, sinkr);
      p[ik]    = phi_of_r;

      dp[ik] = std::complex<float>(-sinkr, coskr) * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[ik](la, la) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[ik](la, lb) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[ik](lb, la) = hess[ik](la, lb);
        }
      }
    }
  }
}

template<>
void FreeOrbitalT<std::complex<double>>::evaluate_notranspose(const ParticleSet& P,
                                                              int first,
                                                              int last,
                                                              ValueMatrix& phi,
                                                              GradMatrix& dphi,
                                                              HessMatrix& d2phi_mat)
{
  RealType sinkr, coskr;
  std::complex<double> phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], this->OrbitalSetSize);
    GradVector dp(dphi[i], this->OrbitalSetSize);
    HessVector hess(d2phi_mat[i], this->OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);

      phi_of_r = std::complex<double>(coskr, sinkr);
      p[ik]    = phi_of_r;

      dp[ik] = std::complex<double>(-sinkr, coskr) * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[ik](la, la) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[ik](la, lb) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[ik](lb, la) = hess[ik](la, lb);
        }
      }
    }
  }
}

template<class T>
void FreeOrbitalT<T>::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& phi,
                                           GradMatrix& dphi,
                                           HessMatrix& d2phi_mat,
                                           GGGMatrix& d3phi_mat)
{}

template<>
void FreeOrbitalT<float>::evaluate_notranspose(const ParticleSet& P,
                                               int first,
                                               int last,
                                               ValueMatrix& phi,
                                               GradMatrix& dphi,
                                               HessMatrix& d2phi_mat,
                                               GGGMatrix& d3phi_mat)
{
  RealType sinkr, coskr;
  ValueType phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], OrbitalSetSize);
    GradVector dp(dphi[i], OrbitalSetSize);
    HessVector hess(d2phi_mat[i], OrbitalSetSize);
    GGGVector ggg(d3phi_mat[i], OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const int j2 = 2 * ik;
      const int j1 = j2 - 1;
      p[j1]        = coskr;
      p[j2]        = sinkr;
      dp[j1]       = -sinkr * kvecs[ik];
      dp[j2]       = coskr * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[j1](la, la)    = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la];
        hess[j2](la, la)    = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[la];
        ggg[j1][la](la, la) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[la] * (kvecs[ik])[la];
        ggg[j2][la](la, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[j1](la, lb)    = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j2](la, lb)    = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j1](lb, la)    = hess[j1](la, lb);
          hess[j2](lb, la)    = hess[j2](la, lb);
          ggg[j1][la](lb, la) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[la];
          ggg[j2][la](lb, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[la];
          ggg[j1][la](la, lb) = ggg[j1][la](lb, la);
          ggg[j2][la](la, lb) = ggg[j2][la](lb, la);
          ggg[j1][lb](la, la) = ggg[j1][la](lb, la);
          ggg[j2][lb](la, la) = ggg[j2][la](lb, la);
          ggg[j1][la](lb, lb) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lb];
          ggg[j2][la](lb, lb) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lb];
          ggg[j1][lb](la, lb) = ggg[j1][la](lb, lb);
          ggg[j2][lb](la, lb) = ggg[j2][la](lb, lb);
          ggg[j1][lb](lb, la) = ggg[j1][la](lb, lb);
          ggg[j2][lb](lb, la) = ggg[j2][la](lb, lb);
          for (int lc = lb + 1; lc < OHMMS_DIM; lc++)
          {
            ggg[j1][la](lb, lc) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lc];
            ggg[j2][la](lb, lc) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lc];
            ggg[j1][la](lc, lb) = ggg[j1][la](lb, lc);
            ggg[j2][la](lc, lb) = ggg[j2][la](lb, lc);
            ggg[j1][lb](la, lc) = ggg[j1][la](lb, lc);
            ggg[j2][lb](la, lc) = ggg[j2][la](lb, lc);
            ggg[j1][lb](lc, la) = ggg[j1][la](lb, lc);
            ggg[j2][lb](lc, la) = ggg[j2][la](lb, lc);
            ggg[j1][lc](la, lb) = ggg[j1][la](lb, lc);
            ggg[j2][lc](la, lb) = ggg[j2][la](lb, lc);
            ggg[j1][lc](lb, la) = ggg[j1][la](lb, lc);
            ggg[j2][lc](lb, la) = ggg[j2][la](lb, lc);
          }
        }
      }
    }

    p[0]    = 1.0;
    dp[0]   = 0.0;
    hess[0] = 0.0;
    ggg[0]  = 0.0;
  }
}

template<>
void FreeOrbitalT<double>::evaluate_notranspose(const ParticleSet& P,
                                                int first,
                                                int last,
                                                ValueMatrix& phi,
                                                GradMatrix& dphi,
                                                HessMatrix& d2phi_mat,
                                                GGGMatrix& d3phi_mat)
{
  RealType sinkr, coskr;
  ValueType phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], OrbitalSetSize);
    GradVector dp(dphi[i], OrbitalSetSize);
    HessVector hess(d2phi_mat[i], OrbitalSetSize);
    GGGVector ggg(d3phi_mat[i], OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const int j2 = 2 * ik;
      const int j1 = j2 - 1;
      p[j1]        = coskr;
      p[j2]        = sinkr;
      dp[j1]       = -sinkr * kvecs[ik];
      dp[j2]       = coskr * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[j1](la, la)    = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la];
        hess[j2](la, la)    = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[la];
        ggg[j1][la](la, la) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[la] * (kvecs[ik])[la];
        ggg[j2][la](la, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[j1](la, lb)    = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j2](la, lb)    = -sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[j1](lb, la)    = hess[j1](la, lb);
          hess[j2](lb, la)    = hess[j2](la, lb);
          ggg[j1][la](lb, la) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[la];
          ggg[j2][la](lb, la) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[la];
          ggg[j1][la](la, lb) = ggg[j1][la](lb, la);
          ggg[j2][la](la, lb) = ggg[j2][la](lb, la);
          ggg[j1][lb](la, la) = ggg[j1][la](lb, la);
          ggg[j2][lb](la, la) = ggg[j2][la](lb, la);
          ggg[j1][la](lb, lb) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lb];
          ggg[j2][la](lb, lb) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lb];
          ggg[j1][lb](la, lb) = ggg[j1][la](lb, lb);
          ggg[j2][lb](la, lb) = ggg[j2][la](lb, lb);
          ggg[j1][lb](lb, la) = ggg[j1][la](lb, lb);
          ggg[j2][lb](lb, la) = ggg[j2][la](lb, lb);
          for (int lc = lb + 1; lc < OHMMS_DIM; lc++)
          {
            ggg[j1][la](lb, lc) = sinkr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lc];
            ggg[j2][la](lb, lc) = -coskr * (kvecs[ik])[la] * (kvecs[ik])[lb] * (kvecs[ik])[lc];
            ggg[j1][la](lc, lb) = ggg[j1][la](lb, lc);
            ggg[j2][la](lc, lb) = ggg[j2][la](lb, lc);
            ggg[j1][lb](la, lc) = ggg[j1][la](lb, lc);
            ggg[j2][lb](la, lc) = ggg[j2][la](lb, lc);
            ggg[j1][lb](lc, la) = ggg[j1][la](lb, lc);
            ggg[j2][lb](lc, la) = ggg[j2][la](lb, lc);
            ggg[j1][lc](la, lb) = ggg[j1][la](lb, lc);
            ggg[j2][lc](la, lb) = ggg[j2][la](lb, lc);
            ggg[j1][lc](lb, la) = ggg[j1][la](lb, lc);
            ggg[j2][lc](lb, la) = ggg[j2][la](lb, lc);
          }
        }
      }
    }

    p[0]    = 1.0;
    dp[0]   = 0.0;
    hess[0] = 0.0;
    ggg[0]  = 0.0;
  }
}

template<>
void FreeOrbitalT<std::complex<float>>::evaluate_notranspose(const ParticleSet& P,
                                                             int first,
                                                             int last,
                                                             ValueMatrix& phi,
                                                             GradMatrix& dphi,
                                                             HessMatrix& d2phi_mat,
                                                             GGGMatrix& d3phi_mat)
{
  RealType sinkr, coskr;
  ValueType phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], OrbitalSetSize);
    GradVector dp(dphi[i], OrbitalSetSize);
    HessVector hess(d2phi_mat[i], OrbitalSetSize);
    GGGVector ggg(d3phi_mat[i], OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const ValueType compi(0, 1);
      phi_of_r = ValueType(coskr, sinkr);
      p[ik]    = phi_of_r;
      dp[ik]   = compi * phi_of_r * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[ik](la, la) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[ik](la, lb) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[ik](lb, la) = hess[ik](la, lb);
        }
      }
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        ggg[ik][la] = compi * (kvecs[ik])[la] * hess[ik];
      }
    }
  }
}

template<>
void FreeOrbitalT<std::complex<double>>::evaluate_notranspose(const ParticleSet& P,
                                                              int first,
                                                              int last,
                                                              ValueMatrix& phi,
                                                              GradMatrix& dphi,
                                                              HessMatrix& d2phi_mat,
                                                              GGGMatrix& d3phi_mat)
{
  RealType sinkr, coskr;
  ValueType phi_of_r;
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], OrbitalSetSize);
    GradVector dp(dphi[i], OrbitalSetSize);
    HessVector hess(d2phi_mat[i], OrbitalSetSize);
    GGGVector ggg(d3phi_mat[i], OrbitalSetSize);

    const PosType& r = P.activeR(iat);
    for (int ik = mink; ik < maxk; ik++)
    {
      sincos(dot(kvecs[ik], r), &sinkr, &coskr);
      const ValueType compi(0, 1);
      phi_of_r = ValueType(coskr, sinkr);
      p[ik]    = phi_of_r;
      dp[ik]   = compi * phi_of_r * kvecs[ik];
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        hess[ik](la, la) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[la];
        for (int lb = la + 1; lb < OHMMS_DIM; lb++)
        {
          hess[ik](la, lb) = -phi_of_r * (kvecs[ik])[la] * (kvecs[ik])[lb];
          hess[ik](lb, la) = hess[ik](la, lb);
        }
      }
      for (int la = 0; la < OHMMS_DIM; la++)
      {
        ggg[ik][la] = compi * (kvecs[ik])[la] * hess[ik];
      }
    }
  }
}

// generic implementation

template<class T>
FreeOrbitalT<T>::~FreeOrbitalT()
{}

template<class T>
void FreeOrbitalT<T>::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& phi,
                                           GradMatrix& dphi,
                                           ValueMatrix& d2phi)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector p(phi[i], this->OrbitalSetSize);
    GradVector dp(dphi[i], this->OrbitalSetSize);
    ValueVector d2p(d2phi[i], this->OrbitalSetSize);
    evaluateVGL(P, iat, p, dp, d2p);
  }
}


template<class T>
void FreeOrbitalT<T>::report(const std::string& pad) const
{
  app_log() << pad << "FreeOrbital report" << std::endl;
  for (int ik = 0; ik < kvecs.size(); ik++)
  {
    app_log() << pad << ik << " " << kvecs[ik] << std::endl;
  }
  app_log() << pad << "end FreeOrbital report" << std::endl;
  app_log().flush();
}

template class FreeOrbitalT<float>;
template class FreeOrbitalT<double>;
template class FreeOrbitalT<std::complex<float>>;
template class FreeOrbitalT<std::complex<double>>;


} // namespace qmcplusplus
