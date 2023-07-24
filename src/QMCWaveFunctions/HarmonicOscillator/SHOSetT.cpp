//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SHOSetT.h"

#include "Utilities/string_utils.h"

namespace qmcplusplus
{
template<typename T>
SHOSetT<T>::SHOSetT(const std::string& my_name, RealType l, PosType c, const std::vector<SHOState*>& sho_states)
    : SPOSetT<T>(my_name), length(l), center(c)
{
  state_info.resize(sho_states.size());
  for (int s = 0; s < sho_states.size(); ++s)
    state_info[s] = *sho_states[s];
  initialize();
}

template<typename T>
void SHOSetT<T>::initialize()
{
  using std::sqrt;

  this->OrbitalSetSize = state_info.size();

  qn_max = -1;
  for (int s = 0; s < state_info.size(); ++s)
    for (int d = 0; d < QMCTraits::DIM; ++d)
      qn_max[d] = std::max(qn_max[d], state_info[s].quantum_number[d]);
  qn_max += 1;

  nmax = -1;
  for (int d = 0; d < QMCTraits::DIM; ++d)
    nmax = std::max(nmax, qn_max[d]);

  prefactors.resize(nmax);
  hermite.resize(QMCTraits::DIM, nmax);
  bvalues.resize(QMCTraits::DIM, nmax);

  if (nmax > 0)
  {
    prefactors[0] = 1.0 / (sqrt(sqrt(M_PI) * length));
    for (int n = 1; n < nmax; ++n)
      prefactors[n] = prefactors[n - 1] / sqrt(2. * n);
  }
}

template<typename T>
SHOSetT<T>::~SHOSetT() = default;

template<typename T>
std::unique_ptr<SPOSetT<T>> SHOSetT<T>::makeClone() const
{
  return std::make_unique<SHOSetT<T>>(*this);
}

template<typename T>
void SHOSetT<T>::report(const std::string& pad) const
{
  app_log() << pad << "SHOSet report" << std::endl;
  app_log() << pad << "  length    = " << length << std::endl;
  app_log() << pad << "  center    = " << center << std::endl;
  app_log() << pad << "  nmax      = " << nmax << std::endl;
  app_log() << pad << "  qn_max    = " << qn_max << std::endl;
  app_log() << pad << "  # states  = " << state_info.size() << std::endl;
  app_log() << pad << "  states" << std::endl;
  for (int s = 0; s < state_info.size(); ++s)
    state_info[s].sho_report(pad + "    " + int2string(s) + " ");
  app_log() << pad << "end SHOSet report" << std::endl;
  app_log().flush();
}

template<typename T>
void SHOSetT<T>::evaluateValue(const ParticleSetT<T>& P, int iat, ValueVector& psi)
{
  const PosType& r(P.activeR(iat));
  ValueVector p(&psi[0], this->size());
  evaluate_v(r, p);
}

template<typename T>
void SHOSetT<T>::evaluateVGL(const ParticleSetT<T>& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  const PosType& r(P.activeR(iat));
  ValueVector p(&psi[0], this->size());
  GradVector dp(&dpsi[0], this->size());
  ValueVector d2p(&d2psi[0], this->size());
  evaluate_vgl(r, p, dp, d2p);
}

template<typename T>
void SHOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P,
                                      int first,
                                      int last,
                                      ValueMatrix& logdet,
                                      GradMatrix& dlogdet,
                                      ValueMatrix& d2logdet)
{
  for (int iat = first, i = 0; iat < last; ++iat, ++i)
  {
    ValueVector p(logdet[i], this->size());
    GradVector dp(dlogdet[i], this->size());
    ValueVector d2p(d2logdet[i], this->size());
    evaluate_vgl(P.R[iat], p, dp, d2p);
  }
}

template<typename T>
void SHOSetT<T>::evaluate_v(PosType r, ValueVector& psi)
{
  PosType x = (r - center) / length;
  evaluate_hermite(x);
  evaluate_d0(x, psi);
}

template<typename T>
void SHOSetT<T>::evaluate_vgl(PosType r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  PosType x = (r - center) / length;
  evaluate_hermite(x);
  evaluate_d0(x, psi);
  evaluate_d1(x, psi, dpsi);
  evaluate_d2(x, psi, d2psi);
}

template<typename T>
void SHOSetT<T>::evaluate_hermite(const PosType& xpos)
{
  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    int nh = qn_max[d];
    if (nh > 0)
    {
      RealType x    = xpos[d];
      hermite(d, 0) = 1.0;
      RealType Hnm2 = 0.0;
      RealType Hnm1 = 1.0;
      for (int n = 1; n < nh; ++n)
      {
        RealType Hn   = 2 * (x * Hnm1 - (n - 1) * Hnm2);
        hermite(d, n) = Hn;
        Hnm2          = Hnm1;
        Hnm1          = Hn;
      }
    }
  }
}

template<typename T>
void SHOSetT<T>::evaluate_d0(const PosType& xpos, ValueVector& psi)
{
  using std::exp;
  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    RealType x = xpos[d];
    RealType g = exp(-.5 * x * x);
    for (int n = 0; n < qn_max[d]; ++n)
    {
      bvalues(d, n) = prefactors[n] * g * hermite(d, n);
    }
  }
  for (int s = 0; s < state_info.size(); ++s)
  {
    const SHOState& state = state_info[s];
    RealType phi          = 1.0;
    for (int d = 0; d < QMCTraits::DIM; ++d)
      phi *= bvalues(d, state.quantum_number[d]);
    psi[s] = phi;
  }
}

template<typename T>
void SHOSetT<T>::evaluate_d1(const PosType& xpos, ValueVector& psi, GradVector& dpsi)
{
  RealType ol = 1.0 / length;
  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    RealType x    = xpos[d];
    RealType Hnm1 = 0.0;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      RealType Hn   = hermite(d, n);
      bvalues(d, n) = (-x + 2 * n * Hnm1 / Hn) * ol;
      Hnm1          = Hn;
    }
  }
  for (int s = 0; s < state_info.size(); ++s)
  {
    const SHOState& state = state_info[s];
    TinyVector<T, QMCTraits::DIM> dphi;
    for (int d = 0; d < QMCTraits::DIM; ++d)
      dphi[d] = bvalues(d, state.quantum_number[d]);
    dphi *= psi[s];
    dpsi[s] = dphi;
  }
}

template<typename T>
void SHOSetT<T>::evaluate_d2(const PosType& xpos, ValueVector& psi, ValueVector& d2psi)
{
  RealType ol2 = 1.0 / (length * length);
  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    RealType x  = xpos[d];
    RealType x2 = x * x;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      bvalues(d, n) = (-1.0 + x2 - 2 * n) * ol2;
    }
  }
  for (int s = 0; s < state_info.size(); ++s)
  {
    const SHOState& state = state_info[s];
    T d2phi               = 0.0;
    for (int d = 0; d < QMCTraits::DIM; ++d)
      d2phi += bvalues(d, state.quantum_number[d]);
    d2phi *= psi[s];
    d2psi[s] = d2phi;
  }
}

template<typename T>
void SHOSetT<T>::evaluate_check(PosType r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  using std::exp;
  using std::sqrt;

  evaluate_vgl(r, psi, dpsi, d2psi);

  const int N = 6;
  RealType H[N], dH[N], d2H[N], pre[N];
  RealType p[N], dp[N], d2p[N];

  pre[0] = 1.0 / (sqrt(sqrt(M_PI) * length));
  for (int n = 1; n < N; ++n)
    pre[n] = pre[n - 1] / sqrt(2. * n);

  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    RealType x  = (r[d] - center[d]) / length;
    RealType x2 = x * x, x3 = x * x * x, x4 = x * x * x * x, x5 = x * x * x * x * x;
    H[0]       = 1;
    dH[0]      = 0;
    d2H[0]     = 0;
    H[1]       = 2 * x;
    dH[1]      = 2;
    d2H[1]     = 0;
    H[2]       = 4 * x2 - 2;
    dH[2]      = 8 * x;
    d2H[2]     = 8;
    H[3]       = 8 * x3 - 12 * x;
    dH[3]      = 24 * x2 - 12;
    d2H[3]     = 48 * x;
    H[4]       = 16 * x4 - 48 * x2 + 12;
    dH[4]      = 64 * x3 - 96 * x;
    d2H[4]     = 192 * x2 - 96;
    H[5]       = 32 * x5 - 160 * x3 + 120 * x;
    dH[5]      = 160 * x4 - 480 * x2 + 120;
    d2H[5]     = 640 * x3 - 960 * x;
    RealType g = exp(-x2 / 2);
    for (int n = 0; n < N; ++n)
    {
      p[n]   = pre[n] * g * H[n];
      dp[n]  = pre[n] * g * (-x * H[n] + dH[n]);
      d2p[n] = pre[n] * g * ((x2 - 1) * H[n] - 2 * x * dH[n] + d2H[n]);
    }
    app_log() << "eval check dim = " << d << "  x = " << x << std::endl;
    app_log() << "  hermite check" << std::endl;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      app_log() << "    " << n << " " << H[n] << std::endl;
      app_log() << "    " << n << " " << hermite(d, n) << std::endl;
    }
    app_log() << "  phi d0 check" << std::endl;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      app_log() << "    " << n << " " << p[n] << std::endl;
      app_log() << "    " << n << " " << d0_values(d, n) << std::endl;
    }
    app_log() << "  phi d1 check" << std::endl;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      app_log() << "    " << n << " " << dp[n] / p[n] << std::endl;
      app_log() << "    " << n << " " << d1_values(d, n) << std::endl;
    }
    app_log() << "  phi d2 check" << std::endl;
    for (int n = 0; n < qn_max[d]; ++n)
    {
      app_log() << "    " << n << " " << d2p[n] / p[n] << std::endl;
      app_log() << "    " << n << " " << d2_values(d, n) << std::endl;
    }
  }
}

template<typename T>
void SHOSetT<T>::test_derivatives()
{
  int n       = 3;
  PosType c   = 5.123;
  PosType L   = 1.0;
  PosType drg = L / n;
  PosType dr  = L / 1000;
  int nphi    = state_info.size();

  PosType o2dr, odr2;

  ValueVector vpsi, vpsitmp;
  GradVector vdpsi, vdpsin;
  ValueVector vd2psi, vd2psin;

  vpsi.resize(nphi);
  vdpsi.resize(nphi);
  vd2psi.resize(nphi);

  vpsitmp.resize(nphi);
  vdpsin.resize(nphi);
  vd2psin.resize(nphi);

  ValueVector psi(&vpsi[0], this->size());
  GradVector dpsi(&vdpsi[0], this->size());
  ValueVector d2psi(&vd2psi[0], this->size());

  ValueVector psitmp(&vpsitmp[0], this->size());
  GradVector dpsin(&vdpsin[0], this->size());
  ValueVector d2psin(&vd2psin[0], this->size());

  app_log() << " loading dr" << std::endl;

  RealType odr2sum = 0.0;
  for (int d = 0; d < QMCTraits::DIM; ++d)
  {
    RealType odr = 1.0 / dr[d];
    o2dr[d]      = .5 * odr;
    odr2[d]      = odr * odr;
    odr2sum += odr2[d];
  }

  app_log() << "SHOSet::test_derivatives" << std::endl;

  const SimulationCellT<T> simulation_cell;
  ParticleSetT<T> Ps(simulation_cell);

  int p = 0;
  PosType r, rtmp;
  for (int i = 0; i < n; ++i)
  {
    r[0] = c[0] + i * drg[0];
    for (int j = 0; j < n; ++j)
    {
      r[1] = c[1] + j * drg[1];
      for (int k = 0; k < n; ++k)
      {
        r[2] = c[2] + k * drg[2];

        evaluate_vgl(r, psi, dpsi, d2psi);

        for (int m = 0; m < nphi; ++m)
          d2psin[m] = -2 * odr2sum * psi[m];
        for (int d = 0; d < QMCTraits::DIM; ++d)
        {
          rtmp = r;
          rtmp[d] += dr[d];
          evaluate_v(rtmp, psitmp);
          for (int m = 0; m < nphi; ++m)
          {
            T phi       = psitmp[m];
            dpsin[m][d] = phi * o2dr[d];
            d2psin[m] += phi * odr2[d];
          }
          rtmp = r;
          rtmp[d] -= dr[d];
          evaluate_v(rtmp, psitmp);
          for (int m = 0; m < nphi; ++m)
          {
            T phi = psitmp[m];
            dpsin[m][d] -= phi * o2dr[d];
            d2psin[m] += phi * odr2[d];
          }
        }

        RealType dphi_diff  = 0.0;
        RealType d2phi_diff = 0.0;
        for (int m = 0; m < nphi; ++m)
          for (int d = 0; d < QMCTraits::DIM; ++d)
            dphi_diff = std::max<RealType>(dphi_diff, std::abs(dpsi[m][d] - dpsin[m][d]) / std::abs(dpsin[m][d]));
        for (int m = 0; m < nphi; ++m)
          d2phi_diff = std::max<RealType>(d2phi_diff, std::abs(d2psi[m] - d2psin[m]) / std::abs(d2psin[m]));
        app_log() << "  " << p << " " << dphi_diff << " " << d2phi_diff << std::endl;
        app_log() << "    derivatives" << std::endl;
        for (int m = 0; m < nphi; ++m)
        {
          std::string qn = "";
          for (int d = 0; d < QMCTraits::DIM; ++d)
            qn += int2string(state_info[m].quantum_number[d]) + " ";
          app_log() << "    " << qn;
          for (int d = 0; d < QMCTraits::DIM; ++d)
            app_log() << real(dpsi[m][d]) << " ";
          app_log() << std::endl;
          app_log() << "    " << qn;
          for (int d = 0; d < QMCTraits::DIM; ++d)
            app_log() << real(dpsin[m][d]) << " ";
          app_log() << std::endl;
        }
        app_log() << "    laplacians" << std::endl;
        PosType x = r / length;
        for (int m = 0; m < nphi; ++m)
        {
          std::string qn = "";
          for (int d = 0; d < QMCTraits::DIM; ++d)
            qn += int2string(state_info[m].quantum_number[d]) + " ";
          app_log() << "    " << qn << real(d2psi[m] / psi[m]) << std::endl;
          app_log() << "    " << qn << real(d2psin[m] / psi[m]) << std::endl;
        }
        p++;
      }
    }
  }

  app_log() << "end SHOSet::test_derivatives" << std::endl;
}

template<typename T>
void SHOSetT<T>::test_overlap()
{
  app_log() << "SHOSet::test_overlap" << std::endl;

  // linear
  int d = 0;

  app_log() << "  length = " << length << std::endl;
  app_log() << "  prefactors" << std::endl;
  for (int n = 0; n < qn_max[d]; ++n)
    app_log() << "    " << n << " " << prefactors[n] << std::endl;

  app_log() << "  1d overlap" << std::endl;

  ValueVector vpsi;
  vpsi.resize(this->size());
  ValueVector psi(&vpsi[0], this->size());

  double xmax = 4.0;
  double dx   = .1;
  double dr   = length * dx;

  int nphi = qn_max[d];
  Array<double, 2> omat;
  omat.resize(nphi, nphi);
  for (int i = 0; i < nphi; ++i)
    for (int j = 0; j < nphi; ++j)
      omat(i, j) = 0.0;

  PosType xp = 0.0;
  for (double x = -xmax; x < xmax; x += dx)
  {
    xp[d] = x;
    evaluate_hermite(xp);
    evaluate_d0(xp, psi);

    for (int i = 0; i < nphi; ++i)
      for (int j = 0; j < nphi; ++j)
        omat(i, j) += bvalues(d, i) * bvalues(d, j) * dr;
  }

  for (int i = 0; i < nphi; ++i)
  {
    app_log() << std::endl;
    for (int j = 0; j < nphi; ++j)
      app_log() << omat(i, j) << " ";
  }
  app_log() << std::endl;

  // volumetric
  app_log() << "  3d overlap" << std::endl;
  double dV = dr * dr * dr;
  nphi      = this->size();
  omat.resize(nphi, nphi);
  for (int i = 0; i < nphi; ++i)
    for (int j = 0; j < nphi; ++j)
      omat(i, j) = 0.0;
  for (double x = -xmax; x < xmax; x += dx)
    for (double y = -xmax; y < xmax; y += dx)
      for (double z = -xmax; z < xmax; z += dx)
      {
        xp[0] = x;
        xp[1] = y;
        xp[2] = z;
        evaluate_hermite(xp);
        evaluate_d0(xp, psi);

        for (int i = 0; i < nphi; ++i)
          for (int j = 0; j < nphi; ++j)
            omat(i, j) += std::abs(psi[i] * psi[j]) * dV;
      }
  for (int i = 0; i < nphi; ++i)
  {
    app_log() << std::endl;
    for (int j = 0; j < nphi; ++j)
      app_log() << omat(i, j) << " ";
  }
  app_log() << std::endl;

  app_log() << "end SHOSet::test_overlap" << std::endl;
}

template<typename T>
void SHOSetT<T>::evaluateThirdDeriv(const ParticleSetT<T>& P, int first, int last, GGGMatrix& grad_grad_grad_logdet)
{
  not_implemented("evaluateThirdDeriv(P,first,last,dddlogdet)");
}

template<typename T>
void SHOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P,
                                      int first,
                                      int last,
                                      ValueMatrix& logdet,
                                      GradMatrix& dlogdet,
                                      HessMatrix& grad_grad_logdet)
{
  not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet)");
}

template<typename T>
void SHOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P,
                                      int first,
                                      int last,
                                      ValueMatrix& logdet,
                                      GradMatrix& dlogdet,
                                      HessMatrix& grad_grad_logdet,
                                      GGGMatrix& grad_grad_grad_logdet)
{
  not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
}

template<typename T>
void SHOSetT<T>::evaluateGradSource(const ParticleSetT<T>& P,
                                    int first,
                                    int last,
                                    const ParticleSetT<T>& source,
                                    int iat_src,
                                    GradMatrix& gradphi)
{
  not_implemented("evaluateGradSource(P,first,last,source,iat,dphi)");
}

template<typename T>
void SHOSetT<T>::evaluateGradSource(const ParticleSetT<T>& P,
                                    int first,
                                    int last,
                                    const ParticleSetT<T>& source,
                                    int iat_src,
                                    GradMatrix& grad_phi,
                                    HessMatrix& grad_grad_phi,
                                    GradMatrix& grad_lapl_phi)
{
  not_implemented("evaluateGradSource(P,first,last,source,iat,dphi,ddphi,dd2phi)");
}

// Class concrete types from ValueType
#ifndef QMC_COMPLEX
template class SHOSetT<double>;
template class SHOSetT<float>;
#else
template class SHOSetT<std::complex<double>>;
template class SHOSetT<std::complex<float>>;
#endif

} // namespace qmcplusplus
