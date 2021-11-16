//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HEGGRID_H
#define QMCPLUSPLUS_HEGGRID_H

#include "Lattice/CrystalLattice.h"
#include <map>
#include <optional>

namespace qmcplusplus
{
template<typename T, unsigned D>
struct kpdata
{
  TinyVector<T, D> k;
  T k2;
  int g;
};


template<typename T, unsigned D>
bool kpdata_comp(const kpdata<T, D>& left, const kpdata<T, D>& right)
{
  return left.k2 < right.k2;
}


template<class T, unsigned D>
struct HEGGrid
{};


//three-d specialization
template<class T>
struct HEGGrid<T, 3>
{
  typedef CrystalLattice<T, 3> PL_t;
  typedef typename PL_t::SingleParticlePos_t PosType;
  typedef typename PL_t::Scalar_t RealType;

  ///number of kpoints of a half sphere excluding gamma
  int NumKptsHalf;
  ///maxmim ksq
  T MaxKsq;
  PL_t& Lattice;
  std::map<int, std::vector<PosType>> rs;
  std::vector<PosType> kpt;
  std::vector<T> mk2;
  std::vector<int> deg;
  std::vector<int> n_within_shell;
  PosType twist;


  typedef kpdata<T, 3> kpdata_t;
  typedef std::vector<kpdata_t> kpoints_t;

  std::optional<kpoints_t> kpoints_grid;
  int nctmp;


  HEGGrid(PL_t& lat) : Lattice(lat), twist(0.0), nctmp(-1)
  {
    n_within_shell.resize(31);
    n_within_shell[0]  = 1;
    n_within_shell[1]  = 7;
    n_within_shell[2]  = 19;
    n_within_shell[3]  = 27;
    n_within_shell[4]  = 33;
    n_within_shell[5]  = 57;
    n_within_shell[6]  = 81;
    n_within_shell[7]  = 93;
    n_within_shell[8]  = 123;
    n_within_shell[9]  = 147;
    n_within_shell[10] = 171;
    n_within_shell[11] = 179;
    n_within_shell[12] = 203;
    n_within_shell[13] = 251;
    n_within_shell[14] = 257;
    n_within_shell[15] = 305;
    n_within_shell[16] = 341;
    n_within_shell[17] = 365;
    n_within_shell[18] = 389;
    n_within_shell[19] = 437;
    n_within_shell[20] = 461;
    n_within_shell[21] = 485;
    n_within_shell[22] = 515;
    n_within_shell[23] = 587;
    n_within_shell[24] = 619;
    n_within_shell[25] = 691;
    n_within_shell[26] = 739;
    n_within_shell[27] = 751;
    n_within_shell[28] = 799;
    n_within_shell[29] = 847;
    n_within_shell[30] = 895;
  }

  ~HEGGrid() = default;

  /** return the estimated number of grid in each direction */
  inline int getNC(int nup) const { return static_cast<int>(std::pow(static_cast<T>(nup), 1.0 / 3.0)) / 2 + 1; }

  /** return the estimated number of grid in each direction (upper bound) */
  inline int get_nc(int nstates) const
  {
    return static_cast<int>(std::pow(static_cast<T>(nstates), 1.0 / 3.0) * .7) + 1;
  }

  //return the number of k-points upto nsh-shell
  inline int getNumberOfKpoints(int nsh) const
  {
    if (nsh < n_within_shell.size())
      return n_within_shell[nsh];
    else
      return -1;
  }


  inline int getShellFromStates(int nst)
  {
    for (int i = 0; i < n_within_shell.size(); i++)
      if (n_within_shell[i] == nst)
        return i;
    return -1;
  }

  //return the shell index for nkpt k-points
  inline int getShellIndex(int nkpt) const
  {
    std::vector<int>::const_iterator loc = std::upper_bound(n_within_shell.begin(), n_within_shell.end(), nkpt);
    if (loc < n_within_shell.end())
      return loc - n_within_shell.begin() - 1;
    else
      return getNC(nkpt);
  }

  /** return the cell size  for the number of particles and rs
   * @param nptcl number of particles
   * @param rs_in rs
   */
  inline T getCellLength(int nptcl, T rs_in) const { return std::pow(4.0 * M_PI * nptcl / 3.0, 1.0 / 3.0) * rs_in; }

  void sortGrid(int nc)
  {
    int first_ix2, first_ix3;
    for (int ix1 = 0; ix1 <= nc; ix1++)
    {
      if (ix1 == 0)
        first_ix2 = 0;
      else
        first_ix2 = -nc;
      for (int ix2 = first_ix2; ix2 <= nc; ix2++)
      {
        if (ix1 == 0 && ix2 == 0)
          first_ix3 = 1;
        else
          first_ix3 = -nc;
        for (int ix3 = first_ix3; ix3 <= nc; ix3++)
        {
          int ih = ix1 * ix1 + ix2 * ix2 + ix3 * ix3;
          if (auto it = rs.find(ih); it == rs.end())
            rs[ih] = {PosType(ix1, ix2, ix3)};
          else
            it->second.push_back(PosType(ix1, ix2, ix3));
        }
      }
    }
  }

  void createGrid(int nc, int nkpts)
  {
    if (rs.empty())
      sortGrid(nc);
    NumKptsHalf = nkpts;
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    int ikpt = 0;
    MaxKsq    = 0.0;
    auto rs_it  = rs.begin();
    auto rs_end = rs.end();
    while (ikpt < nkpts && rs_it != rs_end)
    {
      auto ns_it  = rs_it->second.begin();
      auto ns_end = rs_it->second.end();
      T minus_ksq = -Lattice.ksq(*ns_it);
      while (ikpt < nkpts && ns_it != ns_end)
      {
        kpt[ikpt] = Lattice.k_cart(*ns_it);
        mk2[ikpt] = minus_ksq;
        ++ikpt;
        ++ns_it;
      }
      ++rs_it;
    }
    MaxKsq = Lattice.ksq(rs_it->second.front());
    app_log() << "List of kpoints (half-sphere) " << std::endl;
    for (int ik = 0; ik < kpt.size(); ik++)
    {
      app_log() << ik << " " << kpt[ik] << " " << mk2[ik] << std::endl;
    }
  }


  void clear_kpoints() { kpoints_grid.reset(); }


  void create_kpoints(int nc, const PosType& tw, T tol = 1e-6)
  {
    if (!kpoints_grid.has_value())
      kpoints_grid = kpoints_t();
    else if (nc <= nctmp)
      return;
    nctmp              = nc;
    kpoints_t& kpoints = *kpoints_grid;
    app_log() << "  resizing kpoint grid" << std::endl;
    app_log() << "  current size = " << kpoints.size() << std::endl;
    // make space for the kpoint grid
    int nkpoints = pow(2 * (nc + 1) + 1, 3);
    kpoints.resize(nkpoints);
    app_log() << "  cubic size = " << kpoints.size() << std::endl;
    typename kpoints_t::iterator kptmp, kp = kpoints.begin(), kp_end = kpoints.end();
    // make the kpoint grid
    T k2max = std::numeric_limits<RealType>::max();
    for (int i0 = -nc - 1; i0 <= nc + 1; ++i0)
      for (int i1 = -nc - 1; i1 <= nc + 1; ++i1)
        for (int i2 = -nc - 1; i2 <= nc + 1; ++i2)
        {
          PosType k(i0 + tw[0], i1 + tw[1], i2 + tw[2]);
          kp->k  = Lattice.k_cart(k);
          kp->k2 = Lattice.ksq(k);
          if (std::abs(i0) == (nc + 1) || std::abs(i1) == (nc + 1) || std::abs(i2) == (nc + 1))
            k2max = std::min(k2max, kp->k2);
          ++kp;
        }
    // sort kpoints by magnitude
    sort(kpoints.begin(), kpoints.end(), kpdata_comp<T, 3>);
    // eliminate kpoints outside of inscribing sphere
    int nkp = 0;
    kp      = kpoints.begin();
    while (kp != kp_end && kp->k2 < k2max + 1e-3)
    {
      nkp++;
      ++kp;
    }
    kpoints.resize(nkp);
    app_log() << "  new spherical size = " << kpoints.size() << std::endl;
    kp_end = kpoints.end();
    // count degeneracies
    kp = kpoints.begin();
    while (kp != kp_end)
    {
      T k2  = kp->k2;
      kptmp = kp;
      int g = 1;
      ++kptmp;
      // look ahead to count
      while (kptmp != kp_end && std::abs(kptmp->k2 - k2) < tol)
      {
        g++;
        ++kptmp;
      }
      kp->g = g;
      // run over degenerate states to assign
      for (int n = 0; n < g - 1; ++n)
        (++kp)->g = g;
      ++kp;
    }
    //app_log()<<"create_kpoints"<< std::endl;
    //app_log()<<"  nkpoints = "<<nkpoints<< std::endl;
    //app_log()<<"  kpoints"<< std::endl;
    //for(kp=kpoints.begin();kp!=kp_end;++kp)
    //  app_log()<<"    "<<kp->k2<<" "<<kp->g<<" "<<kp->k<< std::endl;
    //APP_ABORT("end create_kpoints");
  }


  void createGrid(int nc, int nkpts, const PosType& twistAngle)
  {
    twist = twistAngle;
    create_kpoints(nc, twistAngle);
    kpoints_t& kpoints = *kpoints_grid;
    if (nkpts > kpoints.size())
      APP_ABORT("HEGGrid::createGrid  requested more kpoints than created");
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    deg.resize(nkpts);
    for (int i = 0; i < nkpts; ++i)
    {
      const kpdata_t& kp = kpoints[i];
      kpt[i]             = kp.k;
      mk2[i]             = -kp.k2;
      deg[i]             = kp.g;
    }
    app_log() << "List of kpoints with twist = " << twistAngle << std::endl;
    for (int ik = 0; ik < kpt.size(); ik++)
      app_log() << ik << " " << kpt[ik] << " " << -mk2[ik] << std::endl;
  }


  void createGrid(const std::vector<int>& states, T tol = 1e-6) { createGrid(states, twist, tol); }


  void createGrid(const std::vector<int>& states, const PosType& twistAngle, T tol = 1e-6)
  {
    int smax = 0;
    for (int i = 0; i < states.size(); ++i)
      smax = std::max(smax, states[i]);
    smax++;
    create_kpoints(get_nc(smax), twistAngle, tol);
    kpoints_t& kpoints = *kpoints_grid;
    if (smax > kpoints.size())
      APP_ABORT("HEGGrid::createGrid(states)  requested more kpoints than created");
    int nkpts = states.size();
    kpt.resize(nkpts);
    mk2.resize(nkpts);
    deg.resize(nkpts);
    for (int i = 0; i < states.size(); ++i)
    {
      const kpdata_t& kp = kpoints[states[i]];
      kpt[i]             = kp.k;
      mk2[i]             = -kp.k2;
      deg[i]             = kp.g;
    }
  }
};
} // namespace qmcplusplus
#endif
