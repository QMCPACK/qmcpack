//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LATTICE_ANALYZER_H
#define QMCPLUSPLUS_LATTICE_ANALYZER_H
#include "OhmmsPETE/TinyVector.h"
namespace qmcplusplus
{
/** enumeration for DTD_BConds specialization
 *
 * G = general cell with image-cell checks
 * S = special treatment of a general cell with Wigner-cell radius == Simulation cell
 * O = orthogonal cell
 * X = exhaustive search (reference implementation)
 */
enum
{
  PPPG = SUPERCELL_BULK,
  PPPS = SUPERCELL_BULK + 1,
  PPPO = SUPERCELL_BULK + 2,
  PPPX = SUPERCELL_BULK + 3,
  PPNG = SUPERCELL_SLAB,
  PPNS = SUPERCELL_SLAB + 1,
  PPNO = SUPERCELL_SLAB + 2,
  PPNX = SUPERCELL_SLAB + 3
};


///generic class to analyze a Lattice
template<typename T, unsigned D>
struct LatticeAnalyzer
{};

/** specialization for  3D lattice
*/
template<typename T>
struct LatticeAnalyzer<T, 3>
{
  using SingleParticlePos = TinyVector<T, 3>;
  using Tensor_t          = Tensor<T, 3>;
  ///SuperCell type
  int mySC;

  inline int operator()(const TinyVector<int, 3>& box) { return mySC = box[0] + 2 * (box[1] + box[2] * 2); }

  inline bool isDiagonalOnly(const Tensor_t& R) const
  {
    T offdiag = std::abs(R(0, 1)) + std::abs(R(0, 2)) + std::abs(R(1, 0)) + std::abs(R(1, 2)) + std::abs(R(2, 0)) +
        std::abs(R(2, 1));
    return (offdiag < std::numeric_limits<T>::epsilon());
  }

  inline SingleParticlePos calcSolidAngles(const TinyVector<SingleParticlePos, 3>& Rv,
                                           const SingleParticlePos& OneOverLength)
  {
    const T rad_to_deg = 180.0 / M_PI;
    return SingleParticlePos(rad_to_deg * std::acos(dot(Rv[0], Rv[1]) * OneOverLength[0] * OneOverLength[1]),
                             rad_to_deg * std::acos(dot(Rv[1], Rv[2]) * OneOverLength[1] * OneOverLength[2]),
                             rad_to_deg * std::acos(dot(Rv[2], Rv[0]) * OneOverLength[2] * OneOverLength[0]));
  }

  inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos, 3>& a)
  {
    T rMin = 0.5 * std::numeric_limits<T>::max();
    if (mySC == SUPERCELL_BULK || mySC == SUPERCELL_BULK + SOA_OFFSET) //bulk type
    {
      for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
          for (int k = -1; k <= 1; k++)
            if (i || j || k)
            {
              SingleParticlePos L = (static_cast<T>(i) * a[0] + static_cast<T>(j) * a[1] + static_cast<T>(k) * a[2]);
              rMin                = std::min(rMin, dot(L, L));
            }
    }
    else if (mySC == SUPERCELL_SLAB || mySC == SUPERCELL_SLAB + SOA_OFFSET) //slab type
    {
      for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
          if (i || j)
          {
            SingleParticlePos L = (static_cast<T>(i) * a[0] + static_cast<T>(j) * a[1]);
            rMin                = std::min(rMin, dot(L, L));
          }
    }
    else if (mySC == SUPERCELL_WIRE || mySC == SUPERCELL_WIRE + SOA_OFFSET) //wire
    {
      rMin = dot(a[0], a[0]);
    }
    return 0.5 * std::sqrt(rMin);
  }

  inline T calcSimulationCellRadius(TinyVector<SingleParticlePos, 3>& a)
  {
    T scr = 0.5 * std::numeric_limits<T>::max();
    //if(mySC == SUPERCELL_BULK)
    //{
    for (int i = 0; i < 3; ++i)
    {
      SingleParticlePos A   = a[i];
      SingleParticlePos B   = a[(i + 1) % 3];
      SingleParticlePos C   = a[(i + 2) % 3];
      SingleParticlePos BxC = cross(B, C);
      T dist                = 0.5 * std::abs(dot(A, BxC)) / std::sqrt(dot(BxC, BxC));
      scr                   = std::min(scr, dist);
    }
    //}
    //else if(mySC == SUPERCELL_SLAB)
    //{
    //  T a0mag  = std::sqrt(dot(a[0],a[0]));
    //  T a1mag  = std::sqrt(dot(a[1],a[1]));
    //  scr=0.5*std::min(a0mag,a1mag);
    //  //T dist = 0.5*std::abs(dot(A,BxC))/std::sqrt(dot(BxC,BxC));
    //  //T theta1 = dot (a[0], a[1])/(a0mag*a1mag);
    //  //T theta2 = M_PI - theta1;
    //  //T theta  = std::min (theta1, theta2);
    //  //T dist   = std::min (a0mag, a1mag);
    //  //scr=0.5*std::sin(theta)*dist;
    //  std::cout << " calcSimulationCellRadius for Slab" << std::endl;
    //}
    //else if(mySC == SUPERCELL_WIRE)
    //{
    //  scr=0.5*std::sqrt(dot(a[0],a[0]));
    //}
    return scr;
  }
};

/** specialization for  2D lattice
*/
template<typename T>
struct LatticeAnalyzer<T, 2>
{
  using SingleParticlePos = TinyVector<T, 2>;
  using Tensor_t          = Tensor<T, 2>;

  /** return supercell enum
   * @param[in] box[2] if box[i]==1, PBC
   * @return SUPERCELL_OPEN or SUPERCELL_BULK
   */
  inline int operator()(const TinyVector<int, 2>& box)
  {
    return (box[0] + 2 * box[1]) ? SUPERCELL_BULK : SUPERCELL_OPEN;
  }

  inline bool isDiagonalOnly(const Tensor<T, 2>& R) const
  {
    T offdiag = std::abs(R(0, 1)) + std::abs(R(1, 0));
    return (offdiag < std::numeric_limits<T>::epsilon());
  }

  inline SingleParticlePos calcSolidAngles(const TinyVector<SingleParticlePos, 2>& Rv,
                                           const SingleParticlePos& OneOverLength)
  {
    const T rad_to_deg = 180.0 / M_PI;
    return SingleParticlePos(rad_to_deg * std::acos(dot(Rv[0], Rv[1]) * OneOverLength[0] * OneOverLength[1]), 0.0);
  }

  inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos, 2>& a)
  {
    T rMin;
    T dotP               = dot(a[0], a[1]);
    SingleParticlePos L0 = a[0] - dotP * a[1];
    SingleParticlePos L1 = a[1] - dotP * a[0];
    rMin                 = 0.5 * std::min(std::sqrt(dot(L0, L0)), std::sqrt(dot(L1, L1)));
    return rMin;
  }

  inline T calcSimulationCellRadius(TinyVector<SingleParticlePos, 2>& a)
  {
    T a0mag  = std::sqrt(dot(a[0], a[0]));
    T a1mag  = std::sqrt(dot(a[1], a[1]));
    T theta1 = std::acos(dot(a[0], a[1]) / (a0mag * a1mag));
    T theta2 = M_PI - theta1;
    T theta  = std::min(theta1, theta2);
    T dist   = std::min(a0mag, a1mag);
    return 0.5 * std::sin(theta) * dist;
    // return calcWignerSeitzRadius(a);
  }
};

/** specialization for 1D lattice
*/
template<typename T>
struct LatticeAnalyzer<T, 1>
{
  using SingleParticlePos = TinyVector<T, 1>;
  inline bool isDiagonalOnly(const Tensor<T, 1>& R) const { return true; }

  inline int operator()(const TinyVector<int, 1>& box) { return (box[0]) ? SUPERCELL_BULK : SUPERCELL_OPEN; }

  inline T calcWignerSeitzRadius(TinyVector<SingleParticlePos, 1>& a) { return a[0] * 0.5; }
};

template<typename T>
inline bool found_shorter_base(TinyVector<TinyVector<T, 3>, 3>& rb)
{
  const T eps = 10.0 * std::numeric_limits<T>::epsilon();
  int imax    = 0;
  T r2max     = dot(rb[0], rb[0]);
  for (int i = 1; i < 3; i++)
  {
    T r2 = dot(rb[i], rb[i]);
    if ((r2 - r2max) > eps)
    {
      r2max = r2;
      imax  = i;
    }
  }

  T rmax = std::sqrt(r2max);
  T tol  = 2.0 * rmax * eps; //Error propagation for x^2

  TinyVector<TinyVector<T, 3>, 4> rb_new;
  rb_new[0] = rb[0] + rb[1] - rb[2];
  rb_new[1] = rb[0] + rb[2] - rb[1];
  rb_new[2] = rb[1] + rb[2] - rb[0];
  rb_new[3] = rb[0] + rb[1] + rb[2];
  for (int i = 0; i < 4; ++i)
  {
    T r2 = dot(rb_new[i], rb_new[i]);
    if ((r2 - r2max) < -tol)
    {
      rb[imax] = rb_new[i];
      return true;
    }
  }
  return false;
}
template<typename T>
inline void find_reduced_basis(TinyVector<TinyVector<T, 3>, 3>& rb)
{
  int maxIter = 10000;

  for (int count = 0; count < maxIter; count++)
  {
    TinyVector<TinyVector<T, 3>, 3> saved(rb);
    bool changed = false;
    for (int i = 0; i < 3; ++i)
    {
      rb[i]   = 0.0;
      changed = found_shorter_base(rb);
      rb[i]   = saved[i];
      if (changed)
        break;
    }
    if (!changed && !found_shorter_base(rb))
      return;
  }

  throw std::runtime_error("Reduced basis not found in allowed number of iterations. "
                           "Check unit cell or contact a developer.");
}

} // namespace qmcplusplus
#endif
