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


/** @file PWBasis.cpp
 * @brief Definition of member functions of Plane-wave basis set
 */
#include "PWBasis.h"

namespace qmcplusplus
{
int PWBasis::readbasis(hdf_archive& h5basisgroup,
                       const ParticleLayout& lat,
                       bool resizeContainer)
{
  ///make a local copy
  Lattice = lat;
  h5basisgroup.read(gvecs, "/electrons/kpoint_0/gvectors");
  NumPlaneWaves = std::max(gvecs.size(), kplusgvecs_cart.size());
  if (NumPlaneWaves == 0)
    throw std::runtime_error("  PWBasis::readbasis Basis is missing.");

  if (kplusgvecs_cart.empty())
  {
    kplusgvecs_cart.resize(NumPlaneWaves);
    for (int i = 0; i < NumPlaneWaves; i++)
      kplusgvecs_cart[i] = Lattice.k_cart(gvecs[i]);
  }
  //app_log() << "  Gx Gy Gz " << std::endl;
  //for(int i=0; i<kplusgvecs_cart.size(); i++)
  //{
  //  app_log() << kplusgvecs_cart[i] << std::endl;
  //}
  if (resizeContainer)
    reset();
  //std::copy(gvecs.begin(),gvecs.end(),std::ostream_iterator<GIndex_t>(std::cout,"\n"));
  return NumPlaneWaves;
}

void PWBasis::setTwistAngle(const PosType& tang)
{
  PosType dang   = twist - tang;
  bool sameTwist = dot(dang, dang) < std::numeric_limits<RealType>::epsilon();
  if (maxmaxg && sameTwist)
    return;
  twist = tang;
  reset();
}

void PWBasis::reset()
{
  rebuildBasis();
  //logC.resize(3,2*maxmaxg+1);
  Z.resize(NumPlaneWaves, 2 + DIM);
  Zv.resize(NumPlaneWaves);
  phi.resize(NumPlaneWaves);
}

void PWBasis::rebuildBasis()
{
  //Convert the twist angle to Cartesian coordinates.
  twist_cart = Lattice.k_cart(twist);
  inputmap.resize(NumPlaneWaves);
  app_log() << "  PWBasis::TwistAngle (unit) =" << twist << std::endl;
  app_log() << "  PWBasis::TwistAngle (cart) =" << twist_cart << std::endl;
  app_log() << "  PWBasis::rebuildBasis NumPlaneWaves (before) =" << NumPlaneWaves << std::endl;
  std::vector<GIndex_t> gvecCopy(gvecs);
  std::vector<PosType> gcartCopy(kplusgvecs_cart);
  gvecs.clear();
  kplusgvecs_cart.clear();
  minusModKplusG2.reserve(NumPlaneWaves);
  int ngIn = NumPlaneWaves;
  for (int ig = 0, newig = 0; ig < ngIn; ig++)
  {
    //PosType tempvec = Lattice.k_cart(gvecCopy[ig]+twist);
    PosType tempvec = gcartCopy[ig] + twist_cart;
    RealType mod2   = dot(tempvec, tempvec);

    gvecs.push_back(gvecCopy[ig]);
    kplusgvecs_cart.push_back(tempvec);
    minusModKplusG2.push_back(-mod2);
    //Remember which position in the HDF5 file this came from for coefficients.
    inputmap[ig] = newig++;
  }
#if defined(PWBASIS_USE_RECURSIVE)
  //Store the maximum number of translations of any reciprocal cell vector.
  for (int ig = 0; ig < NumPlaneWaves; ig++)
    for (int i = 0; i < OHMMS_DIM; i++)
      if (std::abs(gvecs[ig][i]) > maxg[i])
        maxg[i] = std::abs(gvecs[ig][i]);
  gvecs_shifted.resize(NumPlaneWaves);
  for (int ig = 0; ig < NumPlaneWaves; ig++)
    gvecs_shifted[ig] = gvecs[ig] + maxg;
  maxmaxg = std::max(maxg[0], std::max(maxg[1], maxg[2]));
  //changes the order???? ok
  C.resize(3, 2 * maxmaxg + 2);
#else
  maxmaxg = 1;
#endif
  app_log() << "                       NumPlaneWaves (after)  =" << NumPlaneWaves << std::endl;
}
} // namespace qmcplusplus
