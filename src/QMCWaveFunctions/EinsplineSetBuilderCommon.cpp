//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * Instantiates the static data
 * Implements member functions of EinsplineSetBuilder
 * - EinsplineSetBuilder
 * -
*/
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include <Message/UniformCommunicateError.h>
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
//std::map<H5OrbSet,SPOSet*,H5OrbSet>  EinsplineSetBuilder::SPOSetMap;
//std::map<TinyVector<int,4>,EinsplineSetBuilder::OrbType*,Int4less> EinsplineSetBuilder::OrbitalMap;
////std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet> EinsplineSetBuilder::ExtendedMap_z;
////std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet> EinsplineSetBuilder::ExtendedMap_d;

EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, const PSetMap& psets, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("spline", comm),
      ParticleSets(psets),
      TargetPtcl(p),
      XMLRoot(cur),
      H5FileID(-1),
      Format(QMCPACK),
      NumBands(0),
      NumElectrons(0),
      NumSpins(0),
      NumTwists(0),
      NumCoreStates(0),
      MeshFactor(1.0),
      MeshSize(0, 0, 0),
      twist_num_(-1),
      TileFactor(1, 1, 1),
      NumMuffinTins(0),
      LastSpinSet(-1),
      NumOrbitalsRead(-1),
      makeRotations(false)
{
  ClassName = "EinsplineSetBuilder";

  MatchingTol = 10 * std::numeric_limits<float>::epsilon();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      TileMatrix(i, j) = 0;

  //invalidate states by the basis class
  states.clear();
  states.resize(p.groups());

  //create vectors with nullptr
  FullBands.resize(p.groups());
}

template<typename T>
inline TinyVector<T, 3> IntPart(const TinyVector<T, 3>& twist)
{
  return TinyVector<T, 3>(round(twist[0] - 1.0e-6), round(twist[1] - 1.0e-6), round(twist[2] - 1.0e-6));
}

template<typename T>
inline TinyVector<T, 3> FracPart(const TinyVector<T, 3>& twist)
{
  return twist - IntPart(twist);
}


EinsplineSetBuilder::~EinsplineSetBuilder()
{
  DEBUG_MEMORY("EinsplineSetBuilder::~EinsplineSetBuilder");
  if (H5FileID >= 0)
    H5Fclose(H5FileID);
}


bool EinsplineSetBuilder::CheckLattice()
{
  double diff = 0.0;
  for (int i = 0; i < OHMMS_DIM; i++)
    for (int j = 0; j < OHMMS_DIM; j++)
    {
      double max_abs =
          std::max(std::abs(SuperLattice(i, j)), static_cast<double>(std::abs(TargetPtcl.getLattice().R(i, j))));
      if (max_abs > MatchingTol)
        diff = std::max(diff, std::abs(SuperLattice(i, j) - TargetPtcl.getLattice().R(i, j)) / max_abs);
    }

  if (diff > MatchingTol)
  {
    std::ostringstream o;
    o.setf(std::ios::scientific, std::ios::floatfield);
    o.precision(6);
    o << "EinsplineSetBuilder::ReadOrbitalInfo_ESHDF \n"
      << "Mismatched supercell lattices.\n";
    o << " Lattice in ESHDF5 " << std::endl;
    o << SuperLattice << std::endl;
    o << " Lattice in xml" << std::endl;
    o << TargetPtcl.getLattice().R << std::endl;
    o << " Difference " << std::endl;
    o << SuperLattice - TargetPtcl.getLattice().R << std::endl;
    o << " Max relative error = " << diff << std::endl;
    o << " Tolerance      = " << MatchingTol << std::endl;
    app_error() << o.str();
    return false;
  }
  return true;
}

void EinsplineSetBuilder::BroadcastOrbitalInfo()
{
  if (myComm->size() == 1)
    return;
  int numIons           = IonTypes.size();
  int numAtomicOrbitals = AtomicOrbitals.size();
  int numDensityGvecs   = TargetPtcl.DensityReducedGvecs.size();
  PooledData<double> abuffer;
  PooledData<int> aibuffer;
  aibuffer.add(Version.begin(), Version.end()); //myComm->bcast(Version);
  aibuffer.add(Format);
  abuffer.add(Lattice.begin(), Lattice.end());           //myComm->bcast(Lattice);
  abuffer.add(RecipLattice.begin(), RecipLattice.end()); //myComm->bcast(RecipLattice);
  abuffer.add(SuperLattice.begin(), SuperLattice.end()); //myComm->bcast(SuperLattice);
  abuffer.add(LatticeInv.begin(), LatticeInv.end());     //myComm->bcast(LatticeInv);
  aibuffer.add(NumBands);                                //myComm->bcast(NumBands);
  aibuffer.add(NumElectrons);                            //myComm->bcast(NumElectrons);
  aibuffer.add(NumSpins);                                //myComm->bcast(NumSpins);
  aibuffer.add(NumTwists);                               //myComm->bcast(NumTwists);
  aibuffer.add(numIons);                                 //myComm->bcast(numIons);
  aibuffer.add(NumMuffinTins);
  aibuffer.add(numAtomicOrbitals);
  aibuffer.add(numDensityGvecs);
  aibuffer.add(HaveOrbDerivs);
  myComm->bcast(abuffer);
  myComm->bcast(aibuffer);
  if (myComm->rank())
  {
    abuffer.rewind();
    aibuffer.rewind();
    aibuffer.get(Version.begin(), Version.end());
    aibuffer.get(Format);
    abuffer.get(Lattice.begin(), Lattice.end());
    abuffer.get(RecipLattice.begin(), RecipLattice.end());
    abuffer.get(SuperLattice.begin(), SuperLattice.end());
    abuffer.get(LatticeInv.begin(), LatticeInv.end());
    aibuffer.get(NumBands);
    aibuffer.get(NumElectrons);
    aibuffer.get(NumSpins);
    aibuffer.get(NumTwists);
    aibuffer.get(numIons);
    aibuffer.get(NumMuffinTins);
    aibuffer.get(numAtomicOrbitals);
    aibuffer.get(numDensityGvecs);
    aibuffer.get(HaveOrbDerivs);
    MT_APW_radii.resize(NumMuffinTins);
    MT_APW_lmax.resize(NumMuffinTins);
    MT_APW_rgrids.resize(NumMuffinTins);
    MT_APW_num_radial_points.resize(NumMuffinTins);
    MT_centers.resize(NumMuffinTins);
    TargetPtcl.DensityReducedGvecs.resize(numDensityGvecs);
    TargetPtcl.Density_G.resize(numDensityGvecs);
    AtomicOrbitals.resize(numAtomicOrbitals);
  }
  std::vector<int> rgrids_sizes(NumMuffinTins);
  for (int tin = 0; tin < NumMuffinTins; tin++)
    rgrids_sizes[tin] = MT_APW_rgrids[tin].size();
  myComm->bcast(rgrids_sizes);
  if (myComm->rank())
    for (int tin = 0; tin < NumMuffinTins; tin++)
      MT_APW_rgrids[tin].resize(rgrids_sizes[tin]);
  if (IonTypes.size() != numIons)
  {
    IonTypes.resize(numIons);
    IonPos.resize(numIons);
  }
  //new buffer
  PooledData<double> bbuffer;
  PooledData<int> bibuffer;
  for (int i = 0; i < numIons; ++i)
    bibuffer.add(IonTypes[i]);
  //myComm->bcast(IonTypes);
  bbuffer.add(&IonPos[0][0], &IonPos[0][0] + OHMMS_DIM * numIons);
  //myComm->bcast(IonPos);
  if (TwistAngles.size() != NumTwists)
    TwistAngles.resize(NumTwists);
  bbuffer.add(&TwistAngles[0][0], &TwistAngles[0][0] + OHMMS_DIM * NumTwists);
  //myComm->bcast(TwistAngles);
  if (TwistSymmetry.size() != NumTwists)
    TwistSymmetry.resize(NumTwists);
  bibuffer.add(&TwistSymmetry[0], &TwistSymmetry[0] + NumTwists);
  if (TwistWeight.size() != NumTwists)
    TwistWeight.resize(NumTwists);
  bibuffer.add(&TwistWeight[0], &TwistWeight[0] + NumTwists);
  bbuffer.add(MT_APW_radii.begin(), MT_APW_radii.end());
  bibuffer.add(MT_APW_lmax.begin(), MT_APW_lmax.end());
  bibuffer.add(MT_APW_num_radial_points.begin(), MT_APW_num_radial_points.end());
  bbuffer.add(&(MT_centers[0][0]), &(MT_centers[0][0]) + OHMMS_DIM * NumMuffinTins);
  for (int i = 0; i < NumMuffinTins; i++)
    bbuffer.add(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
  bibuffer.add(&(TargetPtcl.DensityReducedGvecs[0][0]),
               &(TargetPtcl.DensityReducedGvecs[0][0]) + numDensityGvecs * OHMMS_DIM);
  bbuffer.add(&(TargetPtcl.Density_G[0]), &(TargetPtcl.Density_G[0]) + numDensityGvecs);
  for (int iat = 0; iat < numAtomicOrbitals; iat++)
  {
    AtomicOrbital<std::complex<double>>& orb = AtomicOrbitals[iat];
    bibuffer.add(orb.SplinePoints);
    bibuffer.add(orb.PolyOrder);
    bibuffer.add(orb.lMax);
    bibuffer.add(orb.Numlm);
    bbuffer.add(&orb.Pos[0], &orb.Pos[0] + OHMMS_DIM);
    bbuffer.add(orb.CutoffRadius);
    bbuffer.add(orb.SplineRadius);
    bbuffer.add(orb.PolyRadius);
  }
  myComm->bcast(bbuffer);
  myComm->bcast(bibuffer);
  if (myComm->rank())
  {
    bbuffer.rewind();
    bibuffer.rewind();
    for (int i = 0; i < numIons; ++i)
      bibuffer.get(IonTypes[i]);
    bbuffer.get(&IonPos[0][0], &IonPos[0][0] + OHMMS_DIM * numIons);
    bbuffer.get(&TwistAngles[0][0], &TwistAngles[0][0] + OHMMS_DIM * NumTwists);
    bibuffer.get(&TwistSymmetry[0], &TwistSymmetry[0] + NumTwists);
    bibuffer.get(&TwistWeight[0], &TwistWeight[0] + NumTwists);
    bbuffer.get(MT_APW_radii.begin(), MT_APW_radii.end());
    bibuffer.get(MT_APW_lmax.begin(), MT_APW_lmax.end());
    bibuffer.get(MT_APW_num_radial_points.begin(), MT_APW_num_radial_points.end());
    bbuffer.get(&(MT_centers[0][0]), &(MT_centers[0][0]) + OHMMS_DIM * NumMuffinTins);
    for (int i = 0; i < NumMuffinTins; i++)
      bbuffer.get(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
    bibuffer.get(&(TargetPtcl.DensityReducedGvecs[0][0]),
                 &(TargetPtcl.DensityReducedGvecs[0][0]) + numDensityGvecs * OHMMS_DIM);
    bbuffer.get(&(TargetPtcl.Density_G[0]), &(TargetPtcl.Density_G[0]) + numDensityGvecs);
    for (int iat = 0; iat < numAtomicOrbitals; iat++)
    {
      AtomicOrbital<std::complex<double>>& orb = AtomicOrbitals[iat];
      bibuffer.get(orb.SplinePoints);
      bibuffer.get(orb.PolyOrder);
      bibuffer.get(orb.lMax);
      bibuffer.get(orb.Numlm);
      bbuffer.get(&orb.Pos[0], &orb.Pos[0] + OHMMS_DIM);
      bbuffer.get(orb.CutoffRadius);
      bbuffer.get(orb.SplineRadius);
      bbuffer.get(orb.PolyRadius);
    }
  }
  //buffer to bcast hybrid representation atomic orbital info
  PooledData<double> cbuffer;
  PooledData<int> cibuffer;
  myComm->bcast(cbuffer);
  myComm->bcast(cibuffer);
  AtomicCentersInfo.resize(numIons);
  Super2Prim.resize(SourcePtcl->R.size());
  cbuffer.add(AtomicCentersInfo.inner_cutoff.begin(), AtomicCentersInfo.inner_cutoff.end());
  cbuffer.add(AtomicCentersInfo.non_overlapping_radius.begin(), AtomicCentersInfo.non_overlapping_radius.end());
  cbuffer.add(AtomicCentersInfo.cutoff.begin(), AtomicCentersInfo.cutoff.end());
  cbuffer.add(AtomicCentersInfo.spline_radius.begin(), AtomicCentersInfo.spline_radius.end());
  cibuffer.add(Super2Prim.begin(), Super2Prim.end());
  cibuffer.add(AtomicCentersInfo.lmax.begin(), AtomicCentersInfo.lmax.end());
  cibuffer.add(AtomicCentersInfo.GroupID.begin(), AtomicCentersInfo.GroupID.end());
  cibuffer.add(AtomicCentersInfo.spline_npoints.begin(), AtomicCentersInfo.spline_npoints.end());
  myComm->bcast(cbuffer);
  myComm->bcast(cibuffer);
  if (myComm->rank())
  {
    cbuffer.rewind();
    cibuffer.rewind();
    cbuffer.get(AtomicCentersInfo.inner_cutoff.begin(), AtomicCentersInfo.inner_cutoff.end());
    cbuffer.get(AtomicCentersInfo.non_overlapping_radius.begin(), AtomicCentersInfo.non_overlapping_radius.end());
    cbuffer.get(AtomicCentersInfo.cutoff.begin(), AtomicCentersInfo.cutoff.end());
    cbuffer.get(AtomicCentersInfo.spline_radius.begin(), AtomicCentersInfo.spline_radius.end());
    cibuffer.get(Super2Prim.begin(), Super2Prim.end());
    cibuffer.get(AtomicCentersInfo.lmax.begin(), AtomicCentersInfo.lmax.end());
    cibuffer.get(AtomicCentersInfo.GroupID.begin(), AtomicCentersInfo.GroupID.end());
    cibuffer.get(AtomicCentersInfo.spline_npoints.begin(), AtomicCentersInfo.spline_npoints.end());
    for (int i = 0; i < numIons; i++)
      AtomicCentersInfo.ion_pos[i] = IonPos[i];
  }
}

////////////////////////////////////////////////////////////////
//// Create the ion ParticleSet from the data in the HDF file //
////////////////////////////////////////////////////////////////
//void
//EinsplineSetBuilder::CreateIonParticleSet( std::string sourceName)
//{
//  //    ParticleSet &pTemp = *(new MCWalkerConfiguration);
//  ParticleSet &pTemp = *(new ParticleSet);
//  pTemp.setName (sourceName);
//  SpeciesSet& tspecies(pTemp.getSpeciesSet());
//  ParticleSets[sourceName] = &pTemp;
//}
//

void EinsplineSetBuilder::TileIons()
{
  //set the primitive lattice
  SourcePtcl->getPrimitiveLattice().set(Lattice);

  for (int j = 0; j < IonPos.size(); ++j)
    IonPos[j] = FracPart(SourcePtcl->getPrimitiveLattice().toUnit(IonPos[j]));

  IonPos.resize(SourcePtcl->getTotalNum());
  IonTypes.resize(SourcePtcl->getTotalNum());
  std::copy(SourcePtcl->R.begin(), SourcePtcl->R.end(), IonPos.begin());
  std::copy(SourcePtcl->GroupID.begin(), SourcePtcl->GroupID.end(), IonTypes.begin());

  //app_log() << "  Primitive Cell\n";
  //SourcePtcl->getPrimitiveLattice().print(app_log());
  //app_log() << "  Super Cell\n";
  //SourcePtcl->Lattice.print(app_log());

  //Don't need to do this, already one by ParticleSetPool.cpp
  //  Vector<TinyVector<double, OHMMS_DIM> > primPos   = IonPos;
  //  Vector<int>                            primTypes = IonTypes;
  //  int numCopies = std::abs(det(TileMatrix));
  //  IonTypes.resize(primPos.size()*numCopies);
  //  IonPos.resize  (primPos.size()*numCopies);
  //  int maxCopies = 10;
  //  using Vec3 = TinyVector<double,3>;
  //  int index=0;
  //  for (int i0=-maxCopies; i0<=maxCopies; i0++)
  //    for (int i1=-maxCopies; i1<=maxCopies; i1++)
  //      for (int i2=-maxCopies; i2<=maxCopies; i2++)
  //        for (int iat=0; iat < primPos.size(); iat++)
  //        {
  //          Vec3 r     = primPos[iat];
  //          Vec3 uPrim = PrimCell.toUnit(r);
  //          for (int i=0; i<3; i++)
  //            uPrim[i] -= std::floor(uPrim[i]);
  //          r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0) +
  //              (double)i1*PrimCell.a(1) + (double)i2*PrimCell.a(2);
  //          Vec3 uSuper = SuperCell.toUnit(r);
  //          if ((uSuper[0] >= -1.0e-4) && (uSuper[0] < 0.9999) &&
  //              (uSuper[1] >= -1.0e-4) && (uSuper[1] < 0.9999) &&
  //              (uSuper[2] >= -1.0e-4) && (uSuper[2] < 0.9999))
  //          {
  //            IonPos[index]= r;
  //            IonTypes[index]= primTypes[iat];
  //            index++;
  //          }
  //        }
  //  if (index != primPos.size()*numCopies)
  //  {
  //    app_error() << "The number of tiled ions, " << IonPos.size()
  //                << ", does not match the expected number of "
  //                << primPos.size()*numCopies << " or the index "<< index <<".  Aborting.\n";
  //    APP_ABORT("EinsplineSetBuilder::TileIons()");
  //  }
  //  if (myComm->rank() == 0)
  //  {
  //    char buf[1000];
  //    snprintf (buf, 1000, "Supercell reduced ion positions = \n");
  //    app_log() << buf;
  //    app_log().flush();
  //    for (int i=0; i<IonPos.size(); i++)
  //    {
  //      PosType u = SuperCell.toUnit(IonPos[i]);
  //      char buf2[1000];
  //      snprintf (buf2, 1000, "   %14.10f %14.10f %14.10f\n",
  //               u[0], u[1], u[2]);
  //      app_log() << buf2;
  //      app_log().flush();
  //      //		 IonPos[i][0], IonPos[i][1], IonPos[i][2]);
  //    }
  //  }
}


bool EinsplineSetBuilder::TwistPair(PosType a, PosType b)
{
  bool pair = true;
  for (int n = 0; n < OHMMS_DIM; n++)
  {
    double d = a[n] + b[n];
    if (std::abs(d - round(d)) > MatchingTol)
      pair = false;
  }
  return pair;
}

void EinsplineSetBuilder::AnalyzeTwists2(const int twist_num_inp, const TinyVector<double, OHMMS_DIM>& twist_inp)
{
  Tensor<double, 3> S;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      S(i, j) = (double)TileMatrix(i, j);
  std::vector<PosType> superFracs;
  std::vector<std::vector<int>> superSets;
  { // build super twists
    // This holds to which supercell kpoint each primitive k-point belongs
    std::vector<int> superIndex;
    const int numPrimTwists = TwistAngles.size();
    for (int ki = 0; ki < numPrimTwists; ki++)
    {
      PosType primTwist  = TwistAngles[ki];
      PosType superTwist = dot(S, primTwist);
      PosType kp         = PrimCell.k_cart(primTwist);
      PosType ks         = SuperCell.k_cart(superTwist);
      if (dot(ks - kp, ks - kp) > 1.0e-6)
      {
        app_error() << "Primitive and super k-points do not agree.  Error in coding.\n";
        APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
      }
      PosType frac = FracPart(superTwist);
      bool found   = false;
      for (int j = 0; j < superFracs.size(); j++)
      {
        PosType diff = frac - superFracs[j];
        if (dot(diff, diff) < 1.0e-6)
        {
          found = true;
          superIndex.push_back(j);
        }
      }
      if (!found)
      {
        superIndex.push_back(superFracs.size());
        superFracs.push_back(frac);
      }
    }
    const int numSuperTwists = superFracs.size();
    app_log() << "Found " << numSuperTwists << " distinct supercell twists.\n";
    // For each supercell twist, create a list of primitive twists which
    // belong to it.
    superSets.resize(numSuperTwists);
    for (int ki = 0; ki < numPrimTwists; ki++)
      superSets[superIndex[ki]].push_back(ki);
    app_log() << "number of things" << std::endl;
    app_log() << TwistSymmetry.size() << std::endl;
    app_log() << TwistWeight.size() << std::endl;
    //     for (int ki=0; ki<TwistSymmetry.size(); ki++)
    //       fprintf (stderr, "%d %d %d\n",ki,TwistSymmetry[ki],TwistWeight[ki]);
    if (myComm->rank() == 0)
    {
      int n_tot_irred(0);
      for (int si = 0; si < numSuperTwists; si++)
      {
        char buf[1000];
        snprintf(buf, 1000, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n", si, superFracs[si][0], superFracs[si][1],
                 superFracs[si][2]);
        app_log() << buf;
        app_log().flush();
      }
    }
  }
  const int numSuperTwists = superFracs.size();

  { // determine twist_num_
    std::function find_twist = [&](const TinyVector<double, OHMMS_DIM>& twist) {
      int twist_num  = -1;
      PosType gtFrac = FracPart(twist);
      float eps      = 1e-5;
      for (int si = 0; si < numSuperTwists; si++)
      {
        PosType locDiff = gtFrac - superFracs[si];
        if (dot(locDiff, locDiff) < eps)
          twist_num = si;
      }

      if (twist_num < 0)
      {
        char buf[1000];
        snprintf(buf, 1000,
                 "AnalyzeTwists2. Input twist [ %9.5f %9.5f %9.5f] not found in the list of super twists above.\n",
                 twist[0], twist[1], twist[2]);
        throw UniformCommunicateError(buf);
      }
      return twist_num;
    };

    if (twist_inp[0] > TWIST_NO_INPUT || twist_inp[1] > TWIST_NO_INPUT || twist_inp[2] > TWIST_NO_INPUT)
    {
      if (twist_num_inp != TWISTNUM_NO_INPUT)
        app_warning() << "twist attribute exists. twistnum attribute ignored. "
                         "To prevent this message, remove twistnum from input."
                      << std::endl;

      twist_num_ = find_twist(twist_inp);
    }
    else if (twist_num_inp != TWISTNUM_NO_INPUT)
    {
      app_warning() << "twist attribute does't exist but twistnum attribute was found. "
                    << "This is potentially ambiguous. Specifying twist attribute is preferred." << std::endl;
      if (twist_num_inp < 0 || twist_num_inp >= numSuperTwists)
      {
        std::ostringstream msg;
        msg << "AnalyzeTwists2. twistnum input value " << twist_num_inp << " is outside the acceptable range [0, "
            << numSuperTwists << ")." << std::endl;
        throw UniformCommunicateError(msg.str());
      }
      twist_num_ = twist_num_inp;
    }
    else
    {
      app_log() << "twist attribte does't exist. Set Gamma point." << std::endl;
      twist_num_ = find_twist({0, 0, 0});
    }

    assert(twist_num_ >= 0 && twist_num_ < numSuperTwists);

    char buf[1000];
    snprintf(buf, 1000, "  Using supercell twist %d:  [ %9.5f %9.5f %9.5f]", twist_num_, superFracs[twist_num_][0],
             superFracs[twist_num_][1], superFracs[twist_num_][2]);
    app_log() << buf << std::endl;
  }

  TargetPtcl.setTwist(superFracs[twist_num_]);
#ifndef QMC_COMPLEX
  // Check to see if supercell twist is okay to use with real wave
  // functions
  for (int dim = 0; dim < OHMMS_DIM; dim++)
  {
    double t = 2.0 * superFracs[twist_num_][dim];
    if (std::abs(t - round(t)) > MatchingTol * 100)
    {
      app_error() << "Cannot use this super twist with real wavefunctions.\n"
                  << "Please recompile with QMC_COMPLEX=1.\n";
      APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
    }
  }
#endif
  // Now check to see that each supercell twist has the right twists
  // to tile the primitive cell orbitals.
  const int numTwistsNeeded = std::abs(det(TileMatrix));
  for (int si = 0; si < numSuperTwists; si++)
  {
    // First make sure we have enough points
    if (superSets[si].size() != numTwistsNeeded)
    {
      char buf[1000];
      snprintf(buf, 1000, "Super twist %d should own %d k-points, but owns %d.\n", si, numTwistsNeeded,
               static_cast<int>(superSets[si].size()));
      app_error() << buf;
      if (si == twist_num_)
      {
        APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
      }
      else
        continue;
    }
    // Now, make sure they are all distinct
    int N = superSets[si].size();
    for (int i = 0; i < N; i++)
    {
      PosType twistPrim_i  = TwistAngles[superSets[si][i]];
      PosType twistSuper_i = dot(S, twistPrim_i);
      PosType superInt_i   = IntPart(twistSuper_i);
      for (int j = i + 1; j < N; j++)
      {
        PosType twistPrim_j  = TwistAngles[superSets[si][j]];
        PosType twistSuper_j = dot(S, twistPrim_j);
        PosType superInt_j   = IntPart(twistSuper_j);
        if (dot(superInt_i - superInt_j, superInt_i - superInt_j) < 1.0e-6)
        {
          app_error() << "Identical k-points detected in super twist set " << si << std::endl;
          APP_ABORT_TRACE(__FILE__, __LINE__, "AnalyzeTwists2");
        }
      }
    }
  }
  app_log().flush();
  // Finally, record which k-points to include on this group of
  // processors, which have been assigned supercell twist twist_num_
  IncludeTwists.clear();
  for (int i = 0; i < superSets[twist_num_].size(); i++)
    IncludeTwists.push_back(superSets[twist_num_][i]);
  // Now, find out which twists are distinct
  DistinctTwists.clear();
#ifndef QMC_COMPLEX
  std::vector<int> copyTwists;
  for (int i = 0; i < IncludeTwists.size(); i++)
  {
    int ti          = IncludeTwists[i];
    PosType twist_i = TwistAngles[ti];
    bool distinct   = true;
    for (int j = i + 1; j < IncludeTwists.size(); j++)
    {
      int tj          = IncludeTwists[j];
      PosType twist_j = TwistAngles[tj];
      PosType sum     = twist_i + twist_j;
      PosType diff    = twist_i - twist_j;
      if (TwistPair(twist_i, twist_j))
        distinct = false;
    }
    if (distinct)
      DistinctTwists.push_back(ti);
    else
      copyTwists.push_back(ti);
  }
  // Now determine which distinct twists require two copies
  MakeTwoCopies.resize(DistinctTwists.size());
  for (int i = 0; i < DistinctTwists.size(); i++)
  {
    MakeTwoCopies[i] = false;
    int ti           = DistinctTwists[i];
    PosType twist_i  = TwistAngles[ti];
    for (int j = 0; j < copyTwists.size(); j++)
    {
      int tj          = copyTwists[j];
      PosType twist_j = TwistAngles[tj];
      if (TwistPair(twist_i, twist_j))
        MakeTwoCopies[i] = true;
    }
    if (myComm->rank() == 0)
    {
      char buf[1000];

      snprintf(buf, 1000, "Using %d copies of twist angle [%6.3f, %6.3f, %6.3f]\n", MakeTwoCopies[i] ? 2 : 1,
               twist_i[0], twist_i[1], twist_i[2]);
      app_log() << buf;
      app_log().flush();
    }
  }
  // Find out if we can make real orbitals
  use_real_splines_ = true;
  for (int i = 0; i < DistinctTwists.size(); i++)
  {
    int ti        = DistinctTwists[i];
    PosType twist = TwistAngles[ti];
    for (int j = 0; j < OHMMS_DIM; j++)
      if (std::abs(twist[j] - 0.0) > MatchingTol && std::abs(twist[j] - 0.5) > MatchingTol &&
          std::abs(twist[j] + 0.5) > MatchingTol)
        use_real_splines_ = false;
  }
  if (use_real_splines_ && (DistinctTwists.size() > 1))
  {
    app_log() << "***** Use of real orbitals is possible, but not currently implemented\n"
              << "      with more than one twist angle.\n";
    use_real_splines_ = false;
  }
  if (use_real_splines_)
    app_log() << "Using real splines.\n";
  else
    app_log() << "Using complex splines.\n";
#else
  DistinctTwists.resize(IncludeTwists.size());
  MakeTwoCopies.resize(IncludeTwists.size());
  for (int i = 0; i < IncludeTwists.size(); i++)
  {
    DistinctTwists[i] = IncludeTwists[i];
    MakeTwoCopies[i]  = false;
  }
  use_real_splines_ = false;
#endif
}


void EinsplineSetBuilder::OccupyBands(int spin, int sortBands, int numOrbs, bool skipChecks)
{
  if (myComm->rank() != 0)
    return;
  if (spin >= NumSpins && !skipChecks)
  {
    app_error() << "To developer: User is requesting for orbitals in an invalid spin group " << spin
                << ". Current h5 file only contains spin groups "
                << "[0.." << NumSpins - 1 << "]." << std::endl;
    app_error() << "To user: Orbital H5 file contains no spin down data and is appropriate only for spin unpolarized "
                   "calculations. "
                << "If this is your intent, please replace 'spindataset=1' with 'spindataset=0' in the input file."
                << std::endl;
    abort();
  }
  if (Format == ESHDF)
  {
    OccupyBands_ESHDF(spin, sortBands, numOrbs);
    return;
  }
  std::string eigenstatesGroup;
  if (Version[0] == 0 && Version[1] == 11)
    eigenstatesGroup = "/eigenstates_3";
  else if (Version[0] == 0 && Version[1] == 20)
    eigenstatesGroup = "/eigenstates";

  if (FullBands[spin]->size())
  {
    app_log() << "  FullBand[" << spin << "] exists. Reuse it. " << std::endl;
    return;
  }

  std::vector<BandInfo>& SortBands(*FullBands[spin]);

  SortBands.clear();
  for (int ti = 0; ti < DistinctTwists.size(); ti++)
  {
    int tindex = DistinctTwists[ti];
    // First, read valence states
    for (int bi = 0; bi < NumBands; bi++)
    {
      BandInfo band;
      band.IsCoreState   = false;
      band.TwistIndex    = tindex;
      band.BandIndex     = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      // Read eigenenergy from file
      std::ostringstream ePath, sPath;
      if ((Version[0] == 0 && Version[1] == 11) || NumTwists > 1)
      {
        ePath << eigenstatesGroup << "/twist_" << tindex << "/band_" << bi << "/eigenvalue";
        sPath << eigenstatesGroup << "/twist_" << tindex << "/band_" << bi << "/spin";
      }
      else if (NumBands > 1)
      {
        ePath << eigenstatesGroup << "/twist/band_" << bi << "/eigenvalue";
        sPath << eigenstatesGroup << "/twist/band_" << bi << "/spin";
      }
      else
      {
        ePath << eigenstatesGroup << "/twist/band/eigenvalue";
        sPath << eigenstatesGroup << "/twist/band/spin";
      }
      HDFAttribIO<double> h_energy(band.Energy);
      HDFAttribIO<int> h_spin(band.Spin);
      band.Energy = -1.01e100;
      h_energy.read(H5FileID, ePath.str().c_str());
      if (band.Energy > -1.0e100)
      {
        h_spin.read(H5FileID, sPath.str().c_str());
        if (band.Spin == spin)
          SortBands.push_back(band);
      }
    }
    // Now, read core states
    for (int cs = 0; cs < NumCoreStates; cs++)
    {
      BandInfo band;
      band.IsCoreState   = true;
      band.TwistIndex    = tindex;
      band.BandIndex     = cs;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      HDFAttribIO<double> h_energy(band.Energy);
      h_energy.read(H5FileID, (CoreStatePath(ti, cs) + "eigenvalue").c_str());
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
    }
  }
  int orbIndex        = 0;
  int numOrbs_counter = 0;
  NumValenceOrbs      = 0;
  NumCoreOrbs         = 0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    if (SortBands[orbIndex].IsCoreState)
      NumCoreOrbs++;
    else
      NumValenceOrbs++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
  app_log() << "There are " << NumCoreOrbs << " core states and " << NumValenceOrbs << " valence states.\n";
}

bool EinsplineSetBuilder::bcastSortBands(int spin, int n, bool root)
{
  std::vector<BandInfo>& SortBands(*FullBands[spin]);

  TinyVector<int, 4> nbands(int(SortBands.size()), n, NumValenceOrbs, NumCoreOrbs);
  mpi::bcast(*myComm, nbands);

  //buffer to serialize BandInfo
  PooledData<OHMMS_PRECISION_FULL> misc(nbands[0] * 5);
  bool isCore = false;
  n = NumDistinctOrbitals = nbands[1];
  NumValenceOrbs          = nbands[2];
  NumCoreOrbs             = nbands[3];

  if (root)
  {
    misc.rewind();
    //misc.put(NumValenceOrbs);
    //misc.put(NumCoreOrbs);
    for (int i = 0; i < n; ++i)
    {
      misc.put(SortBands[i].TwistIndex);
      misc.put(SortBands[i].BandIndex);
      misc.put(SortBands[i].Energy);
      misc.put(SortBands[i].MakeTwoCopies);
      misc.put(SortBands[i].IsCoreState);

      isCore |= SortBands[i].IsCoreState;
    }

    for (int i = n; i < SortBands.size(); ++i)
    {
      misc.put(SortBands[i].TwistIndex);
      misc.put(SortBands[i].BandIndex);
      misc.put(SortBands[i].Energy);
      misc.put(SortBands[i].MakeTwoCopies);
      misc.put(SortBands[i].IsCoreState);
    }
  }
  myComm->bcast(misc);

  if (!root)
  {
    SortBands.resize(nbands[0]);
    misc.rewind();
    //misc.get(NumValenceOrbs);
    //misc.get(NumCoreOrbs);
    for (int i = 0; i < n; ++i)
    {
      misc.get(SortBands[i].TwistIndex);
      misc.get(SortBands[i].BandIndex);
      misc.get(SortBands[i].Energy);
      misc.get(SortBands[i].MakeTwoCopies);
      misc.get(SortBands[i].IsCoreState);

      isCore |= SortBands[i].IsCoreState;
    }
    for (int i = n; i < SortBands.size(); ++i)
    {
      misc.get(SortBands[i].TwistIndex);
      misc.get(SortBands[i].BandIndex);
      misc.get(SortBands[i].Energy);
      misc.get(SortBands[i].MakeTwoCopies);
      misc.get(SortBands[i].IsCoreState);
    }
  }

  //char fname[64];
  //sprintf(fname,"debug.%d",myComm->rank());
  //ofstream fout(fname);
  //fout.setf(std::ios::scientific, std::ios::floatfield);
  //fout.precision(12);
  //for(int i=0; i<misc.size();++i)
  //  fout << misc[i] << std::endl;
  return isCore;
}

} // namespace qmcplusplus
