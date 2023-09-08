
#include "QMCWaveFunctions/EinsplineSetBuilderT.h"

#include "CPU/SIMD/vmath.hpp"
#include "CPU/e2iphi.h"
#include "CPU/math.hpp"
#include "Message/CommOperators.h"
#include "Message/Communicate.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTableT.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBaseT.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSetT.h"
#include "QMCWaveFunctions/BsplineFactory/createBsplineReaderT.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/einspline_helper.hpp"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/Timer.h"
#include "Utilities/qmc_common.h"
#include <Message/UniformCommunicateError.h>
#include <PlatformSelector.hpp>
#include <fftw3.h>

#include <array>
#include <string_view>
#include <vector>

namespace qmcplusplus
{
// std::map<H5OrbSet,SPOSet*,H5OrbSet>  EinsplineSetBuilder::SPOSetMap;
// std::map<TinyVector<int,4>,EinsplineSetBuilder::OrbType*,Int4less>
// EinsplineSetBuilder::OrbitalMap;
////std::map<H5OrbSet,multi_UBspline_3d_z*,H5OrbSet>
/// EinsplineSetBuilder::ExtendedMap_z;
////std::map<H5OrbSet,multi_UBspline_3d_d*,H5OrbSet>
/// EinsplineSetBuilder::ExtendedMap_d;

template<typename T>
EinsplineSetBuilderT<T>::EinsplineSetBuilderT(ParticleSetT<T>& p,
                                              const PSetMap& psets,
                                              Communicate* comm,
                                              xmlNodePtr cur)
    : SPOSetBuilderT<T>("spline", comm),
      ParticleSets(psets),
      TargetPtcl(p),
      XMLRoot(cur),
      Format(QMCPACK),
      NumBands(0),
      NumElectrons(0),
      NumSpins(0),
      NumTwists(0),
      MeshFactor(1.0),
      MeshSize(0, 0, 0),
      twist_num_(-1),
      LastSpinSet(-1),
      NumOrbitalsRead(-1),
      makeRotations(false)
{
  this->ClassName = "EinsplineSetBuilder";

  MatchingTol = 10 * std::numeric_limits<float>::epsilon();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      TileMatrix(i, j) = 0;

  // invalidate states by the basis class
  this->states.clear();
  this->states.resize(p.groups());

  // create vectors with nullptr
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

template<typename T>
EinsplineSetBuilderT<T>::~EinsplineSetBuilderT()
{
  DEBUG_MEMORY("EinsplineSetBuilder::~EinsplineSetBuilder");
}

template<typename T>
bool EinsplineSetBuilderT<T>::CheckLattice()
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

template<typename T>
void EinsplineSetBuilderT<T>::BroadcastOrbitalInfo()
{
  if (this->myComm->size() == 1)
    return;
  int numIons         = IonTypes.size();
  int numDensityGvecs = TargetPtcl.DensityReducedGvecs.size();
  PooledData<double> abuffer;
  PooledData<int> aibuffer;
  aibuffer.add(Version.begin(), Version.end()); // myComm->bcast(Version);
  aibuffer.add(Format);
  abuffer.add(Lattice.begin(), Lattice.end()); // myComm->bcast(Lattice);
  abuffer.add(RecipLattice.begin(),
              RecipLattice.end()); // myComm->bcast(RecipLattice);
  abuffer.add(SuperLattice.begin(),
              SuperLattice.end());                   // myComm->bcast(SuperLattice);
  abuffer.add(LatticeInv.begin(), LatticeInv.end()); // myComm->bcast(LatticeInv);
  aibuffer.add(NumBands);                            // myComm->bcast(NumBands);
  aibuffer.add(NumElectrons);                        // myComm->bcast(NumElectrons);
  aibuffer.add(NumSpins);                            // myComm->bcast(NumSpins);
  aibuffer.add(NumTwists);                           // myComm->bcast(NumTwists);
  aibuffer.add(numIons);                             // myComm->bcast(numIons);
  aibuffer.add(numDensityGvecs);
  aibuffer.add(HaveOrbDerivs);
  this->myComm->bcast(abuffer);
  this->myComm->bcast(aibuffer);
  if (this->myComm->rank())
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
    aibuffer.get(numDensityGvecs);
    aibuffer.get(HaveOrbDerivs);
    TargetPtcl.DensityReducedGvecs.resize(numDensityGvecs);
    TargetPtcl.Density_G.resize(numDensityGvecs);
  }
  if (IonTypes.size() != numIons)
  {
    IonTypes.resize(numIons);
    IonPos.resize(numIons);
  }
  // new buffer
  PooledData<double> bbuffer;
  PooledData<int> bibuffer;
  for (int i = 0; i < numIons; ++i)
    bibuffer.add(IonTypes[i]);
  // myComm->bcast(IonTypes);
  bbuffer.add(&IonPos[0][0], &IonPos[0][0] + OHMMS_DIM * numIons);
  // myComm->bcast(IonPos);
  if (primcell_kpoints.size() != NumTwists)
    primcell_kpoints.resize(NumTwists);
  bbuffer.add(&primcell_kpoints[0][0], &primcell_kpoints[0][0] + OHMMS_DIM * NumTwists);
  bibuffer.add(&(TargetPtcl.DensityReducedGvecs[0][0]),
               &(TargetPtcl.DensityReducedGvecs[0][0]) + numDensityGvecs * OHMMS_DIM);
  bbuffer.add(&(TargetPtcl.Density_G[0]), &(TargetPtcl.Density_G[0]) + numDensityGvecs);
  this->myComm->bcast(bbuffer);
  this->myComm->bcast(bibuffer);
  if (this->myComm->rank())
  {
    bbuffer.rewind();
    bibuffer.rewind();
    for (int i = 0; i < numIons; ++i)
      bibuffer.get(IonTypes[i]);
    bbuffer.get(&IonPos[0][0], &IonPos[0][0] + OHMMS_DIM * numIons);
    bbuffer.get(&primcell_kpoints[0][0], &primcell_kpoints[0][0] + OHMMS_DIM * NumTwists);
    bibuffer.get(&(TargetPtcl.DensityReducedGvecs[0][0]),
                 &(TargetPtcl.DensityReducedGvecs[0][0]) + numDensityGvecs * OHMMS_DIM);
    bbuffer.get(&(TargetPtcl.Density_G[0]), &(TargetPtcl.Density_G[0]) + numDensityGvecs);
  }
  // buffer to bcast hybrid representation atomic orbital info
  PooledData<double> cbuffer;
  PooledData<int> cibuffer;
  this->myComm->bcast(cbuffer);
  this->myComm->bcast(cibuffer);
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
  this->myComm->bcast(cbuffer);
  this->myComm->bcast(cibuffer);
  if (this->myComm->rank())
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
// void
// EinsplineSetBuilder::CreateIonParticleSet( std::string sourceName)
//{
//   //    ParticleSet &pTemp = *(new MCWalkerConfiguration);
//   ParticleSet &pTemp = *(new ParticleSet);
//   pTemp.setName (sourceName);
//   SpeciesSet& tspecies(pTemp.getSpeciesSet());
//   ParticleSets[sourceName] = &pTemp;
// }
//

template<typename T>
void EinsplineSetBuilderT<T>::TileIons()
{
  // set the primitive lattice
  SourcePtcl->getPrimitiveLattice().set(Lattice);

  for (int j = 0; j < IonPos.size(); ++j)
    IonPos[j] = FracPart(SourcePtcl->getPrimitiveLattice().toUnit(IonPos[j]));

  IonPos.resize(SourcePtcl->getTotalNum());
  IonTypes.resize(SourcePtcl->getTotalNum());
  std::copy(SourcePtcl->R.begin(), SourcePtcl->R.end(), IonPos.begin());
  std::copy(SourcePtcl->GroupID.begin(), SourcePtcl->GroupID.end(), IonTypes.begin());

  // app_log() << "  Primitive Cell\n";
  // SourcePtcl->getPrimitiveLattice().print(app_log());
  // app_log() << "  Super Cell\n";
  // SourcePtcl->Lattice.print(app_log());

  // Don't need to do this, already one by ParticleSetPool.cpp
  //   Vector<TinyVector<double, OHMMS_DIM> > primPos   = IonPos;
  //   Vector<int>                            primTypes = IonTypes;
  //   int numCopies = std::abs(det(TileMatrix));
  //   IonTypes.resize(primPos.size()*numCopies);
  //   IonPos.resize  (primPos.size()*numCopies);
  //   int maxCopies = 10;
  //   using Vec3 = TinyVector<double,3>;
  //   int index=0;
  //   for (int i0=-maxCopies; i0<=maxCopies; i0++)
  //     for (int i1=-maxCopies; i1<=maxCopies; i1++)
  //       for (int i2=-maxCopies; i2<=maxCopies; i2++)
  //         for (int iat=0; iat < primPos.size(); iat++)
  //         {
  //           Vec3 r     = primPos[iat];
  //           Vec3 uPrim = PrimCell.toUnit(r);
  //           for (int i=0; i<3; i++)
  //             uPrim[i] -= std::floor(uPrim[i]);
  //           r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0) +
  //               (double)i1*PrimCell.a(1) + (double)i2*PrimCell.a(2);
  //           Vec3 uSuper = SuperCell.toUnit(r);
  //           if ((uSuper[0] >= -1.0e-4) && (uSuper[0] < 0.9999) &&
  //               (uSuper[1] >= -1.0e-4) && (uSuper[1] < 0.9999) &&
  //               (uSuper[2] >= -1.0e-4) && (uSuper[2] < 0.9999))
  //           {
  //             IonPos[index]= r;
  //             IonTypes[index]= primTypes[iat];
  //             index++;
  //           }
  //         }
  //   if (index != primPos.size()*numCopies)
  //   {
  //     app_error() << "The number of tiled ions, " << IonPos.size()
  //                 << ", does not match the expected number of "
  //                 << primPos.size()*numCopies << " or the index "<< index
  //                 <<".  Aborting.\n";
  //     APP_ABORT("EinsplineSetBuilder::TileIons()");
  //   }
  //   if (myComm->rank() == 0)
  //   {
  //     char buf[1000];
  //     snprintf (buf, 1000, "Supercell reduced ion positions = \n");
  //     app_log() << buf;
  //     app_log().flush();
  //     for (int i=0; i<IonPos.size(); i++)
  //     {
  //       PosType u = SuperCell.toUnit(IonPos[i]);
  //       char buf2[1000];
  //       snprintf (buf2, 1000, "   %14.10f %14.10f %14.10f\n",
  //                u[0], u[1], u[2]);
  //       app_log() << buf2;
  //       app_log().flush();
  //       //		 IonPos[i][0], IonPos[i][1], IonPos[i][2]);
  //     }
  //   }
}

template<typename T>
bool EinsplineSetBuilderT<T>::TwistPair(PosType a, PosType b) const
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

template<typename T>
void EinsplineSetBuilderT<T>::AnalyzeTwists2(const int twist_num_inp, const TinyVector<double, OHMMS_DIM>& twist_inp)
{
  Tensor<double, 3> S;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      S(i, j) = (double)TileMatrix(i, j);

  const int num_prim_kpoints = primcell_kpoints.size();

  // build a list of unique super twists that all the primitive cell k-point
  // correspond to.
  std::vector<PosType> superFracs; // twist super twist coordinates
  std::vector<int> superIndex;     // the indices of the super twists that correpsond to all
                                   // the primitive cell k-points in the unique list.
  {
    // scan all the primitive cell k-points
    for (int ki = 0; ki < num_prim_kpoints; ki++)
    {
      PosType primTwist  = primcell_kpoints[ki];
      PosType superTwist = dot(S, primTwist);
      PosType kp         = PrimCell.k_cart(primTwist);
      PosType ks         = SuperCell.k_cart(superTwist);
      // check the consistency of tiling, primitive and super cells.
      if (dot(ks - kp, ks - kp) > 1.0e-6)
      {
        app_error() << "Primitive and super k-points do not agree.  "
                       "Error in coding.\n";
        APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
      }
      PosType frac = FracPart(superTwist);
      // verify if the super twist that correpsonds to this primitive cell
      // k-point exists in the unique list or not.
      bool found = false;
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
    assert(superIndex.size() == num_prim_kpoints);
  }

  const int numSuperTwists = superFracs.size();
  {
    app_log() << "Found " << numSuperTwists << " distinct supercell twist" << (numSuperTwists > 1 ? "s" : "")
              << " based on " << num_prim_kpoints << " primitive cell k-point" << (num_prim_kpoints > 1 ? "s" : "")
              << std::endl;
    if (this->myComm->rank() == 0)
    {
      int n_tot_irred(0);
      for (int si = 0; si < numSuperTwists; si++)
      {
        std::array<char, 1000> buf;
        int length = std::snprintf(buf.data(), buf.size(), "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n", si,
                                   superFracs[si][0], superFracs[si][1], superFracs[si][2]);
        if (length < 0)
          throw std::runtime_error("Error converting Super twist to a string");
        app_log() << std::string_view(buf.data(), length);
        app_log().flush();
      }
    }
  }

  // For each supercell twist, create a list of primitive twists which
  // correspond to it.
  std::vector<std::vector<int>> superSets;
  {
    superSets.resize(numSuperTwists);
    for (int ki = 0; ki < num_prim_kpoints; ki++)
      superSets[superIndex[ki]].push_back(ki);
  }

  { // look up a super cell twist and return its index in the unique list of
    // super cell twists.
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
        std::array<char, 1000> buf;
        int length = std::snprintf(buf.data(), buf.size(),
                                   "AnalyzeTwists2. Input twist [ %9.5f %9.5f %9.5f] not "
                                   "found in the list of super twists above.\n",
                                   twist[0], twist[1], twist[2]);
        if (length < 0)
          throw std::runtime_error("Error generating error message");
        throw UniformCommunicateError(buf.data());
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
      app_warning() << "twist attribute does't exist but twistnum "
                       "attribute was found. "
                    << "This is potentially ambiguous. Specifying twist "
                       "attribute is preferred."
                    << std::endl;
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

    std::array<char, 1000> buf;
    int length = std::snprintf(buf.data(), buf.size(), "  Using supercell twist %d:  [ %9.5f %9.5f %9.5f]", twist_num_,
                               superFracs[twist_num_][0], superFracs[twist_num_][1], superFracs[twist_num_][2]);
    if (length < 0)
      throw std::runtime_error("Error converting supercell twist to a string");
    app_log() << std::string_view(buf.data(), length) << std::endl;
  }

  TargetPtcl.setTwist(superFracs[twist_num_]);
  if constexpr (!IsComplex_t<T>{}())
  {
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
  }
  // Now check to see that each supercell twist has the right twists
  // to tile the primitive cell orbitals.
  const int numTwistsNeeded = std::abs(det(TileMatrix));
  for (int si = 0; si < numSuperTwists; si++)
  {
    // First make sure we have enough points
    if (superSets[si].size() != numTwistsNeeded)
    {
      std::array<char, 1000> buf;
      int length = std::snprintf(buf.data(), buf.size(), "Super twist %d should own %d k-points, but owns %d.\n", si,
                                 numTwistsNeeded, static_cast<int>(superSets[si].size()));
      if (length < 0)
        throw std::runtime_error("Error generating Super twist string");
      app_error() << std::string_view(buf.data(), length);
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
      PosType twistPrim_i  = primcell_kpoints[superSets[si][i]];
      PosType twistSuper_i = dot(S, twistPrim_i);
      PosType superInt_i   = IntPart(twistSuper_i);
      for (int j = i + 1; j < N; j++)
      {
        PosType twistPrim_j  = primcell_kpoints[superSets[si][j]];
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
  if constexpr (!IsComplex_t<T>{}())
  {
    std::vector<int> copyTwists;
    for (int i = 0; i < IncludeTwists.size(); i++)
    {
      int ti          = IncludeTwists[i];
      PosType twist_i = primcell_kpoints[ti];
      bool distinct   = true;
      for (int j = i + 1; j < IncludeTwists.size(); j++)
      {
        int tj          = IncludeTwists[j];
        PosType twist_j = primcell_kpoints[tj];
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
      PosType twist_i  = primcell_kpoints[ti];
      for (int j = 0; j < copyTwists.size(); j++)
      {
        int tj          = copyTwists[j];
        PosType twist_j = primcell_kpoints[tj];
        if (TwistPair(twist_i, twist_j))
          MakeTwoCopies[i] = true;
      }
      if (this->myComm->rank() == 0)
      {
        std::array<char, 1000> buf;
        int length = std::snprintf(buf.data(), buf.size(), "Using %d copies of twist angle [%6.3f, %6.3f, %6.3f]\n",
                                   MakeTwoCopies[i] ? 2 : 1, twist_i[0], twist_i[1], twist_i[2]);
        if (length < 0)
          throw std::runtime_error("Error generating string");
        app_log() << std::string_view(buf.data(), length);
        app_log().flush();
      }
    }
    // Find out if we can make real orbitals
    use_real_splines_ = true;
    for (int i = 0; i < DistinctTwists.size(); i++)
    {
      int ti        = DistinctTwists[i];
      PosType twist = primcell_kpoints[ti];
      for (int j = 0; j < OHMMS_DIM; j++)
        if (std::abs(twist[j] - 0.0) > MatchingTol && std::abs(twist[j] - 0.5) > MatchingTol &&
            std::abs(twist[j] + 0.5) > MatchingTol)
          use_real_splines_ = false;
    }
    if (use_real_splines_ && (DistinctTwists.size() > 1))
    {
      app_log() << "***** Use of real orbitals is possible, but not "
                   "currently implemented\n"
                << "      with more than one twist angle.\n";
      use_real_splines_ = false;
    }
    if (use_real_splines_)
      app_log() << "Using real splines.\n";
    else
      app_log() << "Using complex splines.\n";
  }
  else
  {
    DistinctTwists.resize(IncludeTwists.size());
    MakeTwoCopies.resize(IncludeTwists.size());
    for (int i = 0; i < IncludeTwists.size(); i++)
    {
      DistinctTwists[i] = IncludeTwists[i];
      MakeTwoCopies[i]  = false;
    }
    use_real_splines_ = false;
  }
}

template<typename T>
void EinsplineSetBuilderT<T>::OccupyBands(int spin, int sortBands, int numOrbs, bool skipChecks)
{
  if (this->myComm->rank() != 0)
    return;
  if (spin >= NumSpins && !skipChecks)
  {
    app_error() << "To developer: User is requesting for orbitals in an "
                   "invalid spin group "
                << spin << ". Current h5 file only contains spin groups "
                << "[0.." << NumSpins - 1 << "]." << std::endl;
    app_error() << "To user: Orbital H5 file contains no spin down data "
                   "and is appropriate only for spin unpolarized "
                   "calculations. "
                << "If this is your intent, please replace 'spindataset=1' "
                   "with 'spindataset=0' in the input file."
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
      band.Energy = -1.01e100;
      H5File.read(band.Energy, ePath.str());
      if (band.Energy > -1.0e100)
      {
        H5File.read(band.Spin, sPath.str());
        if (band.Spin == spin)
          SortBands.push_back(band);
      }
    }
  }
  int orbIndex        = 0;
  int numOrbs_counter = 0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
}

template<typename T>
void EinsplineSetBuilderT<T>::bcastSortBands(int spin, int n, bool root)
{
  std::vector<BandInfo>& SortBands(*FullBands[spin]);

  TinyVector<int, 2> nbands(int(SortBands.size()), n);
  mpi::bcast(*this->myComm, nbands);

  // buffer to serialize BandInfo
  PooledData<OHMMS_PRECISION_FULL> misc(nbands[0] * 4);
  n = NumDistinctOrbitals = nbands[1];

  if (root)
  {
    misc.rewind();
    for (int i = 0; i < n; ++i)
    {
      misc.put(SortBands[i].TwistIndex);
      misc.put(SortBands[i].BandIndex);
      misc.put(SortBands[i].Energy);
      misc.put(SortBands[i].MakeTwoCopies);
    }

    for (int i = n; i < SortBands.size(); ++i)
    {
      misc.put(SortBands[i].TwistIndex);
      misc.put(SortBands[i].BandIndex);
      misc.put(SortBands[i].Energy);
      misc.put(SortBands[i].MakeTwoCopies);
    }
  }
  this->myComm->bcast(misc);

  if (!root)
  {
    SortBands.resize(nbands[0]);
    misc.rewind();
    for (int i = 0; i < n; ++i)
    {
      misc.get(SortBands[i].TwistIndex);
      misc.get(SortBands[i].BandIndex);
      misc.get(SortBands[i].Energy);
      misc.get(SortBands[i].MakeTwoCopies);
    }
    for (int i = n; i < SortBands.size(); ++i)
    {
      misc.get(SortBands[i].TwistIndex);
      misc.get(SortBands[i].BandIndex);
      misc.get(SortBands[i].Energy);
      misc.get(SortBands[i].MakeTwoCopies);
    }
  }
}

inline bool sortByIndex(BandInfo leftB, BandInfo rightB)
{
  if (leftB.BandIndex == rightB.BandIndex)
  {
    if ((leftB.Energy < rightB.Energy + 1e-6) && (leftB.Energy > rightB.Energy - 1e-6))
      return leftB.TwistIndex < rightB.TwistIndex;
    else
      return leftB.Energy < rightB.Energy;
  }
  else
    return (leftB.BandIndex < rightB.BandIndex);
};

template<typename T>
bool EinsplineSetBuilderT<T>::ReadOrbitalInfo_ESHDF(bool skipChecks)
{
  app_log() << "  Reading orbital file in ESHDF format.\n";
  H5File.read(Version, "/version");
  app_log() << "  ESHDF orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << std::endl;
  H5File.read(Lattice, "/supercell/primitive_vectors");
  RecipLattice = 2.0 * M_PI * inverse(Lattice);
  SuperLattice = dot(TileMatrix, Lattice);
  std::array<char, 1000> buff;
  int length = std::snprintf(buff.data(), buff.size(),
                             "  Lattice = \n    [ %9.6f %9.6f %9.6f\n"
                             "      %9.6f %9.6f %9.6f\n"
                             "      %9.6f %9.6f %9.6f ]\n",
                             Lattice(0, 0), Lattice(0, 1), Lattice(0, 2), Lattice(1, 0), Lattice(1, 1), Lattice(1, 2),
                             Lattice(2, 0), Lattice(2, 1), Lattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error converting lattice to a string");
  app_log() << std::string_view(buff.data(), length);
  length =
      std::snprintf(buff.data(), buff.size(),
                    "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
                    "      %9.6f %9.6f %9.6f\n"
                    "      %9.6f %9.6f %9.6f ]\n",
                    SuperLattice(0, 0), SuperLattice(0, 1), SuperLattice(0, 2), SuperLattice(1, 0), SuperLattice(1, 1),
                    SuperLattice(1, 2), SuperLattice(2, 0), SuperLattice(2, 1), SuperLattice(2, 2));
  if (length < 0)
    throw std::runtime_error("Error converting SuperLattice to a string");
  app_log() << std::string_view(buff.data(), length) << std::endl;
  if (!CheckLattice())
    throw std::runtime_error("CheckLattice failed");
  PrimCell.set(Lattice);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      LatticeInv(i, j) = RecipLattice(i, j) / (2.0 * M_PI);
  int have_dpsi = false;
  NumTwists = NumSpins = NumBands = 0;
  NumElectrons                    = TargetPtcl.getTotalNum();
  H5File.read(NumBands, "/electrons/kpoint_0/spin_0/number_of_states");
  H5File.readEntry(NumSpins, "/electrons/number_of_spins");
  H5File.read(NumTwists, "/electrons/number_of_kpoints");
  H5File.readEntry(have_dpsi, "/electrons/have_dpsi");
  HaveOrbDerivs = have_dpsi;
  app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons << ", spins=" << NumSpins << ", twists=" << NumTwists
            << std::endl;
  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  Vector<int> species_ids;
  H5File.read(species_ids, "/atoms/species_ids");
  int num_species;
  H5File.read(num_species, "/atoms/number_of_species");
  std::vector<int> atomic_numbers(num_species);
  for (int isp = 0; isp < num_species; isp++)
  {
    std::ostringstream name;
    name << "/atoms/species_" << isp << "/atomic_number";
    H5File.readEntry(atomic_numbers[isp], name.str());
  }
  IonTypes.resize(species_ids.size());
  for (int i = 0; i < species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
  H5File.read(IonPos, "/atoms/positions");
  for (int i = 0; i < IonTypes.size(); i++)
    app_log() << "Atom type(" << i << ") = " << IonTypes[i] << std::endl;
  /////////////////////////////////////
  // Read atom orbital info from xml //
  /////////////////////////////////////
  // construct Super2Prim mapping.
  if (Super2Prim.size() == 0)
  {
    // SourcePtcl->convert2Cart(SourcePtcl->R);
    Super2Prim.resize(SourcePtcl->R.size(), -1);
    std::vector<int> prim_atom_counts;
    prim_atom_counts.resize(IonPos.size(), 0);
    for (int i = 0; i < SourcePtcl->R.size(); i++)
    {
      PosType ref = PrimCell.toUnit_floor(SourcePtcl->R[i]);
      for (int j = 0; j < IonPos.size(); j++)
      {
        PosType dr = PrimCell.toUnit_floor(IonPos[j]) - ref;
        for (int k = 0; k < OHMMS_DIM; k++)
          dr[k] -= round(dr[k]);
        if (dot(dr, dr) < MatchingTol)
        {
          if (Super2Prim[i] < 0)
          {
            Super2Prim[i] = j;
            prim_atom_counts[j]++;
          }
          else
          {
            app_error() << "Supercell ion " << i << " at " << SourcePtcl->R[j]
                        << " was found twice in the primitive cell as ion " << Super2Prim[i] << " and " << j
                        << std::endl;
            if (!skipChecks)
              abort();
          }
        }
      }
      if (Super2Prim[i] < 0)
      {
        app_error() << "Supercell ion " << i << " not found in the primitive cell" << std::endl;
        if (!skipChecks)
          abort();
      }
      else
      {
        // app_log() << "Supercell ion " << i << " mapped to primitive
        // cell ion " << Super2Prim[i] << std::endl;
      }
    }
    const int tiling_size = std::abs(det(TileMatrix));
    for (int i = 0; i < IonPos.size(); i++)
      if (prim_atom_counts[i] != tiling_size)
      {
        app_error() << "Primitive cell ion " << i << " was found only " << prim_atom_counts[i]
                    << " times in the supercell rather than " << tiling_size << std::endl;
        if (!skipChecks)
          abort();
      }
    // construct AtomicCentersInfo
    AtomicCentersInfo.resize(IonPos.size());
    for (int i = 0; i < IonPos.size(); i++)
      AtomicCentersInfo.ion_pos[i] = IonPos[i];
    const auto& source_species = SourcePtcl->getSpeciesSet();
    int Zind                   = source_species.findAttribute("atomicnumber");
    const int table_id         = SourcePtcl->addTable(*SourcePtcl);
    const auto& ii_table       = SourcePtcl->getDistTable(table_id);
    SourcePtcl->update(true);
    for (int i = 0; i < IonPos.size(); i++)
    {
      AtomicCentersInfo.non_overlapping_radius[i] = std::numeric_limits<RealType>::max();
      // should only call get_first_neighbor to set non_overlapping_radius
      // if there are more than one atom  in the cell
      if (Super2Prim.size() == 1)
        continue;
      for (int j = 0; j < Super2Prim.size(); j++)
        if (Super2Prim[j] == i)
        {
          // set GroupID for each ion in primitive cell
          if ((Zind < 0) || (source_species(Zind, SourcePtcl->GroupID[j]) == IonTypes[i]))
            AtomicCentersInfo.GroupID[i] = SourcePtcl->GroupID[j];
          else
          {
            app_error() << "Primitive cell ion " << i << " vs supercell ion " << j
                        << " atomic number not matching: " << IonTypes[i] << " vs "
                        << source_species(Zind, SourcePtcl->GroupID[j]) << std::endl;
            if (!skipChecks)
              abort();
          }
          // set non_overlapping_radius for each ion in primitive cell
          RealType r(0);
          PosType dr;
          ii_table.get_first_neighbor(j, r, dr, false);
          if (r < 1e-3)
            APP_ABORT("EinsplineSetBuilder::ReadOrbitalInfo_ESHDF "
                      "too close ions <1e-3 bohr!");
          AtomicCentersInfo.non_overlapping_radius[i] = 0.5 * r;
          break;
        }
    }

    // load cutoff_radius, spline_radius, spline_npoints, lmax if exists.
    const int inner_cutoff_ind   = source_species.findAttribute("inner_cutoff");
    const int cutoff_radius_ind  = source_species.findAttribute("cutoff_radius");
    const int spline_radius_ind  = source_species.findAttribute("spline_radius");
    const int spline_npoints_ind = source_species.findAttribute("spline_npoints");
    const int lmax_ind           = source_species.findAttribute("lmax");

    for (int center_idx = 0; center_idx < AtomicCentersInfo.Ncenters; center_idx++)
    {
      const int my_GroupID = AtomicCentersInfo.GroupID[center_idx];
      if (inner_cutoff_ind >= 0)
        AtomicCentersInfo.inner_cutoff[center_idx] = source_species(inner_cutoff_ind, my_GroupID);
      if (cutoff_radius_ind >= 0)
        AtomicCentersInfo.cutoff[center_idx] = source_species(cutoff_radius_ind, my_GroupID);
      if (spline_radius_ind >= 0)
        AtomicCentersInfo.spline_radius[center_idx] = source_species(spline_radius_ind, my_GroupID);
      if (spline_npoints_ind >= 0)
        AtomicCentersInfo.spline_npoints[center_idx] = source_species(spline_npoints_ind, my_GroupID);
      if (lmax_ind >= 0)
        AtomicCentersInfo.lmax[center_idx] = source_species(lmax_ind, my_GroupID);
    }
  }
  ///////////////////////////
  // Read the twist angles //
  ///////////////////////////
  primcell_kpoints.resize(NumTwists);
  for (int ti = 0; ti < NumTwists; ti++)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    TinyVector<double, OHMMS_DIM> primcell_kpoints_DP;
    H5File.read(primcell_kpoints_DP, path.str());
    primcell_kpoints[ti] = primcell_kpoints_DP;
  }
  if (qmc_common.use_density)
  {
    //////////////////////////////////////////////////////////
    // Only if it is bulk: If the density has not been set in TargetPtcl,
    // and   // the density is available, read it in and save it     // in
    // TargetPtcl.                                       //
    //////////////////////////////////////////////////////////
    if (TargetPtcl.getLattice().SuperCellEnum == SUPERCELL_BULK)
    {
      // FIXME:  add support for more than one spin density
      if (TargetPtcl.Density_G.empty())
      {
        Array<double, OHMMS_DIM> Density_r_DP;
        TinyVector<int, 3> mesh;
        H5File.read(TargetPtcl.DensityReducedGvecs, "/electrons/density/gvectors");
        int numG = TargetPtcl.DensityReducedGvecs.size();
// Convert primitive G-vectors to supercell G-vectors
// Also, flip sign since ESHDF format uses opposite sign convention
#pragma omp parallel for
        for (int iG = 0; iG < numG; iG++)
          TargetPtcl.DensityReducedGvecs[iG] = -1 * dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
        app_log() << "  Read " << numG << " density G-vectors.\n";
        for (int ispin = 0; ispin < NumSpins; ispin++)
        {
          std::ostringstream density_r_path, density_g_path;
          density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
          density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
          H5File.readEntry(Density_r_DP, density_r_path.str());
          TargetPtcl.Density_r = Density_r_DP;
          if (TargetPtcl.DensityReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found density in "
                         "the HDF5 file.\n";
            std::vector<ComplexType> density_G;
            std::vector<std::complex<double>> Density_G_DP;
            H5File.read(Density_G_DP, density_g_path.str());
            density_G.assign(Density_G_DP.begin(), Density_G_DP.end());
            if (!density_G.size())
            {
              app_error() << "  Density reduced G-vectors "
                             "defined, but not the"
                          << " density.\n";
              abort();
            }
            else
            {
              if (ispin == 0)
                TargetPtcl.Density_G = density_G;
              else
                for (int iG = 0; iG < density_G.size(); iG++)
                  TargetPtcl.Density_G[iG] += density_G[iG];
            }
          }
        }
      }
      //////////////////////////////////////////////////////////
      // If the density has not been set in TargetPtcl, and   //
      // the density is available, read it in and save it     //
      // in TargetPtcl.                                       //
      //////////////////////////////////////////////////////////
      // FIXME:  add support for more than one spin potential
      if (!TargetPtcl.VHXC_r[0].size())
      {
        TinyVector<int, 3> mesh;
        H5File.readEntry(TargetPtcl.VHXCReducedGvecs, "/electrons/VHXC/gvectors");
        int numG = TargetPtcl.VHXCReducedGvecs.size();
// Convert primitive G-vectors to supercell G-vectors
// Also, flip sign since ESHDF format uses opposite sign convention
#pragma omp parallel for
        for (int iG = 0; iG < numG; iG++)
          TargetPtcl.VHXCReducedGvecs[iG] = -1 * dot(TileMatrix, TargetPtcl.VHXCReducedGvecs[iG]);
        app_log() << "  Read " << numG << " VHXC G-vectors.\n";
        for (int ispin = 0; ispin < NumSpins; ispin++)
        {
          Array<double, OHMMS_DIM> VHXC_r_DP;
          std::ostringstream VHXC_r_path, VHXC_g_path;
          VHXC_r_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_r";
          VHXC_g_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_g";
          H5File.readEntry(VHXC_r_DP, VHXC_r_path.str());
          TargetPtcl.VHXC_r[ispin] = VHXC_r_DP;
          if (TargetPtcl.VHXCReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found VHXC in the "
                         "HDF5 file.\n";
            std::vector<std::complex<double>> VHXC_G_DP;
            std::vector<ComplexType> VHXC_G;
            H5File.read(VHXC_G_DP, VHXC_g_path.str());
            VHXC_G.assign(VHXC_G_DP.begin(), VHXC_G_DP.end());
            if (!VHXC_G.size())
            {
              app_error() << "  VHXC reduced G-vectors defined, "
                             "but not the"
                          << " VHXC.\n";
              abort();
            }
            else
              TargetPtcl.VHXC_G[ispin] = VHXC_G;
          }
        }
      }
    }
  }
  else
  {
    app_log() << "   Skip initialization of the density" << std::endl;
  }
  return true;
}

template<typename T>
void EinsplineSetBuilderT<T>::OccupyBands_ESHDF(int spin, int sortBands, int numOrbs)
{
  if (this->myComm->rank() != 0)
    return;

  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  SortBands.clear(); //??? can exit if SortBands is already made?
  int maxOrbs(0);
  for (int ti = 0; ti < DistinctTwists.size(); ti++)
  {
    int tindex = DistinctTwists[ti];
    // First, read valence states
    std::ostringstream ePath;
    ePath << "/electrons/kpoint_" << tindex << "/spin_" << spin << "/eigenvalues";
    std::vector<double> eigvals;
    H5File.read(eigvals, ePath.str());
    for (int bi = 0; bi < NumBands; bi++)
    {
      BandInfo band;
      band.TwistIndex    = tindex;
      band.BandIndex     = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      band.Energy        = eigvals[bi];
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs += 2;
      else
        maxOrbs++;
    }
  }

  app_log() << SortBands.size() << " complex-valued orbitals supplied by h5 can be expanded up to " << maxOrbs
            << " SPOs." << std::endl;
  if (maxOrbs < numOrbs)
    this->myComm->barrier_and_abort("EinsplineSetBuilder::OccupyBands_ESHDF user input requests "
                                    "more orbitals than what the h5 file supplies.");

  // Now sort the bands by energy
  if (sortBands == 2)
  {
    app_log() << "Sorting the bands by index now:\n";
    sort(SortBands.begin(), SortBands.end(), sortByIndex);
  }
  else if (sortBands == 1)
  {
    app_log() << "Sorting the bands now:\n";
    sort(SortBands.begin(), SortBands.end());
  }

  std::vector<int> gsOcc(maxOrbs);
  int N_gs_orbs = numOrbs;
  int nocced(0);
  for (int ti = 0; ti < SortBands.size(); ti++)
  {
    if (nocced < N_gs_orbs)
    {
      if (SortBands[ti].MakeTwoCopies && (N_gs_orbs - nocced > 1))
      {
        nocced += 2;
        gsOcc[ti] = 2;
      }
      else if ((SortBands[ti].MakeTwoCopies && (N_gs_orbs - nocced == 1)) || !SortBands[ti].MakeTwoCopies)
      {
        nocced += 1;
        gsOcc[ti] = 1;
      }
    }
  }
  if (occ_format == "energy")
  {
    app_log() << "  Occupying bands based on energy in mode " << (Occ.size() > 0 ? "\"excited\"" : "\"ground\"")
              << std::endl;
    // To get the occupations right.
    std::vector<int> Removed(0, 0);
    std::vector<int> Added(0, 0);
    for (int ien = 0; ien < Occ.size(); ien++)
    {
      if (Occ[ien] < 0)
        Removed.push_back(-Occ[ien]);
      else if (Occ[ien] > 0)
        Added.push_back(Occ[ien]);
    }
    if (Added.size() - Removed.size() != 0)
    {
      app_log() << "need to add and remove same number of orbitals. " << Added.size() << " " << Removed.size()
                << std::endl;
      APP_ABORT("ChangedOccupations");
    }
    std::vector<int> DiffOcc(maxOrbs, 0);
    // Probably a cleaner way to do this.
    for (int i = 0; i < Removed.size(); i++)
      DiffOcc[Removed[i] - 1] -= 1;
    for (int i = 0; i < Added.size(); i++)
      DiffOcc[Added[i] - 1] += 1;
    std::vector<int> SumOrb(SortBands.size(), 0);
    int doi(0);
    for (int i = 0; i < SumOrb.size(); i++)
    {
      if (SortBands[i].MakeTwoCopies)
      {
        SumOrb[i] = gsOcc[i] + DiffOcc[doi++];
        SumOrb[i] += DiffOcc[doi++];
      }
      else
        SumOrb[i] = gsOcc[i] + DiffOcc[doi++];
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for (int i = 0; i < SumOrb.size(); i++)
    {
      if (SumOrb[i] == 2)
      {
        SortBands[i].MakeTwoCopies = true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i] == 1)
      {
        SortBands[i].MakeTwoCopies = false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i] == 0)
      {
        SortBands[i].MakeTwoCopies = false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log() << " Trying to add the same orbital (" << i << ") less than zero or more than 2 times." << std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(), RejectedBands.begin(), RejectedBands.end());
    SortBands = ReOrderedBands;
  }
  else if (occ_format == "band")
  {
    app_log() << "  Occupying bands based on (ti,bi) data." << std::endl;
    if (Occ.size() != particle_hole_pairs * 4)
    {
      app_log() << " Need Occ = pairs*4. Occ is (ti,bi) of removed, then added." << std::endl;
      app_log() << Occ.size() << " " << particle_hole_pairs << std::endl;
      APP_ABORT("ChangedOccupations");
    }
    int cnt(0);
    for (int ien = 0; ien < SortBands.size(); ien++)
    {
      if ((Occ[cnt] == SortBands[ien].TwistIndex) && (Occ[cnt + 1] == SortBands[ien].BandIndex))
      {
        if (cnt < particle_hole_pairs * 2)
        {
          gsOcc[ien] -= 1;
          cnt += 2;
          app_log() << "removing orbital " << ien << std::endl;
        }
        else
        {
          gsOcc[ien] += 1;
          app_log() << "adding orbital " << ien << std::endl;
          cnt += 2;
        }
      }
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for (int i = 0; i < SortBands.size(); i++)
    {
      if (gsOcc[i] == 2)
      {
        SortBands[i].MakeTwoCopies = true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i] == 1)
      {
        SortBands[i].MakeTwoCopies = false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i] == 0)
      {
        SortBands[i].MakeTwoCopies = false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log() << " Trying to add the same orbital (" << i << ") less than zero or more than 2 times." << std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(), RejectedBands.begin(), RejectedBands.end());
    SortBands = ReOrderedBands;
  }
  // for(int sw=0;sw<Removed.size();sw++){
  //   app_log()<<" Swapping two orbitals "<<Removed[sw]<<" and "<<Added[sw]<<
  //   std::endl; BandInfo tempband(SortBands[Removed[sw]-1]);
  //   SortBands[Removed[sw]-1] = SortBands[Added[sw]-1];
  //   SortBands[Added[sw]-1] = tempband;
  // }
  int orbIndex        = 0;
  int numOrbs_counter = 0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct complex-valued orbitals from h5.\n";
}

template<typename T>
void EinsplineSetBuilderT<T>::set_metadata(int numOrbs,
                                           int twist_num_inp,
                                           const TinyVector<double, OHMMS_DIM>& twist_inp,
                                           bool skipChecks)
{
  // 1. set a lot of internal parameters in the EinsplineSetBuilder class
  //  e.g. TileMatrix, use_real_splines_, DistinctTwists, MakeTwoCopies.
  // 2. this is also where metadata for the orbitals are read from the
  // wavefunction hdf5 file
  //  and broadcast to MPI groups. Variables broadcasted are listed in
  //  EinsplineSetBuilderCommon.cpp
  //  EinsplineSetBuilder::BroadcastOrbitalInfo()
  //

  Timer orb_info_timer;
  // The tiling can be set by a simple vector, (e.g. 2x2x2), or by a
  // full 3x3 matrix of integers.  If the tilematrix was not set in
  // the input file...
  bool matrixNotSet = true;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      matrixNotSet = matrixNotSet && (TileMatrix(i, j) == 0);
  // then set the matrix to identity.
  if (matrixNotSet)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        TileMatrix(i, j) = (i == j) ? 1 : 0;
  if (this->myComm->rank() == 0)
  {
    std::array<char, 1000> buff;
    int length = std::snprintf(buff.data(), buff.size(),
                               "  TileMatrix = \n [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d "
                               "]\n",
                               TileMatrix(0, 0), TileMatrix(0, 1), TileMatrix(0, 2), TileMatrix(1, 0), TileMatrix(1, 1),
                               TileMatrix(1, 2), TileMatrix(2, 0), TileMatrix(2, 1), TileMatrix(2, 2));
    if (length < 0)
      throw std::runtime_error("Error converting TileMatrix to a string");
    app_log() << std::string_view(buff.data(), length);
  }
  if (numOrbs == 0)
    this->myComm->barrier_and_abort("EinsplineSetBuilder::createSPOSet You must specify the number of "
                                    "orbitals in the input file.");
  else
    app_log() << "  Reading " << numOrbs << " orbitals from HDF5 file.\n";

  /////////////////////////////////////////////////////////////////
  // Read the basic orbital information, without reading all the //
  // orbitals themselves.                                        //
  /////////////////////////////////////////////////////////////////
  orb_info_timer.restart();
  if (this->myComm->rank() == 0)
    if (!ReadOrbitalInfo(skipChecks))
      throw std::runtime_error("EinsplineSetBuilder::set_metadata Error "
                               "reading orbital info from HDF5 file.");
  app_log() << "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  this->myComm->barrier();

  orb_info_timer.restart();
  BroadcastOrbitalInfo();
  app_log() << "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  app_log().flush();

  // setup primitive cell and supercell
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt = dot(transpose(PrimCell.G), PrimCell.G);

  // Now, analyze the k-point mesh to figure out the what k-points  are needed
  AnalyzeTwists2(twist_num_inp, twist_inp);
}

template<typename T>
std::unique_ptr<SPOSetT<T>> EinsplineSetBuilderT<T>::createSPOSetFromXML(xmlNodePtr cur)
{
  // use 2 bohr as the default when truncated orbitals are used based on the
  // extend of the ions
  int numOrbs = 0;
  int sortBands(1);
  int spinSet       = 0;
  bool skipChecks   = false;
  int twist_num_inp = TWISTNUM_NO_INPUT;
  TinyVector<double, OHMMS_DIM> twist_inp(TWIST_NO_INPUT);

  std::string sourceName;
  std::string spo_prec("double");
  std::string truncate("no");
  std::string hybrid_rep("no");
  std::string skip_checks("no");
  std::string use_einspline_set_extended("no"); // use old spline library for high-order derivatives, e.g. needed
                                                // for backflow optimization
  std::string useGPU;
  std::string GPUsharing = "no";
  std::string spo_object_name;

  ScopedTimer spo_timer_scope(createGlobalTimer("einspline::CreateSPOSetFromXML", timer_level_medium));

  {
    TinyVector<int, OHMMS_DIM> TileFactor_do_not_use;
    OhmmsAttributeSet a;
    a.add(H5FileName, "href");
    a.add(TileFactor_do_not_use, "tile", {}, TagStatus::DELETED);
    a.add(sortBands, "sort");
    a.add(TileMatrix, "tilematrix");
    a.add(twist_num_inp, "twistnum");
    a.add(twist_inp, "twist");
    a.add(sourceName, "source");
    a.add(MeshFactor, "meshfactor");
    a.add(hybrid_rep, "hybridrep");
    a.add(useGPU, "gpu", CPUOMPTargetSelector::candidate_values);
    a.add(GPUsharing,
          "gpusharing"); // split spline across GPUs visible per rank
    a.add(spo_prec, "precision");
    a.add(truncate, "truncate");
    a.add(this->myName, "tag");
    a.add(skip_checks, "skip_checks");

    a.put(XMLRoot);
    a.add(numOrbs, "size");
    a.add(numOrbs, "norbs");
    a.add(spinSet, "spindataset");
    a.add(spinSet, "group");
    a.put(cur);

    if (this->myName.empty())
      this->myName = "einspline";
  }

  if (skip_checks == "yes")
    skipChecks = true;

  auto pit(ParticleSets.find(sourceName));
  if (pit == ParticleSets.end())
    this->myComm->barrier_and_abort("Einspline needs the source particleset");
  else
    SourcePtcl = pit->second.get();

  ///////////////////////////////////////////////
  // Read occupation information from XML file //
  ///////////////////////////////////////////////
  const std::vector<int> last_occ(Occ);
  Occ.resize(0, 0); // correspond to ground
  bool NewOcc(false);

  {
    OhmmsAttributeSet oAttrib;
    oAttrib.add(spinSet, "spindataset");
    oAttrib.add(spo_object_name, "name");
    oAttrib.add(spo_object_name, "id");
    oAttrib.put(cur);
  }

  xmlNodePtr spo_cur = cur;
  cur                = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "occupation")
    {
      std::string occ_mode("ground");
      occ_format          = "energy";
      particle_hole_pairs = 0;
      OhmmsAttributeSet oAttrib;
      oAttrib.add(occ_mode, "mode");
      oAttrib.add(spinSet, "spindataset");
      oAttrib.add(occ_format, "format");
      oAttrib.add(particle_hole_pairs, "pairs");
      oAttrib.put(cur);
      if (occ_mode == "excited")
        putContent(Occ, cur);
      else if (occ_mode != "ground")
        this->myComm->barrier_and_abort("EinsplineSetBuilder::createSPOSet Only ground state "
                                        "occupation "
                                        "currently supported in EinsplineSetBuilder.");
    }
    cur = cur->next;
  }
  if (Occ != last_occ)
  {
    NewOcc = true;
  }
  else
    NewOcc = false;
#if defined(MIXED_PRECISION)
  app_log() << "\t  MIXED_PRECISION=1 Overwriting the einspline storage to "
               "single precision.\n";
  spo_prec = "single"; // overwrite
#endif
  H5OrbSet aset(H5FileName, spinSet, numOrbs);
  const auto iter = SPOSetMap.find(aset);
  if ((iter != SPOSetMap.end()) && (!NewOcc))
    app_warning() << "!!!!!!! Identical SPOSets are detected by EinsplineSetBuilder! "
                     "Implicit sharing one SPOSet for spin-up and spin-down "
                     "electrons has been removed. "
                     "Each determinant creates its own SPOSet with dedicated memory "
                     "for spline coefficients. "
                     "To avoid increasing the memory footprint of spline "
                     "coefficients, "
                     "create a single SPOset outside the determinantset using "
                     "'sposet_collection' "
                     "and reference it by name on the determinant line."
                  << std::endl;

  if (FullBands[spinSet] == 0)
    FullBands[spinSet] = std::make_unique<std::vector<BandInfo>>();

  // Ensure the first SPO set must be spinSet==0
  // to correctly initialize key data of EinsplineSetBuilder
  if (SPOSetMap.size() == 0 && spinSet != 0)
    this->myComm->barrier_and_abort("The first SPO set must have spindataset=\"0\"");

  // set the internal parameters
  if (spinSet == 0)
    set_metadata(numOrbs, twist_num_inp, twist_inp, skipChecks);

  //////////////////////////////////
  // Create the OrbitalSet object
  //////////////////////////////////
  Timer mytimer;
  mytimer.restart();
  OccupyBands(spinSet, sortBands, numOrbs, skipChecks);
  if (spinSet == 0)
    TileIons();

  bool use_single = (spo_prec == "single" || spo_prec == "float");

  // safeguard for a removed feature
  if (truncate == "yes")
    this->myComm->barrier_and_abort("The 'truncate' feature of spline SPO has been removed. Please use "
                                    "hybrid orbital representation.");

  createBsplineReader(use_single, hybrid_rep == "yes", useGPU);

  MixedSplineReader->setCommon(XMLRoot);
  // temporary disable the following function call, Ye Luo
  // RotateBands_ESHDF(spinSet,
  // dynamic_cast<EinsplineSetExtended<std::complex<double> >*>(OrbitalSet));
  bcastSortBands(spinSet, NumDistinctOrbitals, this->myComm->rank() == 0);
  auto OrbitalSet = MixedSplineReader->create_spline_set(spinSet, spo_cur);
  if (!OrbitalSet)
    this->myComm->barrier_and_abort("Failed to create SPOSet*");
  app_log() << "Time spent in creating B-spline SPOs " << mytimer.elapsed() << "sec" << std::endl;
  OrbitalSet->finalizeConstruction();
  SPOSetMap[aset] = OrbitalSet.get();
  return OrbitalSet;
}

template<typename T>
void EinsplineSetBuilderT<T>::createBsplineReader(bool useSingle, bool hybridRep, const std::string& useGPU)
{
  if (use_real_splines_)
  {
    // if(TargetPtcl.Lattice.SuperCellEnum != SUPERCELL_BULK &&
    // truncate=="yes")
    if (MixedSplineReader == 0)
    {
      if (useSingle)
        MixedSplineReader = createBsplineRealSingleT(this, hybridRep, useGPU);
      else
        MixedSplineReader = createBsplineRealDoubleT(this, hybridRep, useGPU);
    }
  }
  else
  {
    if (MixedSplineReader == 0)
    {
      if (useSingle)
        MixedSplineReader = createBsplineComplexSingleT(this, hybridRep, useGPU);
      else
        MixedSplineReader = createBsplineComplexDoubleT(this, hybridRep, useGPU);
    }
  }
}

#ifdef QMC_COMPLEX
template<>
void EinsplineSetBuilderT<std::complex<float>>::createBsplineReader(bool useSingle,
                                                                    bool hybridRep,
                                                                    const std::string& useGPU)
{
  if (MixedSplineReader == 0)
  {
    if (useSingle)
      MixedSplineReader = createBsplineComplexSingleT(this, hybridRep, useGPU);
    else
      MixedSplineReader = createBsplineComplexDoubleT(this, hybridRep, useGPU);
  }
}

template<>
void EinsplineSetBuilderT<std::complex<double>>::createBsplineReader(bool useSingle,
                                                                     bool hybridRep,
                                                                     const std::string& useGPU)
{
  if (MixedSplineReader == 0)
  {
    if (useSingle)
      MixedSplineReader = createBsplineComplexSingleT(this, hybridRep, useGPU);
    else
      MixedSplineReader = createBsplineComplexDoubleT(this, hybridRep, useGPU);
  }
}
#endif

template<typename T>
std::unique_ptr<SPOSetT<T>> EinsplineSetBuilderT<T>::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info)
{
  if (MixedSplineReader == 0)
    this->myComm->barrier_and_abort("EinsplineSetExtended<T> cannot create a SPOSet");

  std::string aname;
  int spinSet(0);
  OhmmsAttributeSet a;
  a.add(spinSet, "spindataset");
  a.add(spinSet, "group");
  a.put(cur);

  // allow only non-overlapping index sets and use the max index as the
  // identifier
  int norb = input_info.max_index();
  H5OrbSet aset(H5FileName, spinSet, norb);

  auto bspline_zd = MixedSplineReader->create_spline_set(spinSet, cur, input_info);
  if (bspline_zd)
    SPOSetMap[aset] = bspline_zd.get();
  return bspline_zd;
}

template<typename T>
bool EinsplineSetBuilderT<T>::ReadOrbitalInfo(bool skipChecks)
{
  if (!H5File.open(H5FileName, H5F_ACC_RDONLY))
  {
    app_error() << "Could not open HDF5 file \"" << H5FileName << "\" in EinsplineSetBuilder::ReadOrbitalInfo.\n";
    return false;
  }

  // Read format
  std::string format;
  H5File.read(format, "/format");
  H5File.read(Version, "/version");
  app_log() << "  HDF5 orbital file version " << Version[0] << "." << Version[1] << "." << Version[2] << "\n";
  if (format.find("ES") < format.size())
  {
    Format = ESHDF;
    return ReadOrbitalInfo_ESHDF(skipChecks);
  }

  app_error() << "EinsplineSetBuilder::ReadOrbitalInfo too old h5 file which "
                 "is not in ESHDF format! Regenerate the h5 file";
  return false;
}

template<typename T>
bool EinsplineSetBuilderT<T>::ReadGvectors_ESHDF()
{
  bool root = this->myComm->rank() == 0;
  // this is always ugly
  MeshSize    = 0;
  int hasPsig = 1;
  if (root)
  {
    H5File.readEntry(MeshSize, "/electrons/psi_r_mesh");
    H5File.readEntry(MeshSize, "/electrons/mesh");
  }
  this->myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
  if (hasPsig)
  {
    int nallowed  = 257;
    int allowed[] = {72,    75,    80,    81,    90,    96,    100,   108,   120,   125,   128,   135,   144,   150,
                     160,   162,   180,   192,   200,   216,   225,   240,   243,   250,   256,   270,   288,   300,
                     320,   324,   360,   375,   384,   400,   405,   432,   450,   480,   486,   500,   512,   540,
                     576,   600,   625,   640,   648,   675,   720,   729,   750,   768,   800,   810,   864,   900,
                     960,   972,   1000,  1024,  1080,  1125,  1152,  1200,  1215,  1250,  1280,  1296,  1350,  1440,
                     1458,  1500,  1536,  1600,  1620,  1728,  1800,  1875,  1920,  1944,  2000,  2025,  2048,  2160,
                     2187,  2250,  2304,  2400,  2430,  2500,  2560,  2592,  2700,  2880,  2916,  3000,  3072,  3125,
                     3200,  3240,  3375,  3456,  3600,  3645,  3750,  3840,  3888,  4000,  4050,  4096,  4320,  4374,
                     4500,  4608,  4800,  4860,  5000,  5120,  5184,  5400,  5625,  5760,  5832,  6000,  6075,  6144,
                     6250,  6400,  6480,  6561,  6750,  6912,  7200,  7290,  7500,  7680,  7776,  8000,  8100,  8192,
                     8640,  8748,  9000,  9216,  9375,  9600,  9720,  10000, 10125, 10240, 10368, 10800, 10935, 11250,
                     11520, 11664, 12000, 12150, 12288, 12500, 12800, 12960, 13122, 13500, 13824, 14400, 14580, 15000,
                     15360, 15552, 15625, 16000, 16200, 16384, 16875, 17280, 17496, 18000, 18225, 18432, 18750, 19200,
                     19440, 19683, 20000, 20250, 20480, 20736, 21600, 21870, 22500, 23040, 23328, 24000, 24300, 24576,
                     25000, 25600, 25920, 26244, 27000, 27648, 28125, 28800, 29160, 30000, 30375, 30720, 31104, 31250,
                     32000, 32400, 32768, 32805, 33750, 34560, 34992, 36000, 36450, 36864, 37500, 38400, 38880, 39366,
                     40000, 40500, 40960, 41472, 43200, 43740, 45000, 46080, 46656, 46875, 48000, 48600, 49152, 50000,
                     50625, 51200, 51840, 52488, 54000, 54675, 55296, 56250, 57600, 58320, 59049, 60000, 60750, 61440,
                     62208, 62500, 64000, 64800, 65536};
    MaxNumGvecs   = 0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int, 3> maxIndex(0, 0, 0);
    Gvecs.resize(NumTwists);
    {
      int numg = 0;
      if (root)
      {
        std::ostringstream Gpath;
        Gpath << "/electrons/kpoint_0/gvectors";
        H5File.read(Gvecs[0], Gpath.str());
        numg = Gvecs[0].size();
      }
      this->myComm->bcast(numg);
      if (!root)
        Gvecs[0].resize(numg);
      this->myComm->bcast(Gvecs[0]);
      MaxNumGvecs = Gvecs[0].size();
      for (int ig = 0; ig < Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } // done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0 * MeshFactor * maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0 * MeshFactor * maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0 * MeshFactor * maxIndex[2]);
    // only use 2^a 3^b 5^c where a>=2  up to 65536
    int* ix     = std::lower_bound(allowed, allowed + nallowed, MeshSize[0]);
    int* iy     = std::lower_bound(allowed, allowed + nallowed, MeshSize[1]);
    int* iz     = std::lower_bound(allowed, allowed + nallowed, MeshSize[2]);
    MeshSize[0] = (MeshSize[0] > 128) ? *ix : (MeshSize[0] + MeshSize[0] % 2);
    MeshSize[1] = (MeshSize[1] > 128) ? *iy : (MeshSize[1] + MeshSize[1] % 2);
    MeshSize[2] = (MeshSize[2] > 128) ? *iz : (MeshSize[2] + MeshSize[2] % 2);
    if (Version[0] < 2)
    {
      // get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. "
                   "Regenerate orbital files using updated QE."
                << std::endl;
      for (int k = 0; k < DistinctTwists.size(); ++k)
      {
        int ik = DistinctTwists[k];
        if (ik == 0)
          continue; // already done
        int numg = 0;
        if (root)
        {
          std::ostringstream Gpath;
          Gpath << "/electrons/kpoint_" << ik << "/gvectors";
          H5File.read(Gvecs[ik], Gpath.str());
          numg = Gvecs[ik].size();
        }
        this->myComm->bcast(numg);
        if (numg == 0)
        {
          // copy kpoint_0, default
          Gvecs[ik] = Gvecs[0];
        }
        else
        {
          if (numg != MaxNumGvecs)
          {
            std::ostringstream o;
            o << "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
              << " This is not supported anymore. Rerun "
                 "pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if (!root)
            Gvecs[ik].resize(numg);
          this->myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}

//#ifndef QMC_COMPLEX
template class EinsplineSetBuilderT<double>;
template class EinsplineSetBuilderT<float>;
#ifdef QMC_COMPLEX
template class EinsplineSetBuilderT<std::complex<double>>;
template class EinsplineSetBuilderT<std::complex<float>>;
#endif
} // namespace qmcplusplus
