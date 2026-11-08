//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWOrbitalSetBuilder.h
 * @brief Build plane-wave SPOs from orbital HDF data and supercell input.
 */
#ifndef QMCPLUSPLUS_PWORBITAL_BUILDER_H
#define QMCPLUSPLUS_PWORBITAL_BUILDER_H

#include "SPOSetBuilder.h"
#include "hdf/hdf_archive.h"
#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/PlaneWave/PWOrbitalSet.h"
#else
#include "QMCWaveFunctions/PlaneWave/PWRealOrbitalSet.h"
#endif

namespace qmcplusplus
{
struct PWParameterSet;
class SlaterDet;

/** Build ordered SPOs from primitive k-points folding to a selected supercell twist.
 *
 * tilematrix maps the primitive lattice to the target cell. twist specifies a
 * reduced supercell vector; twistnum and twistIndex select a folded-group index,
 * with Gamma as the default. Ordering uses the same sort modes as the spline
 * builder. Real builds retain time-reversal representatives and project real
 * and imaginary orbitals as needed.
 *
 * Consecutive selected bands from one primitive k-point form a component SPO.
 * Multiple components are returned in a CompositeSPOSet without changing band
 * order. HDF access and construction run on every rank of the supplied communicator.
 */
class PWOrbitalSetBuilder : public SPOSetBuilder
{
private:
#if defined(QMC_COMPLEX)
  using SPOSetType = PWOrbitalSet;
#else
  using SPOSetType = PWRealOrbitalSet;
#endif

  /// Target particle set, whose twist is updated to the selected supercell twist.
  ParticleSet& targetPtcl;
  /// Root collection XML used when reapplying legacy HDF name parameters.
  xmlNodePtr rootNode{nullptr};
  /// Reduced primitive twist of the first contributing k-point.
  PosType TwistAngle;
  /// Selected canonical reduced supercell twist.
  TinyVector<double, OHMMS_DIM> superTwist;
  /// Explicit reduced supercell vector from the twist attribute.
  TinyVector<double, OHMMS_DIM> requestedTwist;
  ///true when the twist XML attribute contains a reduced-coordinate twist
  bool hasRequestedTwist{false};
  /// Band ordering: 0 preserves input order, 1 uses energy, 2 uses band index.
  int sortBands{1};
  ///number of spin channels in the orbital HDF file
  int numSpins{0};
  /// Integer primitive-to-supercell tiling matrix, defaulting to the identity.
  Tensor<int, OHMMS_DIM> tileMatrix;
  ///primitive-cell k-points contributing to the selected supercell twist
  std::vector<int> includedKPoints;
  ///primitive-cell twist angles read from the HDF file
  std::vector<TinyVector<double, OHMMS_DIM>> twistAngles;
  /// Legacy HDF naming parameters and the selected twist-group index.
  std::unique_ptr<PWParameterSet> myParam;
  /// Cached basis for the first contributing primitive k-point.
  std::unique_ptr<PWBasis> myBasisSet;
  /// Communicator-associated, independently read orbital archive.
  hdf_archive hfile;

public:
  /** Read orbital metadata, validate tiling, and prepare the first primitive basis.
   * @param p target particle set with the simulation-cell lattice
   * @param comm communicator whose ranks construct the same collection
   * @param cur collection XML containing href and optional tiling/twist attributes
   * @throws UniformCommunicateError for communicator-wide input validation failures
   * @throws std::runtime_error for lower-level HDF read failures
   */
  PWOrbitalSetBuilder(ParticleSet& p, Communicate& comm, xmlNodePtr cur);
  ~PWOrbitalSetBuilder() override;

  /** Construct the SPO described by a child sposet element.
   * size must be positive and spindataset must index an available HDF spin channel.
   * @throws UniformCommunicateError for invalid size, spin, occupation, or band counts
   */
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

private:
  /// Open the named XML HDF attribute and apply version-dependent naming parameters.
  bool getH5(xmlNodePtr cur, const char* aname);
  /// Read twist groups and lattice data, validate them, and set targetPtcl's twist.
  bool createPWBasis();
  /** Order and select bands, returning one primitive SPO or a composite of runs.
   * Transform and excited-occupation requests require a single contributing k-point.
   */
  std::unique_ptr<SPOSet> createPW(xmlNodePtr cur, const std::string& objname, int spinIndex, int orbital_count);
#if defined(QMC_COMPLEX)
  /// Legacy complex-orbital grid sampling path that writes spline data through hfile.
  void transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc);
#endif
};
} // namespace qmcplusplus
#endif
