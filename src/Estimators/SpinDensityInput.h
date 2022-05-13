//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SpinDensity.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SPINDENSITYINPUT_H
#define QMCPLUSPLUS_SPINDENSITYINPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Containers/OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{

class SpinDensityNew;

/** Native representation for Spin Density Estimators inputs
 *
 *  This class servers three purposes all related to properly handling
 *  and verifying the spin density input.
 *  1. An immutable representation of actual user input
 *  2. Parse the xml node of SpinDensityNew input.
 *  3. Hold the logic of calculating derived parameters.
 *
 */
class SpinDensityInput
{
public:
  using Real               = QMCTraits::RealType;
  using POLT               = PtclOnLatticeTraits;
  using Lattice            = POLT::ParticleLayout;
  using PosType            = QMCTraits::PosType;
  using Consumer = SpinDensityNew;
  static constexpr int DIM = QMCTraits::DIM;

public:
  SpinDensityInput(){}
  SpinDensityInput(xmlNodePtr node);
  void readXML(xmlNodePtr cur);
  Lattice get_cell() const { return cell_; }
  PosType get_corner() const { return corner_; }
  TinyVector<int, DIM> get_grid() const { return grid_; }
  int get_npoints() const { return npoints_; }
  bool get_write_report() const { return write_report_; }
  bool get_save_memory() const { return save_memory_; }

  struct DerivedParameters
  {
    PosType corner;
    TinyVector<int, DIM> grid;
    TinyVector<int, DIM> gdims;
    size_t npoints;
  };

  /** Derived parameters of SpinDensity
   *
   *  These require the cell the SpinDensity is evaluated over,
   *  the caller (SpinDensityNew) either gets this from the input and
   *  passes it back or passes in the cell from the relevant ParticleSet.
   *
   */
  DerivedParameters calculateDerivedParameters(const Lattice& lattice) const;

private:
  ///name of this Estimator
  std::string myName_;

  Lattice cell_;
  PosType corner_;
  PosType dr_;
  PosType center_;
  TinyVector<int, DIM> grid_;
  int npoints_;
  bool write_report_;
  bool save_memory_;
  /** these are necessary for calculateDerivedParameters
   *  
   *  If we are going to later write out a canonical input for
   *  this input then they are needed as well.
   */
  bool have_dr_     = false;
  bool have_grid_   = false;
  bool have_center_ = false;
  bool have_corner_ = false;
  bool have_cell_   = false;
};

} // namespace qmcplusplus
#endif /* SPINDENSITYINPUT_H */
