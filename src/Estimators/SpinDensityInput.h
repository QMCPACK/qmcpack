//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
/** Native representation for Spin Density Estimators input parameters
 */
class SpinDensityInput
{
public:
  using Real               = QMCTraits::RealType;
  using POLT               = PtclOnLatticeTraits;
  using Lattice            = POLT::ParticleLayout_t;
  using PosType            = QMCTraits::PosType;
  static constexpr int DIM = QMCTraits::DIM;

public:
  SpinDensityInput(){};
  SpinDensityInput(Lattice lattice);
  void readXML(xmlNodePtr cur);
  Lattice get_cell() const { return cell_; }
  PosType get_corner() const { return corner_; }
  TinyVector<int, DIM> get_grid() const { return grid_; }
  TinyVector<int, DIM> get_gdims() const { return gdims_; }
  int get_npoints() const { return npoints_; }
  bool get_write_report() const { return write_report_; }

private:
  ///name of this Estimator
  std::string myName_;

  Lattice cell_;
  PosType corner_;
  TinyVector<int, DIM> grid_;
  // Striding of elements of the grid.
  TinyVector<int, DIM> gdims_;
  int npoints_;
  bool write_report_;
};

} // namespace qmcplusplus
#endif /* SPINDENSITYINPUT_H */
