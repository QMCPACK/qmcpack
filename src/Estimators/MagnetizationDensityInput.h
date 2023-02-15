//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MAGNETIZATION_DENSITY_INPUT_H
#define QMCPLUSPLUS_MAGNETIZATION_DENSITY_INPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Containers/OhmmsPETE/TinyVector.h"
#include "Estimators/InputSection.h"
namespace qmcplusplus
{

class MagnetizationDensity;

class MagnetizationDensityInput
{
public:
  enum class Integrator
  {
    SIMPSONS,
    MONTECARLO
  };

  inline static const std::unordered_map<std::string, std::any>
      lookup_input_enum_value{{"integrator-simpsons", Integrator::SIMPSONS},
                              {"integrator-montecarlo", Integrator::MONTECARLO}};

  using Real               = QMCTraits::RealType;
  using POLT               = PtclOnLatticeTraits;
  using Lattice            = POLT::ParticleLayout;
  using PosType            = QMCTraits::PosType;
  using Consumer           = MagnetizationDensity;
  static constexpr int DIM = QMCTraits::DIM;

public:
  MagnetizationDensityInput(xmlNodePtr node);
  /** default copy constructor
   *  This is required due to SDI being part of a variant used as a vector element.
   */
  MagnetizationDensityInput(){};
  MagnetizationDensityInput(const MagnetizationDensityInput&) = default;
  Lattice get_cell() const { return cell_; }
  PosType get_corner() const { return corner_; }
  TinyVector<int, DIM> get_grid() const { return grid_; }
  int get_npoints() const { return npoints_; }
  int get_nsamples() const { return nsamples_; }
  Integrator get_integrator() const { return integrator_; }
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

  class MagnetizationDensityInputSection : public InputSection
  {
  public:
    /** parse time definition of input parameters */
    MagnetizationDensityInputSection()
    {
      // clang-format off
      section_name  = "MagnetizationDensity";
      attributes    = {"name", "type"};
      parameters    = {"integrator", 
                       "corner", "center", "samples", "grid", "dr"
                      };
      bools         = {};
      enums         = {"integrator"};
      strings       = {"name", "type"};
      multi_strings = {};
      integers      = {"samples"};
      reals         = {};
      positions     = {"center", "corner", "grid", "dr"};
      required      = {"name"};
      // I'd much rather see the default defined in simple native c++ as below
      // clang-format on

    }
    std::any assignAnyEnum(const std::string& name) const override;
  };
private:
  MagnetizationDensityInputSection input_section_;
  ///name of this Estimator
  std::string myName_;
  Integrator integrator_;
  Lattice cell_;
  PosType corner_;
  PosType dr_;
  PosType center_;
  TinyVector<int, DIM> grid_;
  int npoints_;
  int nsamples_; //Number of grid points for spin quadrature, or samples for Monte Carlo. 
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
#endif  /* QMCPLUSPLUS_MAGNETIZATION_DENSITY_INPUT_H */

