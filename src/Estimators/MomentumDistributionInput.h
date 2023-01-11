//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: MomentumEstimator.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_MOMENTUMDISTRIBUTIONINPUT_H
#define QMCPLUSPLUS_MOMENTUMDISTRIBUTIONINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

class MomentumDistribution;

/** Native representation for Momentum Distribution Estimators inputs
 */
class MomentumDistributionInput
{
public:
  using Consumer = MomentumDistribution;
  using Real = QMCTraits::RealType;

  class MomentumDistributionInputSection : public InputSection
  {
  public:
    // clang-format: off
    MomentumDistributionInputSection()
    {
      section_name   = "MomentumDistribution";
      attributes     = {"type", "name", "samples", "kmax", "kmax0", "kmax1", "kmax2"};
      strings        = {"type", "name"};
      integers       = {"samples"};
      reals          = {"kmax", "kmax0", "kmax1", "kmax2"};
      // default_values = {{"name", std::string("nofk")}, {"samples", int(40)}, {"kmax", Real(0.0)},
      //                   {"kmax0", Real(0.0)},          {"kmax1", Real(0.0)}, {"kmax2", Real(0.0)}};
    }
    // clang-format: on
  };

  MomentumDistributionInput(xmlNodePtr cur);
  /** default copy constructor
   *  This is required due to MDI being part of a variant used as a vector element.
   */
  MomentumDistributionInput(const MomentumDistributionInput&) = default;
private:
  MomentumDistributionInputSection input_section_;

  std::string name_{"nofk"};
  std::string type_;
  ///number of samples
  int samples_ = 40;
  //maximum k-value in the k-grid in cartesian coordinates
  Real kmax_ = 0.0;
  //maximum k-values in the k-grid along the reciprocal cell axis
  Real kmax0_ = 0.0;
  Real kmax1_ = 0.0;
  Real kmax2_ = 0.0;

public:
  const std::string& get_name() const { return name_; }
  const std::string& get_type() const { return type_; }
  const int& get_samples() const { return samples_; }
  const Real& get_kmax() const { return kmax_; }
  const Real& get_kmax0() const { return kmax0_; }
  const Real& get_kmax1() const { return kmax1_; }
  const Real& get_kmax2() const { return kmax2_; }  
};

} // namespace qmcplusplus
#endif /* MOMENTUMDISTRIBUTIONINPUT_H */
