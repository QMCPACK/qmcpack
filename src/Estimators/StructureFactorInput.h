//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/{SkEstimator.h, SkAllEstimator.h}
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_STRUCTUREFACTORINPUT_H
#define QMCPLUSPLUS_STRUCTUREFACTORINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

class StructureFactorEstimator;

/** Native representation for SK Estimator inputs
 *
 *  This is used for both "structure factor estimators" which are merged in the port to
 *  the unified driver.
 */
class StructureFactorInput
{
public:
  using Consumer = StructureFactorEstimator;
  using Real     = QMCTraits::FullPrecRealType;

  class StructureFactorInputSection : public InputSection
  {
  public:
    // clang-format: off
    StructureFactorInputSection()
    {
      section_name = "StructureFactor";
      attributes   = {"type", "name", "source", "target", "hdf5", "writerho", "writeionion"};
      strings      = {"type", "name", "source", "target"};
      bools        = {"hdf5", "writerho", "writeionion"};
    }
    // clang-format: on
    void checkParticularValidity() override;
  };

  StructureFactorInput(xmlNodePtr cur);
  /** default copy constructor
   *  This is required due to MDI being part of a variant used as a vector element.
   */
  StructureFactorInput(const StructureFactorInput&) = default;

private:
  StructureFactorInputSection input_section_;

  std::string name_{"sk"};
  std::string type_{"sk"};
  std::string source_;
  std::string target_;
  bool write_hdf5_{false};
  bool write_rho_{false};
  bool write_ion_ion_{false};

public:
  std::string get_name() const { return name_; }
  std::string get_type() const { return type_; }
  std::string get_source() const { return source_; }
  std::string get_target() const { return target_; }
  bool get_write_hdf5() const { return write_hdf5_; }
  bool get_write_rho() const { return write_rho_; }
  bool get_write_ion_ion() const { return write_ion_ion_; }
};

} // namespace qmcplusplus
#endif /* QMCPLUSPLUS_STRUCTUREFACTORINPUT_H */
