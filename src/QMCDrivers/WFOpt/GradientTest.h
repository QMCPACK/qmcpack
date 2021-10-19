//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_GRADIENTTEST_H
#define QMCPLUSPLUS_GRADIENTTEST_H

#include "GradientTestInput.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include <fstream>

/**
 * @file GradientTest.h
 * Tests for variational parameter derivatives
 *
 * This optimization type will compare finite difference derivatives with the analytic derivatives.
 * It can also output a file that can be used with 'qmca' to get error bars.
 *
 * The input is placed under the batched linear optimizer (method="opt_batch").
 * It uses the 'optimize' tag with the 'method' attribute set to 'gradient_test'
 * For example:
 * \code
 *  <loop max="4">
 *   <qmc method="opt_batch" move="pbyp">
 *     <optimize method="gradient_test">
 *       <parameter name="output_param_file">yes</parameter>
 *     </optimize>
 *     ...
 *   </qmc>
 *  </loop>
 * \endcode
 *
 * The output_param_file parameter defaults to 'no'.
 * If 'yes', a file named "<project id>.param.s000.scalar.dat" is created.  Each iteration of the optimizer
 * loop outputs one line in the file. (The above example will produce 4 entries in the file).
 * This is a hack to enable computing error bars on the parameter gradients.
 *
 */


namespace qmcplusplus
{


class GradientTest
{
public:
  using Return_t = QMCTraits::RealType;

  GradientTest(GradientTestInput&& input) : input_(std::move(input)), first_(true), param_deriv_index_(0) {}

  // Compute and compare the numerical and analytic gradients
  void run(QMCCostFunctionBase& costFunc, const std::string& root_name);

private:
  GradientTestInput input_;

  std::ofstream param_deriv_file_;
  bool first_;
  int param_deriv_index_;
};

} // namespace qmcplusplus

#endif
