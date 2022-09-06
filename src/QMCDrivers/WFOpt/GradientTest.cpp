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


#include "GradientTest.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"

namespace qmcplusplus
{

void GradientTest::run(QMCCostFunctionBase& costFunc, const std::string& root_name)
{
  int num_params = costFunc.getNumParams();
  std::vector<Return_t> numeric_grad(num_params);
  std::vector<Return_t> params(num_params);
  std::vector<Return_t> analytic_grad(num_params);

  for (int i = 0; i < num_params; i++)
  {
    params[i] = std::real(costFunc.Params(i));
  }


  if (input_.do_param_output() && first_)
  {
    int namelen = root_name.length();
    // Assume that the root_name has the series suffix (".s000").
    // Remove the series suffix to get the project id
    std::string fname = root_name.substr(0, namelen - 5) + ".param.s000.scalar.dat";

    param_deriv_file_.open(fname);
    param_deriv_file_ << "# Index ";
    for (int i = 0; i < num_params; i++)
    {
      param_deriv_file_ << " " << costFunc.getParamName(i);
    }
    param_deriv_file_ << std::endl;
    first_ = false;
  }

  // Numerical gradient
  costFunc.GradCost(numeric_grad, params, input_.get_finite_diff_delta());

  // Analytic gradient
  costFunc.GradCost(analytic_grad, params);


  if (input_.do_param_output())
  {
    param_deriv_file_ << param_deriv_index_ << " ";
    param_deriv_index_++;
  }

  std::string_view param_name_header("Param_Name");
  size_t max_name_len = param_name_header.size();
  for (int k = 0; k < num_params; k++)
  {
    std::string vname = costFunc.getParamName(k);
    max_name_len      = std::max(vname.size(), max_name_len);
  }
  max_name_len += 2; // add some padding

  // clang-format off
  app_log() << std::setw(max_name_len) << std::left << param_name_header
            << std::setw(14) << std::right << " Value "
            << std::setw(20) << std::right << " Numeric "
            << std::setw(20) << std::right << " Analytic "
            << std::setw(14) << " Percent" << std::endl;
  // clang-format on

  for (int k = 0; k < num_params; k++)
  {
    std::string vname = costFunc.getParamName(k);
    std::ostringstream rel_diff_str;
    std::string over_threshold;
    if (numeric_grad[k] != 0)
    {
      double rel_diff_percent = 100 * (numeric_grad[k] - analytic_grad[k]) / numeric_grad[k];
      rel_diff_str << std::scientific << std::setprecision(2) << rel_diff_percent;
      // Highlight problematic differences
      // The thresholds are arbitrary.
      if (std::abs(rel_diff_percent) > 1e-3)
        over_threshold = " !";
      if (std::abs(rel_diff_percent) > 1e-2)
        over_threshold = " !!";
      if (std::abs(rel_diff_percent) > 1e-1)
        over_threshold = " !!!";
    }
    else
    {
      rel_diff_str << "inf";
      over_threshold = " !!!";
    }


    // clang-format off
    app_log() << std::setw(max_name_len) << std::left << vname
              << std::setprecision(6)  << std::setw(14) << std::right << params[k]
              << std::setprecision(10) << std::setw(20) << std::right << numeric_grad[k]
              << std::setw(20) << std::right << analytic_grad[k]
              << std::setw(14) << std::right << rel_diff_str.str() << over_threshold << std::endl;
    // clang-format on

    if (input_.do_param_output())
      param_deriv_file_ << std::setprecision(10) << analytic_grad[k] << " ";
  }

  // Reset precision to the default
  app_log() << std::setprecision(6);

  if (input_.do_param_output())
    param_deriv_file_ << std::endl;

  app_log() << std::endl;
}

} // namespace qmcplusplus
