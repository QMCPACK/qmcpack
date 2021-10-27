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
  std::vector<Return_t> a_xi(num_params);
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
  double finite_diff_delta = 1e-5;
  costFunc.GradCost(numeric_grad, params, finite_diff_delta);

  // Analytic gradient
  costFunc.GradCost(analytic_grad, params);


  if (input_.do_param_output())
  {
    param_deriv_file_ << param_deriv_index_ << " ";
    param_deriv_index_++;
  }

  app_log() << "Param_Name  Value    Numeric            Analytic       Percent" << std::endl;
  for (int k = 0; k < num_params; k++)
  {
    std::string vname = costFunc.getParamName(k);
    if (numeric_grad[k] != 0)
      app_log() << vname << " " << params[k] << "  " << numeric_grad[k] << "  " << analytic_grad[k] << "  "
                << 100 * (numeric_grad[k] - analytic_grad[k]) / numeric_grad[k] << std::endl;
    else
      app_log() << vname << " " << params[k] << "  " << numeric_grad[k] << "  " << analytic_grad[k] << "   inf" << std::endl;
    if (input_.do_param_output())
      param_deriv_file_ << std::setprecision(10) << analytic_grad[k] << " ";
  }

  if (input_.do_param_output())
    param_deriv_file_ << std::endl;

  app_log() << std::endl;
}

} // namespace qmcplusplus
