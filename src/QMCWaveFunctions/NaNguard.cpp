//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <NaNguard.h>

namespace qmcplusplus
{

void NaNguard::checkOneParticleRatio(const PsiValue& ratio, const std::string_view info)
{
  if (qmcplusplus::isnan(std::norm(ratio)))
  {
    std::ostringstream error_message;
    error_message << "NaNguard::checkOneParticleRatio error message: " << info << std::endl
                  << "  ratio = " << ratio << std::endl;
    throw std::runtime_error(error_message.str());
  }
}

void NaNguard::checkOneParticleGradients(const GradType& grads, const std::string_view info)
{
  if (qmcplusplus::isnan(std::norm(dot(grads, grads))))
  {
    std::ostringstream error_message;
    error_message << "NaNguard::checkOneParticleGradients error message: " << info << std::endl;
    for (int i = 0; i < grads.size(); ++i)
      if (qmcplusplus::isnan(std::norm(grads[i])))
        error_message << "  grads[" << i << "] = " << grads[i] << std::endl;
    throw std::runtime_error(error_message.str());
  }
}

} // namespace qmcplusplus
