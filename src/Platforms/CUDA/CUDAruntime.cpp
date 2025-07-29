//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "CUDAruntime.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

void cudaAssert(cudaError_t code, const std::string& cause, const char* filename, int line, bool abort)
{
  if (code != cudaSuccess)
  {
    std::ostringstream err;
    err << "cudaAssert: " << cudaGetErrorName(code) << " " << cudaGetErrorString(code) << ", file " << filename
        << ", line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    if (abort)
      throw std::runtime_error(cause);
  }
}
