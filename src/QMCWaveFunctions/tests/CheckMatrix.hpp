//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CHECKMATRIX_HPP
#define QMCPLUSPLUS_CHECKMATRIX_HPP

#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus {
template<typename T1, typename ALLOC1, typename T2, typename ALLOC2>
void checkMatrix(Matrix<T1, ALLOC1>& a, Matrix<T2, ALLOC2>& b, const std::string & desc = "", int line = 0)
{
  REQUIRE(a.rows() >= b.rows());
  REQUIRE(a.cols() >= b.cols());
  auto matrixElementError = [line](int i, int j, auto& a, auto& b, const std::string& desc) -> std::string {
                      	          std::stringstream error_msg;
				  error_msg << "In " << desc << ":" << line <<  "\nbad element at " << i << ":" << j
					    <<"  " << a(i,j) << " != " << b(i,j) << '\n';
				  return error_msg.str();
			    };
  for (int i = 0; i < b.rows(); i++)
    for (int j = 0; j < b.cols(); j++)
    {
      CHECKED_ELSE(a(i, j) == ValueApprox(b(i, j))) {
	FAIL( matrixElementError(i,j,a,b,desc) );
	      }
    }
}
}

#endif
