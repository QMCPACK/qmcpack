//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    




#include "Numerics/OhmmsBlas.h"

const int BLAS::INCX;
const int BLAS::INCY;
const char BLAS::UPLO;
const char BLAS::TRANS;
const char BLAS::NOTRANS;

const float  BLAS::sone = 1.0e0;
const float  BLAS::szero = 0.0e0;
const double BLAS::done = 1.0e0;
const double BLAS::dzero = 0.0e0;
const std::complex<float> BLAS::cone = 1.0e0;
const std::complex<float> BLAS::czero = 0.0e0;
const std::complex<double> BLAS::zone = 1.0e0;
const std::complex<double> BLAS::zzero = 0.0e0;

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

