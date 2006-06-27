//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_MATRIXOPERATORS_H
#define QMCPLUSPLUS_MATRIXOPERATORS_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Numerics/OhmmsBlas.h"
namespace qmcplusplus {
  struct MatrixOperators {
    /** static function to perform C=AB for real matrices
     *
     * Call dgemm
     */
    inline static void product(const Matrix<double>& A,
        const Matrix<double>& B, Matrix<double>& C) {
      const char transa = 'N';
      const char transb = 'N';
      const double one=1.0;
      const double zero=0.0;
      dgemm(transa, transb, B.cols(), A.rows(), B.rows(),
          one, B.data(), B.cols(), A.data(), A.cols(),
          zero, C.data(), C.cols());
    }

    /** static function to perform C=AB for complex matrices
     *
     * Call zgemm
     */
    inline static void product(const Matrix<std::complex<double> >& A,
        const Matrix<std::complex<double> >& B, Matrix<std::complex<double> >& C) {
      const char transa = 'N';
      const char transb = 'N';
      const std::complex<double> zone(1.0,0.0);
      const std::complex<double> zero(0.0,0.0);
      zgemm(transa, transb, B.cols(), A.rows(), B.rows(),
          zone, B.data(), B.cols(), A.data(), A.cols(),
          zero, C.data(), C.cols());
    }

    /** static function to perform y=Ax for generic matrix and vector
     */
    inline static void product(const Matrix<double>& A, const Vector<double>& x, double* restrict yptr) {
      const char transa = 'T';
      const double one=1.0;
      const double zero=0.0;
      dgemv(transa, A.cols(), A.rows(), one, A.data(), A.cols(), x.data(), 1, zero, yptr, 1);
    }

    /** static function to perform y=Ax for generic matrix and vector
     */
    inline static void product(const Matrix<std::complex<double> >& A, 
        const Vector<std::complex<double> >& x, 
        std::complex<double>* restrict yptr) {
      const char transa = 'T';
      const std::complex<double> zone(1.0,0.0);
      const std::complex<double> zero(0.0,0.0);
      zgemv(transa, A.cols(), A.rows(), zone, A.data(), A.cols(), x.data(), 1, zero, yptr, 1);
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
