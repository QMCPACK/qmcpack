//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ORBITALTRAITS_H
#define QMCPLUSPLUS_ORBITALTRAITS_H
#include <complex>
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus {

  template<class T> struct OrbitalTraits {};

  template<>
    struct OrbitalTraits<double> {
      typedef double          real_type;
      typedef double          value_type;
      typedef std::complex<double> complex_type;
    };

  template<>
    struct OrbitalTraits<std::complex<double> > {
      typedef double          real_type;
      typedef std::complex<double> value_type;
      typedef std::complex<double> complex_type;
    };

  inline double real(double a) {
    return a;
  }
  
  inline TinyVector<double,3> real(const TinyVector<double,3>& a) { 
    return a;
  }

  inline TinyVector<double,3> real(const TinyVector<std::complex<double>,3>& a) { 
    return TinyVector<double,3>(a[0].real(),a[1].real(),a[2].real());
  }

  inline TinyVector<double,2> real(const TinyVector<double,2>& a) { 
    return a;
  }

  inline TinyVector<double,2> real(const TinyVector<std::complex<double>,2>& a) { 
    return TinyVector<double,2>(a[0].real(),a[1].real());
  }

  inline TinyVector<double,1> real(const TinyVector<double,1>& a) { 
    return a;
  }

  inline TinyVector<double,1> real(const TinyVector<std::complex<double>,1>& a) { 
    return TinyVector<double,1>(a[0].real());
  }

  template<typename T> inline double evaluatePhase(T sign_v)
  {
    return (sign_v>0)?0.0:M_PI;
  }

  template<>
  inline double evaluatePhase(const std::complex<double>& psi)
  {
    return std::arg(psi);
  }

  /** evaluate the log(|psi|) and phase
   * @param psi real/complex value
   * @param phase phase of psi
   * @return log(|psi|)
   */
  template<class T>
    inline T evaluateLogAndPhase(const T psi, T& phase) {
      if(psi<0.0) {
        phase= M_PI;
        return std::log(-psi);
      } else {
        phase = 0.0;
        return std::log(psi);
      }
    }

  template<class T>
    inline T
    evaluateLogAndPhase(const std::complex<T>& psi, T& phase) {
      phase = std::arg(psi);
      if(phase<0.0) phase += 2.0*M_PI;
      return 0.5*std::log(psi.real()*psi.real()+psi.imag()*psi.imag());
      //return std::log(psi);
    }

    inline double evaluatePhase(const double psi) {
      return (psi<numeric_limits<double>::epsilon())?M_PI:0.0;
    }

    inline double evaluatePhase(const std::complex<double>& psi) {
      return std::arg(psi);
    }
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
