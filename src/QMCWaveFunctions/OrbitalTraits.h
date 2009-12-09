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
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>

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

  template<>
    struct OrbitalTraits<float> {
      typedef float          real_type;
      typedef float          value_type;
      typedef std::complex<float> complex_type;
    };

  template<>
    struct OrbitalTraits<std::complex<float> > {
      typedef float          real_type;
      typedef std::complex<float> value_type;
      typedef std::complex<float> complex_type;
    };

  /** generic conversion from type T1 to type T2 using implicit conversion
  */
  template<typename T1, typename T2>
    inline void convert(const T1& in, T2& out)
    {
      out=in;
    }

  /** specialization of conversion from complex to real
  */
  template<typename T1>
    inline void convert(const std::complex<T1>& in, double out)
    {
      out=in.real();
    }

  template<typename T1>
    inline void convert(const std::complex<T1>& in, float out)
    {
      out=in.real();
    }

  /* specialization of D-dim vectors
   *
   */
  template<typename T1, typename T2, unsigned D>
    inline void convert(const TinyVector<T1,D>& in, TinyVector<T2,D>& out)
    {
      for(int i=0; i<D;++i) convert(in[i],out[i]);
    }

  /** specialization for 3D */
  template<typename T1, typename T2>
    inline void convert(const TinyVector<T1,3>& in, TinyVector<T2,3>& out)
    {
      convert(in[0],out[0]);
      convert(in[1],out[1]);
      convert(in[2],out[2]);
    }

  /** generic function to convert arrays
   * @param in starting address of type T1
   * @param out starting address of type T2
   * @param n size of in/out
   */
  template<typename T1, typename T2>
    inline void convert(const T1* restrict in, T2* restrict out, std::size_t n)
    {
      for(int i=0; i<n;++i) convert(in[i],out[i]);
    }

  /** specialization for a vector */
  template<typename T1, typename T2>
    inline void convert(const Vector<T1>& in, Vector<T2>& out)
    {
      convert(in.data(),out.data(),in.size());
    }

  /** specialization for a vector */
  template<typename T1, typename T2>
    inline void convert(const Matrix<T1>& in, Matrix<T2>& out)
    {
      convert(in.data(),out.data(),in.size());
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
