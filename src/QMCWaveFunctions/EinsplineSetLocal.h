//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EINSPLINE_SET_LOCAL_H
#define QMCPLUSPLUS_EINSPLINE_SET_LOCAL_H

#include "QMCWaveFunctions/EinsplineOrb.h"
#include "QMCWaveFunctions/EinsplineSet.h"

namespace qmcplusplus {

  class EinsplineSetLocal : public EinsplineSet
  {
    friend class EinsplineSetBuilder;
  protected:
    /////////////////////
    // Orbital storage //
    /////////////////////
    /// Store the orbital objects.  Using template class allows us to
    /// avoid making separate real and complex versions of this class.
    //std::vector<EinsplineOrb<ValueType,OHMMS_DIM>*> Orbitals;
    std::vector<EinsplineOrb<complex<double>,OHMMS_DIM>*> Orbitals;

  public:
    SPOSetBase* makeClone() const;
    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
    void evaluate(const ParticleSet& P, int iat, 
		  ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& psi, GradMatrix_t& dpsi, 
		  ValueMatrix_t& d2psi);

    void resetParameters(const opt_variables_type& active);

    EinsplineSetLocal() 
    {
      className = "EinsplineSetLocal";
    }
  };
}
#endif
