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
/** @file BasisSetBase.h
 * @brief Declaration of a base class of BasisSet
 */
#ifndef QMCPLUSPLUS_ORBITALSETTRAITS_H
#define QMCPLUSPLUS_ORBITALSETTRAITS_H

#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Optimize/VarList.h"

namespace qmcplusplus {
  /** dummy class for templated classes
   */
  struct DummyGrid {
    inline void locate(double r) {}
  };

  typedef TinyVector<int,4> QuantumNumberType;

  enum {q_n=0,q_l,q_m, q_s};

  /** trait class to handel a set of Orbitals
   */
  template<typename T>
  struct OrbitalSetTraits: public OrbitalTraits<T> {
    enum {DIM=OHMMS_DIM};
    typedef typename OrbitalTraits<T>::real_type RealType;
    typedef typename OrbitalTraits<T>::value_type ValueType;
    typedef int                            IndexType;
    typedef TinyVector<RealType,DIM>       PosType;
    typedef TinyVector<ValueType,DIM>      GradType;
    typedef Tensor<RealType,DIM>           TensorType;
    typedef Vector<IndexType>     IndexVector_t;
    typedef Vector<ValueType>     ValueVector_t;
    typedef Matrix<ValueType>     ValueMatrix_t;
    typedef Vector<GradType>      GradVector_t;
    typedef Matrix<GradType>      GradMatrix_t;
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
