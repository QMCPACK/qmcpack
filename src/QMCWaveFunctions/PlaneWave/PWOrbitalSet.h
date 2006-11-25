//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
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
/** @file PWOrbitalSet.h
 * @brief Definition of member functions of Plane-wave basis set
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus {

  class PWOrbitalSet: public SPOSetBase {

    public:

      typedef PWBasis                    BasisSet_t;
#if defined(ENABLE_SMARTPOINTER)
      typedef boost::shared_ptr<PWBasis> PWBasisPtr;
#else
      typedef PWBasis*                   PWBasisPtr;
#endif

      /** inherit the enum of BasisSet_t */
      enum {PW_VALUE=   BasisSet_t::PW_VALUE, 
        PW_LAP=     BasisSet_t::PW_LAP, 
        PW_GRADX=   BasisSet_t::PW_GRADX, 
        PW_GRADY=   BasisSet_t::PW_GRADY, 
        PW_GRADZ=   BasisSet_t::PW_GRADZ,
        PW_MAXINDEX=BasisSet_t::PW_MAXINDEX
      };

      /** default constructor
      */
      PWOrbitalSet(): OwnBasisSet(false) {
      }

      /** delete BasisSet only it owns this
       *
       * Builder takes care of who owns what
       */
      ~PWOrbitalSet();

      void resize(PWBasisPtr bset, int nbands, bool cleanup=false) {
        myBasisSet=bset;
        OrbitalSetSize=nbands;
        OwnBasisSet=cleanup;
        BasisSetSize=myBasisSet->NumPlaneWaves;
        C.resize(OrbitalSetSize,BasisSetSize);
        Temp.resize(OrbitalSetSize,PW_MAXINDEX);
      }


      /** Builder class takes care of the assertion
      */
      inline void addVector(const std::vector<ValueType>& coefs,int jorb) {
        std::copy(coefs.begin(),coefs.end(),C[jorb]);
      }

      void reset();

      void setOrbitalSetSize(int norbs);

      void resetTargetParticleSet(ParticleSet& P);

      //inline RealType
      //  evaluate(const ParticleSet& P, int iat, int jorb) {
      //  LOGMSG("PWOSet: this should not be used");
      //  OHMMS::Controller->abort();
      //  return 0.0;
      //}

      void 
        evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

      void 
        evaluate(const ParticleSet& P, int iat, 
            ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

      void evaluate(const ParticleSet& P, int first, int last,
          ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

      /** boolean
       *
       * If true, this has to delete the BasisSet
       */
      bool OwnBasisSet;
      ///TwistAngle of this PWOrbitalSet
      PosType TwistAngle;
      ///My basis set
      PWBasisPtr myBasisSet;
      /////Plane-wave coefficients: (iband,g-vector)
      //Matrix<ValueType> Coefs;
      /** temporary array to perform gemm operation */
      Matrix<ValueType> Temp;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
