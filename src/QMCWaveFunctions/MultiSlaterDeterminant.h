//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANT_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANT_ORBITAL_H
#include <Configuration.h>
#include <QMCWaveFunctions/OrbitalBase.h>
#include <QMCWaveFunctions/Fermion/SlaterDet.h>

namespace qmcplusplus
  {

  /** @ingroup OrbitalComponent
   *  @brief An AntiSymmetric OrbitalBase composed of a linear combination of SlaterDeterminants.
   *
   *\f[
   *MS({\bf R}) = \sum_n c_n S_n({\bf R})
   *\f].
   *
   *The (S)ingle(P)article(O)rbitalSet template parameter is an
   *engine which fills in single-particle orbital terms.
   *
   \f[
   \frac{\nabla_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n \nabla_i D_n}
   {\sum_{n=1}^M c_n D_n}
   \f]
   \f[
   \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N(\nabla_i
   S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
   \f]
   The Laplacian
   \f[
   \frac{\nabla^2_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N
   (\nabla_i^2S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
   \f]
   */
  class MultiSlaterDeterminant: public OrbitalBase
    {

    public:

      typedef SlaterDet DeterminantSet_t;

      ///constructor
      MultiSlaterDeterminant();

      ///destructor
      ~MultiSlaterDeterminant();

      void checkInVariables(opt_variables_type& active);
      void checkOutVariables(const opt_variables_type& active);
      void resetParameters(const opt_variables_type& active);
      void reportStatus(ostream& os);

      void resetTargetParticleSet(ParticleSet& P);

      ValueType
      evaluate(ParticleSet& P
               ,ParticleSet::ParticleGradient_t& G
               ,ParticleSet::ParticleLaplacian_t& L);

      ValueType
      evaluateLog(ParticleSet& P //const DistanceTableData* dtable,
                  , ParticleSet::ParticleGradient_t& G
                  , ParticleSet::ParticleLaplacian_t& L);

      GradType evalGrad(ParticleSet& P, int iat);
      ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
      ValueType ratio(ParticleSet& P, int iat
                      , ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL);

      ValueType ratio(ParticleSet& P, int iat);
      void acceptMove(ParticleSet& P, int iat);
      void restore(int iat);

      void update(ParticleSet& P
                  , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL
                  , int iat);
      RealType evaluateLog(ParticleSet& P,BufferType& buf);
      RealType registerData(ParticleSet& P, BufferType& buf);
      RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);
      void copyFromBuffer(ParticleSet& P, BufferType& buf);

      OrbitalBasePtr makeClone(ParticleSet& tqp) const;
      void evaluateDerivatives(ParticleSet& P,
                               const opt_variables_type& optvars,
                               vector<RealType>& dlogpsi,
                               vector<RealType>& dhpsioverpsi);

      /**
        add a new SlaterDeterminant with coefficient c to the
        list of determinants
        */
      void add(DeterminantSet_t* sdet, RealType c);

      void add(DeterminantSet_t* sdet, RealType c, const string& id);

      vector<DeterminantSet_t*> SDets;
      vector<RealType> C;
      opt_variables_type myVars;
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
