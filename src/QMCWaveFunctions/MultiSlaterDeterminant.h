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

namespace qmcplusplus 
{

  class SlaterDet;

  /** @ingroup OrbitalComponent
   *  @brief An AntiSymmetric OrbitalBase composed of a linear combination of SlaterDeterminants. 
   *
   *\f$ MS({\bf R}) = \sum_n c_n S_n({\bf R})\f$.
   *
   *The (S)ingle(P)article(O)rbitalSet template parameter is an 
   *engine which fills in single-particle orbital terms.  
   * 
   \f[
   \frac{\nabla_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n \nabla_i D_n} {\sum_{n=1}^M c_n D_n}
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

      ///constructor
      MultiSlaterDeterminant(const ParticleSet& p);

      ///destructor
      ~MultiSlaterDeterminant();

      /** add a new SlaterDeterminant 
       * @param sdet Slater-Determinant
       * @param c coefficient
       * @param id name
       */
      void add(SlaterDet* sdet, RealType c, const string& id) ;

      void checkInVariables(opt_variables_type& active);
      void checkOutVariables(const opt_variables_type& active);
      void resetParameters(const opt_variables_type& active);
      void reportStatus(ostream& os);

      void resetTargetParticleSet(ParticleSet& P);

      ValueType evaluate(ParticleSet& P
            ,ParticleSet::ParticleGradient_t& G
            .ParticleSet::ParticleLaplacian_t& L);

      ValueType evaluateLog(ParticleSet& P
            ,ParticleSet::ParticleGradient_t& G
            .ParticleSet::ParticleLaplacian_t& L);

      ValueType ratio(ParticleSet& P, int iat 
          ,ParticleSet::ParticleGradient_t& dG 
          ,ParticleSet::ParticleLaplacian_t& dL);

      void acceptMove(ParticleSet& P, int iat);

      void restore(int iat);

      ValueType ratio(ParticleSet& P, int iat);

      void update(ParticleSet& P
          ,ParticleSet::ParticleGradient_t& dG
          ,ParticleSet::ParticleLaplacian_t& dL
          ,int iat);

      RealType evaluateLog(ParticleSet& P,BufferType& buf);
      RealType registerData(ParticleSet& P, BufferType& buf);
      RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);

      OrbitalBasePtr makeClone(ParticleSet& tqp) const;

    private:
      /** CI sum
       *
       * my_psi \f$ = \sum_i my\_c[i] my\_sdet[i]({\bf R})\f$.
       */
      ValueType my_psi;
      /** my_psiinv=1/my_psi
       */
      ValueType my_psiinv;
      /** cur_sdet[i] = value of the i-th SlaterDetermiant */
      vector<ValueType> cur_sdet;
      /** cur_ratio[i]= ratio of the i-th SlaterDetermiant */
      vector<ValueType> cur_ratio;
      /** temporary gradients */
      ParticleSet::ParticleGradient_t  my_grad, t_grad;
      /** temporary laplacians */
      ParticleSet::ParticleLaplacian_t my_lap, t_lap;
      /** SlaterDeterminant components */
      vector<SlaterDet_t*> my_sdet;
      /** Coefficients of the SlaterDeterminant components */
      vector<RealType> my_c;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
