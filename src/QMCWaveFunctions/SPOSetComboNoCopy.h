//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SPOSETCOMBO_NOCOPY_H
#define QMCPLUSPLUS_SPOSETCOMBO_NOCOPY_H

#include <QMCWaveFunctions/SPOSetBase.h>

namespace qmcplusplus
{

  /** a composite SPO without temporary containers
   *
   * Each SPO component knows where to put their data
   */
  struct SPOSetComboNoCopy: public SPOSetBase
  {
    std::vector<SPOSetBase*> mySPOs;

    /** default constructor */
    SPOSetComboNoCopy() { }

    inline SPOSetBase* makeClone() const
    {
      SPOSetComboNoCopy* aclone=new SPOSetComboNoCopy(*this);
      for(int k=0; k<mySPOs.size(); ++k)
        aclone->mySPOs[k]=mySPOs[k]->makeClone();
      aclone->setOrbitalSetSize(OrbitalSetSize);
      return aclone;
    }

    inline void add(SPOSetBase* a)
    {
      mySPOs.push_back(a);
    }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
    {
      for(int k=0; k<mySPOs.size(); ++k)
        mySPOs[k]->evaluate(P,iat,psi);
    }

    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      for(int k=0; k<mySPOs.size(); ++k)
        mySPOs[k]->evaluate(P,iat,psi,dpsi,d2psi);
    }
    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
    {
      for(int k=0; k<mySPOs.size(); ++k)
        mySPOs[k]->evaluate(P,iat,psi,dpsi,grad_grad_psi);
    }

    inline void resetParameters(const opt_variables_type& active)
    { }

    inline void resetTargetParticleSet(ParticleSet& e)
    { }

    inline void setOrbitalSetSize(int norbs)
    {
      OrbitalSetSize = norbs;
      BasisSetSize=norbs;
      for(int i=0; i<mySPOs.size(); ++i)
        mySPOs[i]->setOrbitalSetSize(norbs);
    }

    inline  void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      for(int k=0; k<mySPOs.size(); ++k)
        mySPOs[k]->evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet);
    }

    void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
    {
      for(int k=0; k<mySPOs.size(); ++k)
        mySPOs[k]->evaluate_notranspose(P,first,last,logdet,dlogdet,grad_grad_logdet);
    }

  };

}
#endif
