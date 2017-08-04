//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_JELLIUMPOTENTIAL_H
#define QMCPLUSPLUS_JELLIUMPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include <numeric>

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *@brief CoulombPotential for the different source and target particle sets.
 *
 * \f[ H = \sum_i \frac{Z(i)q}{r} \f]
 * where \f$ Z(i) \f$ is the effective charge of the Ith
 * ion and \f$ q \f$ is the charge of the set of quantum particles.
 * For instance, \f$ q = -1 \f$ for electrons and
 * \f$ q = 1 \f$ for positrons.
 *
 * @warning need to be generalized by checking visitor.Species.
 */
struct JelliumPotential: public QMCHamiltonianBase
{

  ///number of ions
  int Centers;
  ParticleSet& sourcePtcl;
  DistanceTableData* d_table;
  ///container for the ion charges
  std::vector<RealType> Z,RC,RS;

  JelliumPotential(ParticleSet& ions, ParticleSet& els):
    sourcePtcl(ions), d_table(0)
  {
    d_table = DistanceTable::add(ions,els,DT_AOS);
    //index for attribute charge
    SpeciesSet& tspecies(ions.getSpeciesSet());
    int iz = tspecies.addAttribute("charge");
    int rs = tspecies.addAttribute("rs");
    Centers = ions.getTotalNum();
    RC.resize(Centers);
    RS.resize(Centers);
    Z.resize(Centers);
    RealType C = -1.0;
    for(int iat=0; iat<Centers; iat++)
    {
      Z[iat] = tspecies( iz ,ions.GroupID[iat])*C;
      RS[iat] = tspecies(rs,ions.GroupID[iat]);
      RC[iat] = std::pow(-Z[iat],1.0/3.0)*RS[iat];
      RS[iat] = 1.0/(RS[iat]*RS[iat]*RS[iat]);
      app_log()<<" rs is "<<rs<< std::endl;
      app_log()<<" iz is "<<iz<< std::endl;
      app_log()<<" RC is "<<RC[iat]<< std::endl;
      app_log()<<" RS^-3 is "<<RS[iat]<< std::endl;
      app_log()<<" Z is "<<Z[iat]<< std::endl;
    }
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(sourcePtcl,P,DT_AOS);
  }

  ~JelliumPotential() { }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value=0.0;
    for(int iat=0; iat<Centers; ++iat)
    {
      for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; ++nn)
        if (d_table->r(nn)<RC[iat])
          Value+= RS[iat]*std::pow(d_table->r(nn),2);
        else
          Value+= Z[iat]*d_table->rinv(nn);
    }
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t
  registerData(ParticleSet& P, BufferType& buffer)
  {
    NewValue=evaluate(P);
    buffer.add(Value);
    return Value;
  }

  inline Return_t
  updateBuffer(ParticleSet& P, BufferType& buffer)
  {
    NewValue=evaluate(P);
    buffer.put(Value);
    return Value;
  }

  inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.get(Value);
  }

  inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.put(Value);
  }

  inline Return_t
  evaluatePbyP(ParticleSet& P, int active)
  {
    APP_ABORT("JelliumPotential::evaluatePbyP should not be used");
    return 0.0;
    //const std::vector<DistanceTableData::TempDistType> &temp(d_table->Temp);
    //Return_t del=0.0;
    //for(int iat=0; iat<Centers; ++iat)
    //{
    //  if (temp[iat].r1<RC[iat]) del+=RS[iat]*std::pow(temp[iat].r1,2);
    //  else del +=Z[iat]*temp[iat].rinv1;
    //  if (temp[iat].r0<RC[iat]) del-=RS[iat]*std::pow(temp[iat].r0,2);
    //  else del -=Z[iat]*temp[iat].rinv0;
    //}
    //
    //return NewValue=Value+del;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "JelliumPotential potential: " << sourcePtcl.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new JelliumPotential(sourcePtcl, qp);
  }

};
}
#endif


