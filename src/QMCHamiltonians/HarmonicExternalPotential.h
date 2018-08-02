//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_HARMONIC_EXTERNAL_POTENTIAL_H
#define QMCPLUSPLUS_HARMONIC_EXTERNAL_POTENTIAL_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>


namespace qmcplusplus
{
  struct HarmonicExternalPotential : public QMCHamiltonianBase
  {
    //data members
    RealType mass;     
    RealType energy;
    RealType length;
    PosType  center;   
    const ParticleSet& Ps;

#if !defined(REMOVE_TRACEMANAGER)
    ///single particle trace sample array
    Array<TraceReal,1>* V_sample;
#endif
    
    //construction/destruction
    HarmonicExternalPotential(ParticleSet& P) : Ps(P)  
    {
      set_energy_domain(potential);
      one_body_quantum_domain(P);
    }

    ~HarmonicExternalPotential() { }

    //unneeded interface functions
    void resetTargetParticleSet(ParticleSet& P) { }

    //standard interface functions
    bool put(xmlNodePtr cur);                        
    bool get(std::ostream& os) const;
    QMCHamiltonianBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);

    //functions for physical (hamiltonian component) estimator
    Return_t evaluate(ParticleSet& P);
    inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
    {
      return evaluate(P);
    }

#if !defined(REMOVE_TRACEMANAGER)
    //traces interface
    virtual void contribute_particle_quantities()
    {
      request.contribute_array(myName);
    }

    virtual void checkout_particle_quantities(TraceManager& tm)
    {
      streaming_particles = request.streaming_array(myName);
      if( streaming_particles)
        V_sample = tm.checkout_real<1>(myName,Ps);
    }

    virtual void delete_particle_quantities()
    {
      if( streaming_particles)
        delete V_sample;
    }

    //  not really for interface, just collects traces
    inline Return_t evaluate_sp(ParticleSet& P);
#endif

  };
}
#endif
