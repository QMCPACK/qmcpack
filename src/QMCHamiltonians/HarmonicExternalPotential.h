//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

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
    ParticleSet* Ps;
    ///single particle trace sample array
    Array<TraceReal,1>* V_sample;
    
    //construction/destruction
    HarmonicExternalPotential(ParticleSet& P) : Ps(&P)  { }
    ~HarmonicExternalPotential() { }

    //unneeded interface functions
    void resetTargetParticleSet(ParticleSet& P) { }

    //standard interface functions
    bool put(xmlNodePtr cur);                        
    bool get(std::ostream& os) const;
    QMCHamiltonianBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);

    //functions for physical (hamiltonian component) estimator
    Return_t evaluate(ParticleSet& P);
    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
    {
      return evaluate(P);
    }

    //traces interface
    virtual void checkout_particle_arrays(TraceManager& tm)
    {
      V_sample = tm.checkout_real<1>(myName,*Ps);
    }
    virtual void delete_particle_arrays()
    {
      delete V_sample;
    }
    //  not really for interface, just implements traces
    inline Return_t spevaluate(ParticleSet& P);
  };
}
#endif
