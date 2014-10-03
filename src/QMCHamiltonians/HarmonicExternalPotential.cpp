//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////


#include <QMCHamiltonians/HarmonicExternalPotential.h>
#include <OhmmsData/AttributeSet.h>


namespace qmcplusplus
{

  bool HarmonicExternalPotential::put(xmlNodePtr cur)                        
  { 
    using std::sqrt;

    mass      = -1.0;
    energy    = -1.0;
    length    = -1.0;
    center    =  0.0;

    OhmmsAttributeSet attrib;
    attrib.add(mass,"mass");
    attrib.add(energy,"frequency");
    attrib.add(energy,"energy");
    attrib.add(length,"length");
    attrib.add(center,"center");
    attrib.put(cur);

    if(energy<0.0)
      energy = 1.0;
    if(mass<0.0 && length<0.0)
      length = 1.0;
    if(mass<0.0)
      mass = 1.0/(energy*length*length);
    else if(length<0.0)
      length = 1.0/sqrt(mass*energy);

    return true; 
  }


  bool HarmonicExternalPotential::get(std::ostream& os) const
  {
    os << "External harmonic potential"<<endl;
    return true;
  } 


  QMCHamiltonianBase* HarmonicExternalPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
  {
    return new HarmonicExternalPotential(*this);
  }


  HarmonicExternalPotential::Return_t 
  HarmonicExternalPotential::evaluate(ParticleSet& P)
  {
    if(streaming_particles)
      Value = evaluate_sp(P);
    else
    {
      Value = 0.0;
      RealType prefactor = .5*energy/(length*length);
      for(int i=0;i<P.getTotalNum();++i)
      {
        PosType r = P.R[i]-center;
        Value += prefactor*dot(r,r);
      }
    }
    return Value;
  }


  HarmonicExternalPotential::Return_t 
  HarmonicExternalPotential::evaluate_sp(ParticleSet& P)
  {
    Array<TraceReal,1>& V_samp = *V_sample;
    Value = 0.0;
    RealType prefactor = .5*energy/(length*length);
    for(int i=0;i<P.getTotalNum();++i)
    {
      PosType r = P.R[i]-center;
      RealType v1 = prefactor*dot(r,r);
      V_samp(i) = v1;
      Value += v1;
    }
    return Value;
  }

}
