#ifndef OHMMS_QMC_QDWFBUILDER_H
#define OHMMS_QMC_QDWFBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "QMCWaveFunctions/QDwf.h"

namespace ohmmsqmc {

  class QDwfBuilder: public OrbitalBuilderBase {

    typedef SingleParticleOrbitalSet<QDwf> SPOSet_t;

  public :

    QDwfBuilder(TrialWaveFunction& a) : OrbitalBuilderBase(a){}

    bool put(xmlNodePtr cur);

  };

}
#endif
