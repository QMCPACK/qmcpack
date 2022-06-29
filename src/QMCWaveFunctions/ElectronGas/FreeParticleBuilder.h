#ifndef QMCPLUSPLUS_FREE_PARTICLE_BUILDER_H
#define QMCPLUSPLUS_FREE_PARTICLE_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{
class FreeParticleBuilder : public SPOSetBuilder
{
public:
  FreeParticleBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur);
  ~FreeParticleBuilder(){}

  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;
private:
  ParticleSet& targetPtcl;
  bool in_list(const int j, const std::vector<int> l);
};
} // qmcplusplus
#endif
