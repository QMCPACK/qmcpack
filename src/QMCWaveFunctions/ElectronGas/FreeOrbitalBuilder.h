#ifndef QMCPLUSPLUS_FREE_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_FREE_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{
class FreeOrbitalBuilder : public SPOSetBuilder
{
public:
  FreeOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur);
  ~FreeOrbitalBuilder() {}

  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;

private:
  ParticleSet& targetPtcl;
  bool in_list(const int j, const std::vector<int> l);
};
} // namespace qmcplusplus
#endif
