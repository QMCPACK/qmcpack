#ifndef QMCPLUSPLUS_FREE_ORBITAL_BUILDERT_H
#define QMCPLUSPLUS_FREE_ORBITAL_BUILDERT_H

#include "QMCWaveFunctions/SPOSetBuilderT.h"

namespace qmcplusplus
{
template <typename T>
class FreeOrbitalBuilderT : public SPOSetBuilderT<T>
{
public:
    using RealType = typename SPOSetBuilderT<T>::RealType;
    using PosType = typename SPOSetBuilderT<T>::PosType;

  FreeOrbitalBuilderT(ParticleSetT<T>& els, Communicate* comm, xmlNodePtr cur);
  ~FreeOrbitalBuilderT() {}

  std::unique_ptr<SPOSetT<T>> createSPOSetFromXML(xmlNodePtr cur) override;

private:
  ParticleSetT<T>& targetPtcl;
  bool in_list(const int j, const std::vector<int> l);
};
} // namespace qmcplusplus
#endif
