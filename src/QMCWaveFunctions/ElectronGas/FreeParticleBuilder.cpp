#include "FreeParticleBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/StructFact.h"
#include "FreeParticle.h"

namespace qmcplusplus
{

FreeParticleBuilder::FreeParticleBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("PW", comm),
    targetPtcl(els)
{}

std::unique_ptr<SPOSet> FreeParticleBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int npw, norb;
  PosType twist(0.0);
  OhmmsAttributeSet attrib;
  attrib.add(norb, "size");
  attrib.add(twist, "twist");
  attrib.put(cur);

#ifdef QMC_COMPLEX
  npw = norb;
#else
  npw = std::ceil((norb+1.0)/2);
#endif
  targetPtcl.setTwist(twist);
  PosType tvec = targetPtcl.getLattice().k_cart(twist);
  app_log() << "twist fraction = " << twist << std::endl;
  app_log() << "twist cartesian = " << tvec << std::endl;

  // extract npw k-points from container
  // kpts_cart is sorted by magnitude
  std::vector<PosType> kpts;
  kpts.resize(npw);
  kpts[0] = tvec;
  for (int ik=1;ik<npw;ik++)
  {
    kpts[ik] = tvec+targetPtcl.getSK().getKLists().kpts_cart[ik-1];
  }
  auto sposet = std::make_unique<FreeParticle>(kpts);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
