#include "OhmmsData/AttributeSet.h"
#include "LongRange/StructFact.h"
#include "LongRange/KContainer.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbital.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbitalBuilder.h"

namespace qmcplusplus
{
FreeOrbitalBuilder::FreeOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("PW", comm), targetPtcl(els)
{}

std::unique_ptr<SPOSet> FreeOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int npw, norb;
  PosType twist(0.0);
  OhmmsAttributeSet attrib;
  attrib.add(norb, "size");
  attrib.add(twist, "twist");
  attrib.put(cur);

  auto lattice = targetPtcl.getLattice();

  PosType tvec = lattice.k_cart(twist);
#ifdef QMC_COMPLEX
  npw = norb;
  targetPtcl.setTwist(twist);
  app_log() << "twist fraction = " << twist << std::endl;
  app_log() << "twist cartesian = " << tvec << std::endl;
#else
  npw = std::ceil((norb + 1.0) / 2);
  for (int ldim = 0; ldim < twist.size(); ldim++)
  {
    if (std::abs(twist[ldim]) > 1e-16)
      throw std::runtime_error("no twist for real orbitals");
  }
#endif

  // extract npw k-points from container
  // kpts_cart is sorted by magnitude
  std::vector<PosType> kpts(npw);
  KContainer klists;
  RealType kcut = lattice.LR_kc; // to-do: reduce kcut to >~ kf
  klists.updateKLists(lattice, kcut, lattice.ndim, twist);

  // k0 is not in kpts_cart
  kpts[0] = tvec;
#ifdef QMC_COMPLEX
  for (int ik = 1; ik < npw; ik++)
  {
    kpts[ik] = klists.kpts_cart[ik - 1];
  }
#else
  const int nktot = klists.kpts.size();
  std::vector<int> mkidx(npw, 0);
  int ik = 1;
  for (int jk = 0; jk < nktot; jk++)
  {
    // check if -k is already chosen
    const int jmk = klists.minusk[jk];
    if (in_list(jk, mkidx))
      continue;
    // if not, then add this kpoint
    kpts[ik]  = klists.kpts_cart[jk];
    mkidx[ik] = jmk; // keep track of its minus
    ik++;
    if (ik >= npw)
      break;
  }
#endif
  auto sposet = std::make_unique<FreeOrbital>(kpts);
  sposet->report("  ");
  return sposet;
}

bool FreeOrbitalBuilder::in_list(const int j, const std::vector<int> l)
{
  for (int i = 0; i < l.size(); i++)
  {
    if (j == l[i])
      return true;
  }
  return false;
}

} // namespace qmcplusplus
