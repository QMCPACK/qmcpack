#include "Particle/LongRange/EwaldHandler2D.h"

namespace qmcplusplus
{

EwaldHandler2D::EwaldHandler2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in)
{
  if (ref.getLattice().ndim != 2)
    APP_ABORT("2D Ewald requires 2D Lattice");
  LR_rc = ref.getLattice().LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.getLattice().LR_kc; // get_kc() is used in QMCFiniteSize
  alpha = std::sqrt(LR_kc/2.0/LR_rc);
  area = ref.getLattice().Volume/ref.getLattice().R(2,2);
  // report
  app_log() << "    alpha = " << alpha << " area = " << area << std::endl;
  fillFk(ref.getSimulationCell().getKLists());
}

EwaldHandler2D::EwaldHandler2D(const EwaldHandler2D& aLR)
  : LRHandlerBase(aLR),
    alpha(aLR.alpha),
    area(aLR.area)
{
}

void EwaldHandler2D::fillFk(const KContainer& KList)
{
  const mRealType knorm = 2.0*M_PI / area;
  const mRealType kalpha = 1.0 / (2.0*alpha);
  mRealType kmag, uk;

  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);

  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    kmag = std::sqrt(KList.ksq[ki]);
    uk = knorm * erfc(kalpha*kmag)/kmag;
    Fk_symm[ks] = uk;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = uk;
  }
}

} // qmcplusplus
