#include "QMCHamiltonians/SpaceWarpTransformation.h"
#include "Particle/DistanceTable.h"
#include "type_traits/scalar_traits.h"
namespace qmcplusplus
{
SpaceWarpTransformation::SpaceWarpTransformation(ParticleSet& elns, const ParticleSet& ions)
    : myTableIndex(elns.addTable(ions)), Nelec(elns.getTotalNum()), Nions(ions.getTotalNum()), swpow(4.0)
{
  warpval.resize(Nelec, Nions);
  gradval.resize(Nelec, Nions);
}

SpaceWarpTransformation::RealType SpaceWarpTransformation::f(RealType r) { return std::pow(r, -swpow); }

SpaceWarpTransformation::RealType SpaceWarpTransformation::df(RealType r) { return -swpow * std::pow(r, -(swpow + 1)); }
//Space warp functions have the form w_I(r_i) = F(|r_i-R_I)/Sum_J F(|r_i-R_J|).  Hence the intermediate we will
//precompute is the matrix "warpval[i][J] = F(|r_i-R_J|) and gradval[i][J]=Grad(F(|r_i-R_J|)).
//This allows the calculation of any space warp value or gradient by a matrix lookup, combined with a sum over columns.
void SpaceWarpTransformation::computeSWTIntermediates(ParticleSet& P, const ParticleSet& ions)
{
  const auto& d_ab(P.getDistTableAB(myTableIndex));
  for (size_t iel = 0; iel < Nelec; ++iel)
  {
    const auto& dist = d_ab.getDistRow(iel);
    const auto& dr   = d_ab.getDisplRow(iel);
    for (size_t ionid = 0; ionid < Nions; ++ionid)
    {
      warpval[iel][ionid] = f(dist[ionid]);
      gradval[iel][ionid] = -dr[ionid] *
          (df(dist[ionid]) / dist[ionid]); //because there's a -1 in distance table displacement definition.  R-r :(.
    }
  }
}

//This function handles parsing of the intermediate matrices and construction of the w_I(r_i) and Grad_i(w_I(r_i)) functions
// that appear in the space warp transformation formulas.
void SpaceWarpTransformation::getSWT(int iat, ParticleScalar_t& w, Force_t& grad_w)
{
  for (size_t iel = 0; iel < Nelec; iel++)
  {
    RealType warpdenom = 0.0;
    PosType denomgrad  = 0.0;
    for (size_t ionid = 0; ionid < Nions; ionid++)
    {
      warpdenom += warpval[iel][ionid];
      denomgrad += gradval[iel][ionid];
    }
    w[iel]      = warpval[iel][iat] / warpdenom;
    grad_w[iel] = gradval[iel][iat] / warpdenom - w[iel] * denomgrad / warpdenom;
  }
}
//This function returns Sum_i w_I(r_i) Grad_i(E_L) (as el_contribution)  and Sum_i[ w_I(r_i)Grad_i(logpsi)+0.5*Grad_i(w_I(r_i)) (as psi_contribution).  See Eq (15) and (16) respectively.
void SpaceWarpTransformation::computeSWT(ParticleSet& P,
                                         const ParticleSet& ions,
                                         Force_t& dEl,
                                         ParticleGradient_t& dlogpsi,
                                         Force_t& el_contribution,
                                         Force_t& psi_contribution)
{
  el_contribution  = 0;
  psi_contribution = 0;
  ParticleScalar_t w;
  Force_t gradw;
  w.resize(Nelec);
  gradw.resize(Nelec);

  PosType gwfn = 0;
  computeSWTIntermediates(P, ions);

  for (size_t iat = 0; iat < Nions; iat++)
  {
    w     = 0;
    gradw = 0;
    getSWT(iat, w, gradw);
    for (size_t iel = 0; iel < Nelec; iel++)
    {
      el_contribution[iat] += w[iel] * dEl[iel];

#if defined(QMC_COMPLEX)
      convert(dlogpsi[iel], gwfn);
#else
      gwfn = dlogpsi[iel];
#endif
      psi_contribution[iat] += w[iel] * gwfn + 0.5 * gradw[iel];
    }
  }
  //REMOVE ME
  //app_log()<<"FINAL SPACE WARP TRANSFORM\n";
  //app_log()<<" Using dEl="<<dEl<<std::endl;
  //app_log()<<" Using dlogpsi="<<dlogpsi<<std::endl;

  //app_log()<<"Space warp EL piece\n";
  //app_log()<<" EL = "<<el_contribution<<std::endl;
  //app_log()<<"Space warp logpsi piece\n";
  //app_log()<<" psipiece = "<<psi_contribution<<std::endl;
}

} // namespace qmcplusplus
