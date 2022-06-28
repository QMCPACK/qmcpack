#ifndef QMCPLUSPLUS_FREE_PARTICLE
#define QMCPLUSPLUS_FREE_PARTICLE

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

class FreeParticle : public SPOSet
{
public:
  FreeParticle(const std::vector<PosType>& kpts_cart);
  ~FreeParticle();

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last)
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix& phi,
    GradMatrix& dphi,
    ValueMatrix& d2phi) override;

  // plug r_i into all orbitals
  void evaluateVGL(
    const ParticleSet& P,
    int i,
    ValueVector& pvec,
    GradVector& dpvec,
    ValueVector& d2pvec
  ) override;
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec) override;

  // hessian matrix is needed by backflow
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix& phi,
    GradMatrix& dphi,
    HessMatrix& d2phi_mat) override;

  // derivative of hessian is needed to optimize backflow
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix& phi,
    GradMatrix& dphi,
    HessMatrix& d2phi_mat,
    GGGMatrix& d3phi_mat) override;

  // ---- begin required overrides
  std::unique_ptr<SPOSet> makeClone() const {return std::make_unique<FreeParticle>(*this);}
  void resetParameters(const opt_variables_type& optVariables) override {} //called by BFTrans}
  void setOrbitalSetSize(int norbs) override {APP_ABORT("not implemented")};
  // required overrides end ----
  void report(const std::string& pad) const override;
private:
  const std::vector<PosType> K; // K vectors
  const int mink; // minimum k index
  const int maxk; // maximum number of K vectors
  std::vector<RealType> mK2; // minus K^2
};

} // qmcplusplus
#endif
