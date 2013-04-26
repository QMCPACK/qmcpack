#ifndef QMCPLUSPLUS_RMCSINGLE_OMP_H
#define QMCPLUSPLUS_RMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "Particle/Reptile.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a RMC using threaded execution.
 */
class RMCSingleOMP: public QMCDriver, public CloneManager
{
public:
  /// Constructor.
  typedef PtclAttribTraits::ParticlePos_t ParticlePos_t;
  RMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
               HamiltonianPool& hpool, WaveFunctionPool& ppool);
  bool run();
  bool put(xmlNodePtr cur);
  //inline vector<RandomGenerator_t*>& getRng() { return Rng;}
private:
  ///period for walker dump
  int myPeriod4WalkerDump;
  ///option to enable/disable drift equation for RMC
  string rescaleDrift;
  ///number of beads on the reptile, beta/tau
  int beads;
  ///rescale for time step studies. some int>2 and new beads are inserted in between the old ones.
  int resizeReptile;

  ///projection time of reptile
  RealType beta;
//       vector of indices for the action and transprob
  std::vector<int> Action;
  std::vector<int> TransProb;

  ///check the run-time environments
  void resetRun();
  ///copy constructor
  RMCSingleOMP(const RMCSingleOMP& a): QMCDriver(a),CloneManager(a) { }
  /// Copy operator (disabled).
  RMCSingleOMP& operator=(const RMCSingleOMP&)
  {
    return *this;
  }

};
}

#endif
