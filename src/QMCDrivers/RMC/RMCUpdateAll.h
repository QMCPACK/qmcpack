#ifndef QMCPLUSPLUS_RMC_UPDATEALL_H
#define QMCPLUSPLUS_RMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 *@brief Implements the RMC algorithm using all electron moves
 */
  class RMCUpdateAllWithDrift:public QMCUpdateBase
  {
  public:
    /// Constructor.

    enum
    { SYM_ACTION, DMC_ACTION };

      RMCUpdateAllWithDrift (MCWalkerConfiguration & w,
			     TrialWaveFunction & psi, QMCHamiltonian & h,
			     RandomGenerator_t & rg, std::vector < int >act,
			     std::vector < int >tp);
     ~RMCUpdateAllWithDrift ();
    void advanceWalkers (WalkerIter_t it, WalkerIter_t it_end, bool measure);
    void advanceWalkersVMC ();
    void advanceWalkersRMC ();
    void checkReptile (WalkerIter_t it, WalkerIter_t it_end);
    void initWalkers (WalkerIter_t it, WalkerIter_t it_end);

    void accumulate (WalkerIter_t it, WalkerIter_t it_end);

    bool put (xmlNodePtr cur);

  private:
    /// Copy Constructor (disabled)
      RMCUpdateAllWithDrift (const RMCUpdateAllWithDrift &
			     a):QMCUpdateBase (a), Action (a.Action),
      TransProb (a.TransProb)
    {
    }
    /// Copy operator (disabled).
    RMCUpdateAllWithDrift & operator= (const RMCUpdateAllWithDrift &)
    {
      return *this;
    }
    std::vector < int >Action, TransProb;

    bool scaleDrift;
    IndexType actionType;

    IndexType vmcSteps;
    IndexType equilSteps;
    IndexType vmcToDoSteps;
    IndexType equilToDoSteps;

  };


}

#endif
